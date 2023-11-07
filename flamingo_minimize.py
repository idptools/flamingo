#!/usr/bin/env python

import mdtraj as md
import os
import argparse

def parse_key_file(file_path):
    key_value_dict = {}
    
    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line:
                    keyword, value = line.split()
                    if not isinstance(keyword, str):
                        raise TypeError(f'Keywords must be strings: {keyword}')

                    key_value_dict[keyword.upper()] = value

    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

    return key_value_dict

def process_keywords(keyword_dict):
    required_keywords = ['PDB_FILE', 'FD_CHAIN_ID', 'MOTIF']
    optional_keywords = {'GROMACS_MINIMIZE':True, 'CAMPARI_MINIMIZE':True, 'CAMPARI_MD':True,
                         'GROMACS_NSTEPS':50000, 'GROMACS_EMTOL':1000.0, 
                         'N_MOTIF_PADDING_RESIDUES':2,
                         'CAMPARI_MIN_NSTEPS':1000000, 
                         'CAMPARI_MD_NSTEPS':1000000, 'CAMPARI_MD_TEMP':100,
                         'MOTIF_START_IDX':None, 'MOTIF_END_IDX':None}

    # Validate the required keywords
    for req_keyword in required_keywords:
        if req_keyword not in keyword_dict:
            raise Exception(f'Missing required keyword: {req_keyword}')

    if isinstance(keyword_dict['FD_CHAIN_ID'], str) and len(keyword_dict['FD_CHAIN_ID']) == 0:
        fd_chain_id = 'AB'.find(keyword_dict['FD_CHAIN_ID'].upper())
    elif isinstance(keyword_dict['FD_CHAIN_ID'], int) and keyword_dict['FD_CHAIN_ID'] in [0,1]:
        fd_chain_id = keyword_dict['FD_CHAIN_ID']
    else:
        raise Exception(f'Invalid chain ID: "{keyword_dict['FD_CHAIN_ID']}"')

    pdb = md.load(keyword_dict['PDB_FILE'])
    topology = pdb.top
    sequence_list = topology.to_fasta()

    if len(sequence_list) != 2:
        raise Exception('PDB file must have two chains.')

    idr_seq = sequence_list[1-chain_id]
    fd_seq = sequence_list[chain_id]

    # Find motif in IDR sequence
    motif_occurences = idr_seq.count(keyword_dict['MOTIF'])
    if motif_occurences == 0:
        raise Exception(f'Motif "{keyword_dict["MOTIF"]}" not found in IDR sequence : {idr_seq}')
    elif motif_occurences == 1:
        motif_start_idx = idr_seq.find(keyword_dict['MOTIF']) + 1
        motif_end_idx = motif_start_idx + len(keyword_dict['MOTIF']) - 1
    else: # More than one occurence
        # Require that MOTIF_START/END_IDX keywords were provided
        if keyword_dict['MOTIF_START_IDX'] is None or keyword_dict['MOTIF_END_IDX'] is None:
            raise Exception(f'Must provide MOTIF_START_IDX and MOTIF_END_IDX\n since there are multiple occurences of motif: "{keyword_dict["MOTIF"]}"')
        else:
            motif_start_idx = keyword_dict['MOTIF_START_IDX']
            motif_end_idx = keyword_dict['MOTIF_END_IDX']

            if idr_seq[motif_start_idx-1:motif_end_idx] != keyword_dict['MOTIF']:
                raise Exception(f'Provided motif idxs: {idr_seq[motif_start_idx-1:motif_end_idx]} does not match motif: {keyword_dict['MOTIF']}')

    # Create and return new dictionary with all of the keywords + defaults
    processed_dict = {}
    # Optional values as defaults first
    for key, value in optional_keywords.items():
        processed_dict[key] = value
    # User provided values to overwrite defaults
    for key, value in keyword_dict.items():
        processed_dict[key] = value

    processed_dict['MOTIF_START_IDX'] = motif_start_idx
    processed_dict['MOTIF_END_IDX'] = motif_end_idx
    processed_dict['IDR_SEQUENCE'] = idr_seq
    processed_dict['FD_SEQUENCE'] = fd_seq
    processed_dict['FD_CHAIN_ID'] = fd_chain_id
    processed_dict['IDR_CHAIN_ID'] = 1 - fd_chain_id
    return processed_dict

def generate_gromacs_mdp_file(fname, nsteps, emtol, coulomb_type):
    outstr=  '; Parameters describing what to do, when to stop and what to save\n'
    outstr+= 'integrator  = steep         ; Algorithm (steep = steepest descent minimization)\n'
    outstr+=f'emtol       = {emtol}       ; Stop minimization when the maximum force < {emtol} kJ/mol/nm\n'
    outstr+= 'emstep      = 0.01          ; Minimization step size\n'
    outstr+=f'nsteps      = {nsteps}         ; Maximum number of (minimization) steps to perform\n\n'
    outstr+= '; Parameters describing how to find the neighbors of each atom and how to calculate the interactions\n'
    outstr+= 'nstlist         = 1         ; Frequency to update the neighbor list and long range forces\n'
    outstr+= 'cutoff-scheme   = Verlet    ; Buffered neighbor searching\n'
    outstr+= 'ns_type         = grid      ; Method to determine neighbor list (simple, grid)\n'
    outstr+=f'coulombtype     = {coulomb_type}      ; Treatment of long range electrostatic interactions\n'
    outstr+= 'rcoulomb        = 1.0       ; Short-range electrostatic cut-off\n'
    outstr+= 'rvdw            = 1.0       ; Short-range Van der Waals cut-off\n'
    outstr+= 'pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions\n'

    with open(fname, 'w') as out:
        out.write(outstr)


# Main execution
if __name__ == "__main__":
    # Process command line arguments
    parser = argparse.ArgumentParser(description="Minimize a bound IDR:FD structure in preparation of a MCMC simulation.")
    parser.add_argument("-k", "--keyfile", required=True, help="Path to the key file")
    # TODO: add logging functionality

    args = parser.parse_args()
    keyfile_path = args.keyfile

    # Parse and process keyfile
    keywords = parse_key_file(keyfile_path)
    processed_keywords = process_keywords(keywords)

    # Range of residues for IDR truncation
    padding_start = max(processed_keywords['MOTIF_START_IDX']-processed_keywords['N_MOTIF_PADDING_RESIDUES'], 1)
    padding_end = min(processed_keywords['MOTIF_END_IDX']+processed_keywords['N_MOTIF_PADDING_RESIDUES'], 
                        len(processed_keywords['IDR_SEQUENCE']))

    # Build minimization script(s)
    minimization_wrapper_script = ['#!/bin/bash', '', 
                                   'mkdir /minimize',
                                   'cd /minimize',
                                   f'cp /mnt/{processed_keywords["PDB_FILE"]} /minimize/start.pdb'

    if processed_keywords['GROMACS_MINIMIZE']:

        # Create GROMACS minimization script & organize file system
        gromacs_minim_script = ['#!/bin/bash', '',
                                 'mkdir /minimize/gromacs',
                                'cd /minimize/gromacs',
                                'mv /minimize/start.pdb /minimize/gromacs/start.pdb', '']

        # Create necessary .mdp files
        generate_gromacs_mdp_file('ions.mdp', processed_keywords['GROMACS_NSTEPS'],
                                  processed_keywords['GROMACS_EMTOL'], 'cutoff')
        generate_gromacs_mdp_file('minim.mdp', processed_keywords['GROMACS_NSTEPS'],
                                  processed_keywords['GROMACS_EMTOL'], 'PME')

        gromacs_minim_script.extend(['# Prepare the system',
            'gmx pdb2gmx -f start.pdb -o complex.gro -ff oplsaa -water tip3p -ignh # convert to gro file and pick ff',
            'gmx editconf -f complex.gro -o complex_newbox.gro -c -d 1.0 -bt cubic # prepare to solvate',
            'gmx solvate -cp complex_newbox.gro -cs spc216.gro -o complex_solv.gro -p topol.top # solvate',
            'gmx grompp -f ions.mdp -c complex_solv.gro -p topol.top -o ions.tpr',
            'echo 13 | gmx genion -s ions.tpr -o complex_solv_ions.gro -p topol.top -pname NA -nname CL -neutral # Add ions',
            '','# Minimize',
            'gmx grompp -f minim.mdp -c complex_solv_ions.gro -p topol.top -o em.tpr',
            'gmx mdrun -v -deffnm em # minimize',
            '','# Post processing',
            'echo 1 1 | gmx trjconv -f em.gro -s em.tpr -o gromacs_minimized.pdb -pbc nojump -center # convert to pdb'])

        # Truncate IDR to the motif + padding residues
        if processed_keywords['IDR_CHAIN_ID']:
            gromacs_minim_script.append(f'/flamingo/setup_scripts/modify_complex_pdb.py --pdb /minimize/gromacs/gromacs_minimized.pdb --out /minimize/truncated.pdb --clip-idx {padding_start} {padding_end} --idr-first')
        else:
            gromacs_minim_script.append(f'/flamingo/setup_scripts/modify_complex_pdb.py --pdb /minimize/gromacs/gromacs_minimized.pdb --out /minimize/truncated.pdb --clip-idx {padding_start} {padding_end}')

    else:
        # Truncate IDR to the motif + padding residues
        if processed_keywords['IDR_CHAIN_ID']:
            gromacs_minim_script.append(f'/flamingo/setup_scripts/modify_complex_pdb.py --pdb /minimize/start.pdb --out /minimize/truncated.pdb --clip-idx {padding_start} {padding_end} --idr-first')
        else:
            gromacs_minim_script.append(f'/flamingo/setup_scripts/modify_complex_pdb.py --pdb /minimize/start.pdb --out /minimize/truncated.pdb --clip-idx {padding_start} {padding_end}')


    if processed_keywords['CAMPARI_MINIMIZE']:
        # Do something
    else:
        # Rename input PDB to match CAMPARI output pdb name

    if processed_keywords['CAMPARI_MD']:
        # Do something
    else:
        # Rename input PDB to match CAMPARI output pdb name


    '''
    Final script:
    #!/bin/bash

    # move files around


    # Minimize
    ./run_gromacs_min.sh
    ./run_campari_min.sh
    ./run_campari_md

    # move files around

    '''

    # Run minimization script
