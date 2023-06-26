#!/usr/bin/env python

import argparse
import numpy as np
import mdtraj as md
import pandas as pd
from collections import Counter
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt

# Parse command line args
parser = argparse.ArgumentParser(description='Extract binding motif residues from AF2 structure.')
parser.add_argument('--pdb', type=str, required=True, 
                    help='input AF2 pdb')
parser.add_argument('--name', type=str, required=True, 
                    help='name of structure (used in output file)')
parser.add_argument('--out', type=str, required=True, 
                    help='output file for motif residues')
parser.add_argument('--idr-first', action='store_true', 
                    help='flag for if the IDR is the first chain in the input PDB')
parser.add_argument('--distance-cutoff', type=float, default=7.0,
                    help='distance in Angstroms used to identify the motif (default=7.0)')
parser.add_argument('--alphafold-cutoff', type=float, default=50,
                    help='pLDDT threshold required to define a motif (default=50)')
parser.add_argument('--silent', '-s', action='store_true')

args = parser.parse_args()

# Load in PDB file to mdtraj
traj = md.load(args.pdb)
top = traj.topology
df = top.to_dataframe()[0]

if args.idr_first:
    idr_idx = 0
    fd_idx = 1
else:
    idr_idx = 1
    fd_idx = 0

# Get CA atoms for IDR and FD chains
idr_seq = top.to_fasta()[idr_idx]
idr_df = df[df['chainID'] == idr_idx]
fd_df = df[df['chainID'] == fd_idx]

idr_CA_idxs = idr_df[idr_df['name'] == 'CA']['serial'].values
fd_CA_idxs = fd_df[fd_df['name'] == 'CA']['serial'].values

pairs = np.array([[i, j] for i in fd_CA_idxs for j in idr_CA_idxs])

# Compute distances (x10 to convert to Angstroms)
distances = md.compute_distances(traj, pairs)[0] * 10

# Find pairs of residues less than cutoff
contact_pairs = pairs[distances <= args.distance_cutoff]

idr_contacts_atoms = Counter([top.atom(p[1]) for p in contact_pairs])

# Get pLDDT scores of AF2 prediction from PDB
# TODO: this won't work if fields in PDB are not separated by whitespace!
with open(args.pdb) as f:
    lines = [x.strip().split() for x in f]
atom_lines = [line for line in lines if line[0] == 'ATOM']
chain_ids = sorted(list(set([line[4] for line in atom_lines])))
idr_chain_id = chain_ids[idr_idx]
idr_plddts = np.array([float(line[-2]) for line in atom_lines if 
        ((line[4] == idr_chain_id) and (line[2] == 'CA'))])
    
confident_binding_region = np.where(idr_plddts > args.alphafold_cutoff)[0]

# Get IDR seq and make density curve over counts
resid_counts = {}
for res in idr_contacts_atoms:
    resid_counts[res.residue.resSeq-1] = idr_contacts_atoms[res]

count_str = ''
plddt_str = ''
for i in range(len(idr_seq)):
    if i in resid_counts:
        count_str += str(resid_counts[i])
    else:
        count_str += '0'

    if i in confident_binding_region:
        plddt_str += '*'
    else:
        plddt_str += ' '

fixed_residues = []
for n_contacts, plddt, aa in zip(count_str, plddt_str, idr_seq):
    if plddt == '*' and int(n_contacts) > 0 and aa not in 'DERKG':
        fixed_residues.append('1')
    else:
        fixed_residues.append('0')

try:
    first_fixed = fixed_residues.index('1')
    last_fixed = len(fixed_residues) - 1 - fixed_residues[::-1].index('1')
    fixed_residues = '0'*first_fixed + '1'*(last_fixed-first_fixed+1) + '0'*(len(idr_seq)-last_fixed-1)
except:
    # No motif identified
    fixed_residues = '0'*len(idr_seq)

if not args.silent:
    print()
    print(args.name + ':')
    print('------------')
    print('Sequence:  | ', idr_seq)
    print('Contacts:  | ', count_str)
    print('AF2 pLDDT: | ', plddt_str)
    print('Motif:     | ', fixed_residues)

with open(args.out, 'w') as outfile:
    outfile.write(f'{args.name} {fixed_residues}\n')

