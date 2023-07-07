#!/bin/zsh

## Version 1
## This script sets up basically EVERYTHING needed for running 
## FD+IDR CAMPARI simulations from a list of sequences on the HTCF cluster.
## It serves both as a functioning script but also as living documentation
## of how to do this. Updated for Holehouselab as of March 5th 2020

## !!!! THESE ARE THE ONLY THINGS YOU SHOULD NEED TO CHANGE !!!!

# Path to FLAMINGO repository
FLAMINGO_DIR="/home/degriffith/repos/flamingo"

# define the input sequence files
idr_filename="IDR_variants.isf"

# define the input structures file (list of PDBs)
### CHANGE THESE and RERUN ###
pdbstructure_filename="pdb_structures.txt"
fd_filename="FD.isf"

# define the fixed IDR residues (must be present in PDB)
fixed_idr_residues="fixed_residues.txt"

# Simulation parameters
MC_MODE=ts          # MUST be either 'ts' (Temperature Sweep), 'hs' (Hamiltonian Switch), 'ev', or 'standard'
TEMPERATURE=360     # 360 Kelvin
SALT=0.015          # 15 mM NaCl
PRE_EQ=2000000      #  2 M - number of steps before main simulation starts
EQ=4000000          #  4 M - number of equilibration steps in main simulation
PROD=60000000       # 60 M - number of steps in main simulation AFTER equilibration
REPS=4              # number of independent replicas (per mode)
XTCOUT=20000        # frequency at which .xtc file is written to (20K is a good number)
SIM_MODE=coil       # MUST be either 'combined', 'helical', or 'coil'
CAMPARI_VERSION=3   # version of CAMPARI to use (2 or 3)
PRIORITY=normal     # Set job priorities

# Parameters specific for TemperatureSweep or HamiltonianSwitch
AUXCHAIN_ENTERPROB=0.8     # Acceptance probability for entering the auxiliary chain (0.8 is good)
AUXCHAIN_ENTERFREQ=25000   # Number of steps between auxiliary chain attempts
AUXCHAIN_NSTEPS=500        # Number of steps for each subchain in the auxiliary chain (10 subchains for TSMC, 1 for HSMC)

# ------------------------------------------------------------------------
# Probably no need to change anything below this line!
# ------------------------------------------------------------------------

# VALIDATE PRIORITY
if [ "$PRIORITY" = "superlow" ]
then
    p=$1
elif [ "$PRIORITY" = "low" ]
then
    p=$1
elif [ "$PRIORITY" = "normal" ]
then
    p=$1
elif [ "$PRIORITY" = "high" ]
then
    p=$1
else
    echo "Invalid priority passed"
    exit
fi

# Validate MC_MODE
if [ "$MC_MODE" = "ts" ]
then
    echo "Temperature Sweep mode"
elif [ "$MC_MODE" = "hs" ] 
then
    echo "Hamiltonian Switch mode"
elif [ "$MC_MODE" = "ev" ]
then
    echo "EV simulation"
elif [ "$MC_MODE" = "standard" ]
then
    echo "Standard MC mode"
else
    echo "Invalid MC mode"
    exit
fi

# submission_list.txt is a convenient file we generate that lists all the directories
# generated, which means when we eventially do submit jobs we can use this as a reference
# we delete it here if it exists as we will then append to it in the while loop
if [ -f submission_list.txt ]
then    
    rm submission_list.txt
fi

# if a launch_all.sh file exists delete it
if [ -f launch_all.sh ]
then    
    rm launch_all.sh
fi

# Read in FD
fd_name=$(cat $fd_filename | awk {'print $1'})
fd_sequence=$(cat $fd_filename | awk {'print $2'})

# this while loop reads in the file (defined at the END) and basically reads EACH
# line one-by-one into the variable $line
while read -r line
do

    # get the sequence name/index from the first column
    idr_name=$(echo "$line" | awk {'print $1'})
    
    # get the sequence from the second column
    idr_sequence=$(echo "$line" | awk {'print $2'})

    # get fixed residues for this idr
    idr_fixed_line=$(grep $idr_name $fixed_idr_residues)
    idr_fixed=$(echo "$idr_fixed_line" | awk '{print $2}')

    # Check flexible IDR region for prolines
    python "${FLAMINGO_DIR}/setup_scripts/check_prolines.py" $idr_sequence $idr_fixed
    if [ "$?" -ne "0" ]
    then
        IDR_prolines=FALSE
    else
        IDR_prolines=TRUE
    fi  

    # get PDB for this idr
    pdbstructure=$(grep $idr_name $pdbstructure_filename | awk {'print $2'})   

    # Clean up HIS->HIE residues in PDB and rename as start.pdb for later steps
    cp $pdbstructure pdbstructure.tmp
    sed -i "s/HIS/HIE/g" pdbstructure.tmp
    mv pdbstructure.tmp start.pdb

    name="${idr_name}_${fd_name}"

    # if the directory called ${name} does not exist
    if [ ! -d "${name}" ]
    then
        # make a directory (defined by the name) and move
        # into that directory
        mkdir $name
        cd $name

        # Temporary files with just single line for this IDR
        echo $line > idr.tmp
        echo $idr_fixed_line > fixed.tmp
        cp ../start.pdb .

        # get keyfiles - the keyfiles linked here are well-defined keyfile for
        # running cABSINTH simulations and can be used without modificiation
        # (autoSim will modify it in all the ways it needs to be modified)

        if [ "${CAMPARI_VERSION}" = "3" ]
        then
            cp "${FLAMINGO_DIR}/keyfiles/build_FD_IDR.key" build.key
        if [ "$MC_MODE" = "ev" ]
            then
                cp "${FLAMINGO_DIR}/keyfiles/run_EV_FD_IDR.key" run.key
            else
                cp "${FLAMINGO_DIR}/keyfiles/run_FD_IDR.key" run.key
        fi
        else
            echo "Invalid option passed for CAMPARI_VERSION: ${CAMPARI_VERSION} (must be 3)"
            exit 1
        fi
        
        ##################################################################
        # New stuff:

        # If no prolines, modify the build and key file FMCSC_PKRFREQ -> 0
        if [ "$IDR_prolines" = "FALSE" ]
        then
            echo "IDR does not contain proline. Changing FMCSC_PKRFREQ -> 0 in key files."
            sed -i "s/FMCSC_PKRFREQ 0.1/FMCSC_PKRFREQ 0/g" build.key
            sed -i "s/FMCSC_PKRFREQ 0.1/FMCSC_PKRFREQ 0/g" run.key
        else
            echo "IDR contains proline. No modification to key files."
        fi


        # Create PSWFILE.psw (very creative name, I know)
        python "${FLAMINGO_DIR}/setup_scripts/create_psw_file.py" $fd_sequence fixed.tmp PSWFILE.psw --idr-first --idr-caps

        # Create seq.in
        python "${FLAMINGO_DIR}/setup_scripts/create_sequence_file.py" $fd_sequence idr.tmp seq.in --idr-first --idr-caps

        # Run 1-step simulation to generate complete PDB structure
        campari3 -k build.key > build.log
        rm __END.pdb
        rm *.int

        # Create dres.in
        python "${FLAMINGO_DIR}/setup_scripts/create_restraint_file.py" __START.pdb fixed.tmp dres.in --idr-first --force-constant 500.0

        # Copy modified autoSim
        cp "${FLAMINGO_DIR}/setup_scripts/autoSim_vFD_IDR.sh" .

        ##############
        ## Now we run autoSim! autoSim takes an input sequence ($sequence), keyfile,
        ## and a bunch of paratemers and constructs the file system and scripts needed
        ## to run all the CAMPARI simulations. In later versions it will also automatically do
        ## analysis but this is not (yet) implemented
        
        # The keywords used are explained below
        # -k  - keyfile
        # -m  - simulation mode (combined means we run simulations that start from a hexlix and
        #       simulations that start from a coil - if we've converged these should end in the
        #       same place.
        # -f  - steps during pre-equilibraion [as helix or as coil]
        # -r  - number of replicas (note for combined this means 2x because we run $n replicas from
        #       a helix and $n from a coil)
        # -e  - steps for true equilibration
        # -p  - production steps
        # -x  - frequency that coordinates are written
        # -t  - temperature (in kelvin)
        # -s  - salt concentration in Molar (so 15 mM)
        # -h  - temperature for pre-equilibration simulation (default is -t + 50)
        
        zsh autoSim_vFD_IDR.sh -i ${idr_sequence} -k run.key -f ${PRE_EQ} -r ${REPS} -e ${EQ} -p ${PROD} -x ${XTCOUT} -t ${TEMPERATURE} -s ${SALT} -m ${SIM_MODE} -v ${CAMPARI_VERSION}

        # the $? variable returns the exit status of autoSim
        if [ "$?" -ne "0" ]
        then
            echo "autoSim has failed"
            exit 1
        fi  

        ##### NEW STUFF HERE - TempSweep/HamSwitch #####

        # Edit temperature in run.key, since we will not be using the keyfiles created by autoSim
        sed -i "s/.*FMCSC_TEMP.*/  FMCSC_TEMP ${TEMPERATURE}/g" run.key

        # Iterate through sub directories (coil/helical_start, Nreps)
        if [ "$SIM_MODE" = "coil" ] || [ "$SIM_MODE" = "combined" ]
        then
            rm pre_eq.key
            rm production.key

            cd coil_start
            rm production.key
            rm pre_eq.key

            for i in {1..$REPS}; do
                cd $i

                # cp ../../dres.in
                rm campari_bash.sh
                rm pre_eq.key
                rm production.key
                cp ../../run.key .
                python "${FLAMINGO_DIR}/setup_scripts/create_tsmc_run_script.py" --out run_sims.py --keyfile run.key --preeq $PRE_EQ --eq $EQ --mc $PROD --mode $MC_MODE --temp $TEMPERATURE --aux-enter-prob $AUXCHAIN_ENTERPROB --aux-chain-freq $AUXCHAIN_ENTERFREQ --aux-chain-steps $AUXCHAIN_NSTEPS   
                cd ..
            done
            cd ..
        fi

        if [ "$SIM_MODE" = "helical" ] || [ "$SIM_MODE" = "combined" ]
        then
            rm pre_eq_helix.key
            rm production.key

            cd helical_start
            rm production.key
            rm pre_eq_helix.key
            
            for i in {1..$REPS}; do
                cd $i
                rm campari_bash.sh
                rm pre_eq_helix.key
                rm production.key
                cp ../../run.key .
                python "${FLAMINGO_DIR}/setup_scripts/create_tsmc_run_script.py" --out run_sims.py --keyfile run.key --preeq $PRE_EQ --preeq-helix --eq $EQ --mc $PROD --mode $MC_MODE --temp $TEMPERATURE --aux-enter-prob $AUXCHAIN_ENTERPROB --aux-chain-freq $AUXCHAIN_ENTERFREQ --aux-chain-steps $AUXCHAIN_NSTEPS 
                cd ..
            done
            cd ..
        fi
            
        rm run_seq.sh

        # copy the script run_seq.sh, which could be improved, but automates job submission
        cp "${FLAMINGO_DIR}/setup_scripts/run_seq_tsmc.sh" .
        
        echo "${name}" > JOB_PREFIX.txt 

        rm *.tmp
        rm tmp.*

        cd ..
        rm start.pdb

    else
        echo "Directory ${name} already exists"
    fi
    
    echo "$name" >> submission_list.txt
    echo "cd ${name}; zsh run_seq_tsmc.sh ${PRIORITY}; cd .." >> launch_all.sh   
    ###############

done < "$idr_filename"
