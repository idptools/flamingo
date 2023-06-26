#!/bin/zsh

##  CHANGELOG
#
# Version 1.0 - initial
#
# Version 1.1
# Added the -h flag for setting the 'heating' temperature in pre-equilibration
#
# Version 1.2
#
# Version 1.3 Added CAMPARI version option
#
# Version 1.4 Added 
#
# Version 1.5 Changed default xtcout to 40000 (from 10000)
#             
# Version 1.6 Added PTM support!
#     
# Version FD_IDR Modified to take in IDRs bound to FD structures        

# SOURCE GLOBAL VARIABLES
source /work/holehouselab/configurations/misc/autoSim_config.sh



########################################################################################
validate_writability(){
    # no input - just checks if we can write in the pwd

    touch tmp 2>/dev/null
    if [ "$?" -ne  "0" ]
    then
	echo "ERROR: Do not have write permission for current directory"
	echo "Please run autoSim from a directory where write permission is enabled"
	rm tmp
	exit 1
    fi
    rm tmp
}
########################################################################################
    
    
########################################################################################
replace_line(){

    # $1 - replace string target
    # $2 - replace string replacement
    # $3 - file upon which replacement should happen

    sed -i.bak "s/.*${1}.*/${2}/" ${3}
}
########################################################################################

########################################################################################
validate_numerical(){

    # $1 - value being tested
    # $2 - label to be printed on error

    # this will fail if float-casting doesn't work...
    tmp=$(python -c "print(float('$1'))" 2>/dev/null)

    if [ "$?" -ne  "0" ]
    then
	echo "ERROR: input value $2 is not numerical [$1]"
	echo "Please ensure this value is numerical and try again"
	exit 1
    fi    
}

########################################################################################


VALID_MODES="coil helix combined"
VERSION=FD_IDR

########################################################################################
##
##
##                            START OF PROGRAM
##
##
########################################################################################

while getopts ":i:m:e:p:t:h:b:s:k:x:r:f:h:v:a:" opt; do
  case $opt in
    i) seq_string="$OPTARG"
    ;;
    m) sim_mode="$OPTARG"
    ;;
    e) eq_steps="$OPTARG"
    ;;
    p) prod_steps="$OPTARG"
    ;;
    t) temperature="$OPTARG"
    ;;
    s) salt="$OPTARG"
    ;;
    k) keyfile="$OPTARG"
    ;;
    x) xtcout="$OPTARG"
    ;;
    r) n_reps="$OPTARG"
    ;;
    f) preeq_steps="$OPTARG"
    ;;
    h) preeq_temp="$OPTARG"
    ;;
    v) campari_version="$OPTARG"
    ;;
    a) ptmfile="$OPTARG"
    ;;
    \?) echo "Invalid option ($OPTARG) discarding..." >&2
    ;;
  esac
done

echo $campari_version
echo "---------------"

echo ""
echo "Welcome to autoSim - a wrapper for building and running CAMPARI simulations"
echo "Version $VERSION"
tmp=$(date)
echo "Launched at $tmp"
echo ""

if [ -z "$seq_string" ]
then
    echo "No sequence string provided. Run simulation takes four possible arguments"
    echo ""
    echo "-i     input amino acid sequence"
    echo "-k     input keyfile [run.key]"    
    echo "-r     number of replicas [5]"    
    echo "-s     NaCl concentration (M) [0.010] "
    echo "-m     mode to run (helix, coil) - reflects pre-equilibration [coil]"
    echo "-f     number of steps for pre-equilibration (i.e. high temp coil or induced helix) [500000] (0.5 M)"
    echo "-h     temperature (K) equilibration simulation will be run at [production + 50]"
    echo "-t     temperature (K) production simulation will be run at [310]"
    echo "-e     number of steps for equilibration in production simulation [2000000] (2 M)"
    echo "-p     total number of steps for production simulation [30000000] (30 M)"
    echo "-x     xtc output frequency [10000]"
    echo "-v     CAMPARI version [2]"
    echo "-a     PTMFile (a[ugmentation]) - see help in simprep for description"
    echo ""
    exit 1    
fi


# check we have write position in the current directory
validate_writability

if [ -d helical_start ]
then
    echo "ERROR: helical_start directory already present"
    echo "Cannot overwrite..."
    exit 1
fi

if [ -d coil_start ]
then
    echo "ERROR: coil_start directory already present"
    echo "Cannot overwrite..."
    exit 1
fi


# set defaults
# 30 M
if [ -z "$prod_steps" ]
then
    prod_steps=30000000
fi


## Mode selection/validation
## 
if [ -z "$sim_mode" ]
then
    sim_mode='coil'
else

    echo ${sim_mode}
    
    # sim_mode is not one of the VALID_MODES
    if [ "${sim_mode}" = "combined" ] || [ "${sim_mode}" = "helix" ] || [ "${sim_mode}" = "coil" ]
    then
	echo ""
    else
	echo "ERROR: Invalid mode provided (must be one of 'coil', 'helix', or 'combined')"
	exit 1
    fi

fi

# 2 M
if [ -z "$eq_steps" ]
then
    eq_steps=2000000
fi

# preeq_steps
if [ -z "$preeq_steps" ]
then
    preeq_steps=500000
fi

# 310 K
if [ -z "$temperature" ]
then
    temperature=310
fi


# 310 K
if [ -z "$preeq_temp" ]
then
    preeq_temp=$(expr ${temperature} + 50)
fi

# 10 mM NaCl
if [ -z "$salt" ]
then
    salt=0.010
fi

if [ -z "$keyfile" ]
then
    keyfile='run.key'
fi

#....................................................
if [ -z "$n_reps" ]
then
    n_reps=5
fi
#....................................................

#....................................................
## start XTC validation/defaults
if [ -z "$xtcout" ]
then
    xtcout=40000
fi

if [ "$xtcout" -gt "$prod_steps" ]
then
    echo "ERROR: Cannot output XTC less frequently than the total number of production steps"
    exit 1
fi

if [ "$xtcout" -lt 1 ]
then 
    echo "ERROR: Cannot output XTC more frequently than once a step"
    exit 1
fi
## end XTC
#....................................................

#....................................................
## validate campari version and the define the CAMPARI
# binary location based on the version
if [ -z "${campari_version}" ]
then
    campari_version=2
fi

# define the CAMPARI_BIN
if [ "${campari_version}" = "2" ]
then
    CAMPARI_BIN=${CAMPARI2_BIN}
elif [ "${campari_version}" = "3" ]
then
    CAMPARI_BIN=${CAMPARI3_BIN}
else
    echo "Error: Invalid CAMPARI selection [${campari_version}] must be 2 or 3"
    exit 1
fi
#....................................................


#....................................................
if [ ! -z "${ptmfile}" ]
then
    if [ ! -e "${ptmfile}" ]
    then
	echo "Error: Not PTMfile (-a flag) found [${ptmfile}]"
	exit 1
    fi
fi
    

# Validate numerical values    
validate_numerical $xtcout "-x "
validate_numerical $n_reps "-n"
validate_numerical $eq_steps "e"
validate_numerical $prod_steps "-p"
validate_numerical $temperature "-t"
validate_numerical $salt "-s"
validate_numerical $preeq_steps "-f"
validate_numerical $preeq_temp "-h"
validate_numerical $campari_version "-v"


# Validate input sequence
################################################################################################
## stage 1 validate input sequence
upper_string=$(python -c "seq='${seq_string}'.upper();print(seq)")
seqvalid_return=$(python -c "VALID='ACDEFGHIKLMNPQRSTVWY'
for i in '${upper_string}':
   if i not in VALID:
      print('X')
")


if [ ! -z "$seqvalid_return" ]
then
    echo "ERROR: Invalid sequence provided."
    echo "Sequence [${upper_string}] must only contain cannonical amino acids"
    echo ""
    exit 1
fi
################################################################################################



# Validate keyfile and input parameters
################################################################################################
if [ ! -f "${keyfile}" ]
then     
    echo "ERROR: Keyfile [$keyfile] does not exist"
    exit 1
fi

# if we get here assumed we're gonna go all the way

# compute numbers of frames
n_frames=$(python -c "print(int(float('${prod_steps}')/float('${xtcout}')))")

# if mode is combined we're going to run both helical and non-helical simulations
if [ "${sim_mode}" = "combined" ]
then
    total_n_frames=$(python -c "print(2*int(float('${n_frames}')*float('${n_reps}')))")
else
    total_n_frames=$(python -c "print(int(float('${n_frames}')*float('${n_reps}')))")
fi

seqlen=$(python -c "print(len('${upper_string}'))")

# build sequence file (seq.in) and extract radius of droplet
# Assume seq.in file already present in current directory
echo "${upper_string}" > tmp.fasta
mv seq.in tmp.seq.in
radius=$(simprep -f tmp.fasta -s ${salt}  --no-log-file | grep "radius --" | awk '{ print $3 }')
rm seq.in
mv tmp.seq.in seq.in


# note this overwrites an existing logfile
herenow=$(pwd)
echo "Logfile genertaed by AutoSim" > log.txt
echo "Working directory: $herenow" | tee -a log.txt
date >> log.txt
echo "" | tee -a log.txt
echo "Running simulation with the following parameters:" | tee -a log.txt
echo "*************************************************" | tee -a log.txt
echo "CAMPARI vers.: ${campari_version}" | tee -a log.txt
echo "Mode         : ${sim_mode}" | tee -a log.txt
echo "input keyfile: ${keyfile}" | tee -a log.txt
echo "" | tee -a log.txt
echo "Sequence     : ${upper_string}" | tee -a log.txt

if [ ! -z "${ptmfile}" ]
then
    echo "ptmfile      : ${ptmfile}" | tee -a log.txt
else
        echo "ptmfile      : None" | tee -a log.txt
fi


echo "Sequence len : ${seqlen} residues" | tee -a log.txt
echo "Salt conc.   : ${salt} M" | tee -a log.txt
echo "Pre-eq. Temp : ${preeq_temp} K" | tee -a log.txt
echo "Pre-eq. Steps: ${preeq_steps} steps" | tee -a log.txt

echo "" | tee -a log.txt

echo "Prod. Temp.  : ${temperature} K" | tee -a log.txt
echo "Equil. steps : ${eq_steps} steps" | tee -a log.txt
echo "Prod.  steps : ${prod_steps} steps" | tee -a log.txt
echo "XTC out freq : Every ${xtcout} steps " | tee -a log.txt
echo "Num. frames  : ${n_frames}" | tee -a log.txt
echo "TOTAL frames : ${total_n_frames}" | tee -a log.txt
echo "Droplet rad. : ${radius} Angstroms" | tee -a log.txt
echo "" | tee -a log.txt
echo "Num. reps:   : ${n_reps}" | tee -a log.txt
echo "*************************************************" | tee -a log.txt

echo ""

################################################################################################

# build keyfiles
################################################################################################    
echo "Building necessary new keyfiles..."
cp ${keyfile} "master_keyfile.key"

# firstly build the main keyfile
replace_line "FMCSC_XYZOUT"  "  FMCSC_XYZOUT ${xtcout}      # **** UPDATED BY AUTOSIM ****" master_keyfile.key
replace_line "FMCSC_TEMP"  "  FMCSC_TEMP ${temperature}     # **** UPDATED BY AUTOSIM ****" master_keyfile.key
replace_line "FMCSC_SIZE"  "  FMCSC_SIZE ${radius}          # **** UPDATED BY AUTOSIM ****" master_keyfile.key
replace_line "FMCSC_SEQFILE"  "  FMCSC_SEQFILE seq.in       # **** UPDATED BY AUTOSIM ****" master_keyfile.key
replace_line "FMCSC_PDBFILE"  "  FMCSC_PDBFILE start.pdb    # **** UPDATED BY AUTOSIM ****" master_keyfile.key
replace_line "FMCSC_PSWFILE"  "  FMCSC_PSWFILE PSWFILE.psw  # **** UPDATED BY AUTOSIM ****" master_keyfile.key

##############################################################################################################################
## Build coil-pre-equilibration
cp master_keyfile.key pre_eq.key
replace_line "FMCSC_NRSTEPS"  "  FMCSC_NRSTEPS ${preeq_steps}     # **** UPDATED BY AUTOSIM ****" pre_eq.key
replace_line "FMCSC_EQUIL"    "  FMCSC_EQUIL 0                 # **** UPDATED BY AUTOSIM ****" pre_eq.key
replace_line "FMCSC_RANDOMIZE"   "  FMCSC_RANDOMIZE 1                # **** UPDATED BY AUTOSIM ****" pre_eq.key
replace_line "FMCSC_SEQFILE"  "  FMCSC_SEQFILE ..\/seq.in       # **** UPDATED BY AUTOSIM ****" pre_eq.key
replace_line "FMCSC_TEMP"  "  FMCSC_TEMP ${preeq_temp}     # **** UPDATED BY AUTOSIM ****" pre_eq.key # note this overwrites master

##############################################################################################################################
## Creating helical equilibration keyfile if needed
if [ "${sim_mode}" = "combined" ] ||  [ "${sim_mode}" = "helix" ]
then
    cp pre_eq.key pre_eq_helix.key
    
    # wipe exising keywords IF they exits
    replace_line "FMCSC_SC_ZSEC"   "" pre_eq_helix.key
    replace_line "FMCSC_ZS_FR_A"   "" pre_eq_helix.key
    replace_line "FMCSC_ZS_FR_B"   "" pre_eq_helix.key
    replace_line "FMCSC_ZS_FR_KA"  "" pre_eq_helix.key

    # add new helical keywords
    echo "  FMCSC_SC_ZSEC   1.0     # **** UPDATED BY AUTOSIM **** [turn on sec. struct. potential]" >> pre_eq_helix.key
    echo "  FMCSC_ZS_FR_A   1.0     # **** UPDATED BY AUTOSIM **** [F-alpha scaling]" >> pre_eq_helix.key
    echo "  FMCSC_ZS_FR_B   0.0     # **** UPDATED BY AUTOSIM **** [F-beta scaling]" >> pre_eq_helix.key
    echo "  FMCSC_ZS_FR_KA  1000    # **** UPDATED BY AUTOSIM **** [spring constant for helical]" >> pre_eq_helix.key
fi

##############################################################################################################################
## create production keyfile
cp master_keyfile.key production.key
replace_line "FMCSC_NRSTEPS"  "  FMCSC_NRSTEPS ${prod_steps}     # **** UPDATED BY AUTOSIM ****" production.key
replace_line "FMCSC_EQUIL"    "  FMCSC_EQUIL ${eq_steps}         # **** UPDATED BY AUTOSIM ****" production.key
replace_line "FMCSC_RANDOMIZE"   "  FMCSC_RANDOMIZE 0                  # **** UPDATED BY AUTOSIM ****" production.key

replace_line "FMCSC_PDBFILE"  "" production.key
echo "   FMCSC_PDBFILE eq/__END.pdb    # **** UPDATED BY AUTOSIM **** [starting_config]" >> production.key
rm *bak


if [ "${sim_mode}" = "combined" ] ||  [ "${sim_mode}" = "helix" ]
then
    mkdir helical_start
    cd helical_start

    echo "Building helical_start directory structure..."

    # build local helical directories
    cp ../production.key .
    cp ../pre_eq_helix.key .
    cp ../seq.in .
    cp ../start.pdb .
    cp ../dres.in .
    cp ../PSWFILE.psw .
    for d in $(seq 1 ${n_reps})
    do
	mkdir $d
	cd $d
	cp ../production.key .
	cp ../pre_eq_helix.key .
	cp ../seq.in .
    cp ../start.pdb .
    cp ../dres.in .
    cp ../PSWFILE.psw .
	here_and_now=$(pwd)
	echo "#!/bin/bash

# script generated by autoSim
cd ${here_and_now}

mkdir eq
cd eq
cp ../pre_eq_helix.key .
cp ../dres.in .
cp ../PSWFILE.psw .
cp ../start.pdb .
${CAMPARI_BIN} -k pre_eq_helix.key | tee logfile_equilibration.out
touch COMPLETE
cd ..
${CAMPARI_BIN} -k production.key  | tee logfile_production.out
touch COMPLETE
"       > campari_bash.sh
	cd .. 
    done
    cd ..
fi

if [ "${sim_mode}" = "combined" ] ||  [ "${sim_mode}" = "coil" ]
then
    mkdir coil_start
    cd coil_start
    echo "Building coil_start directory structure..."

    # build local helical directories
    cp ../production.key .
    cp ../pre_eq.key .
    cp ../seq.in .
    cp ../start.pdb .
    cp ../dres.in .
    cp ../PSWFILE.psw .
    for d in $(seq 1 ${n_reps})
    do
	mkdir $d
	cd $d
	cp ../production.key .
	cp ../pre_eq.key .
	cp ../seq.in .
    cp ../start.pdb .
    cp ../dres.in .
    cp ../PSWFILE.psw .
	here_and_now=$(pwd)
	echo "#!/bin/bash

# script generated by autoSim
cd ${here_and_now}

mkdir eq
cd eq
cp ../pre_eq.key .
cp ../dres.in .
cp ../PSWFILE.psw .
cp ../start.pdb .
${CAMPARI_BIN} -k pre_eq.key | tee logfile_equilibration.out
touch COMPLETE
cd ..
${CAMPARI_BIN} -k production.key | tee logfile_production.out
touch COMPLETE
"       > campari_bash.sh
	cd .. 
    done
    cd ..
fi



####### CLEAN UP

if [ "${sim_mode}" = "helix" ]
then
    rm pre_eq.key
fi


rm master_keyfile.key
echo ""
echo "COMPLETE"
echo "*************************************************"

# added in 1.4
jn=$(basename $(pwd))
echo "${jn}" > JOB_PREFIX.txt
cp /work/holehouselab/tools/campari/setup_sim/run_seq.sh .

    
