#!/bin/zsh
source ~/.zshrc
source /work/bnovak/temp_sweep_test.zsh
##  Version 1.1 
##  --- Updated Sept 2020 - added priorty keyword
##
##


if [ "$#" = "1" ]
then
    if [ "$1" = "superlow" ]
    then
	p=$1
    elif [ "$1" = "low" ]
    then
	p=$1
    elif [ "$1" = "normal" ]
    then
	p=$1
    elif [ "$1" = "high" ]
    then
	p=$1
    elif [ "$1" = "--help" ] || [ "$1" = '-h' ]
    then
	echo "Run as  ./run_seq.sh <priority>"
	#echo -e "If provided priority must be one of:superlow\nlow\normal\nhigh"
	exit
    else
	
	echo "Invalid priority passed"
	exit
    fi
else
    p='normal'
fi



if [ ! -e JOB_PREFIX.txt ]
then
    echo "No JOB_PREFIX.txt file"
    exit
fi


JOB_PREFIX=$(cat JOB_PREFIX.txt)

for m in coil_start helical_start
do
    if [ -d ${m} ]
    then
	cd $m

	# get directory with biggest numerical starting value
	max_dir=$(ls -l | grep '^d' | awk {'print $9'} | sort -n | tail -n 1)

	# check this is a number..
	tst=$(expr $max_dir + 1)
	if [ $? != 0 ]
	then
	    here=$(pwd)
	    echo "Something went wrong figuring out directories in $here"
	    echo "Could there be a directory that starts with a number?"
	    exit
	fi
	
	for i in $(seq 1 $max_dir)
	do
	    if [ -d $i ]
	    then
		cd $i		    
		lsf_run_campari_v3_temperature_sweep ${JOB_PREFIX}_${m:0:4}_${i} ${p}
		cd ..
	    fi
	done
	cd ..
    fi
done
