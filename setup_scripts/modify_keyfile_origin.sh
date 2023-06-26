#!/bin/zsh

replace_line()
{
    # $1 - replace string target
    # $2 - replace string replacement
    # $3 - file upon which replacement should happen
    sed -i.bak "s/.*${1}.*/${2}/" ${3}
}

helpFunction()
{
   echo ""
   echo "Usage: modify_keyfile_origin.sh -p PDB -k keyfile"
   echo -e "\t-p reference PDB file"
   echo -e "\t-k CAMPARI keyfile to modify the FMCSC_ORIGIN keyword"
   exit 1 # Exit script after printing help
}

while getopts "p:k:" opt
do
   case "$opt" in
      p ) pdbfile="$OPTARG" ;;
      k ) keyfile="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$pdbfile" ] || [ -z "$keyfile" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct
echo "PDB: $pdbfile"
echo "KEYFILE : $keyfile"

xfloat=$(awk '$1 == "ATOM" {sum += $7; n++ } END {if (n > 0) print sum / n; }' $pdbfile)
yfloat=$(awk '$1 == "ATOM" {sum += $8; n++ } END {if (n > 0) print sum / n; }' $pdbfile)
zfloat=$(awk '$1 == "ATOM" {sum += $9; n++ } END {if (n > 0) print sum / n; }' $pdbfile)

printf -v xint %.0f "$xfloat"  
printf -v yint %.0f "$yfloat" 
printf -v zint %.0f "$zfloat" 

replace_line "FMCSC_ORIGIN"  "  FMCSC_ORIGIN $xint $yint $zint     # === UPDATED BY modify_keyfile_origin.sh ===" $keyfile
rm "${keyfile}.bak"

