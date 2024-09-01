#!/bin/bash
#before running script do python, hit enter, then do pip install numpy, and then exit.
#ONLY manual thing is selecting charges with -nc flag. use it if you want to change it to something other than -1
#make sure your in a conda environment and this should work otherwise try uncommenting the conda activate line


#now we are converting to amber format using antechamber
#then generating an orca file
module purge
module load orca
module load ambertools
# module load python/anaconda3.6

# conda activate researchQ

#call fheadlig.py
python fheadlig.py

# Initialize variables
nc_value=-1  # Default value if -nc is not provided

# sets our nc value by taking an optional -nc flag when calling script.sh
while [[ "$#" -gt 0 ]]; do
  case $1 in
    -nc)
      if [[ $2 =~ ^- ]]; then
        nc_value="$2"
      else
        nc_value=-1
      fi
      shift 2
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
done
echo "nc_value is $nc_value"
# Example function that uses the parsed value
pdb4amber -i LIG.pdb -o LIG_amber.pdb
antechamber -i LIG_amber.pdb -o LIG.mol2 -fi pdb -fo mol2 -c bcc -pf yes -nc "$nc_value" -at gaff2 -j 5
parmchk2 -i LIG.mol2 -f mol2 -o LIG.frcmod -s 2

# Run tleap and capture its output
tleap_output=$(tleap -f tleap_script.in | tee /dev/tty)

# Use echo to pass the output to grep, then extract and round the charge
# Extract the charge value and remove parentheses
charge=$(echo "$tleap_output" | grep "The unperturbed charge of the unit" | awk '{print $7}' | tr -d '()')

# Print the charge without parentheses
echo "Charge: $charge"
# Use echo to pass the output to grep, then extract and correctly round the charge
rounded_charge=$(awk -v charge="$charge" 'BEGIN {rounded = sprintf("%.0f", charge); print rounded}')

# Print the rounded charge
echo "Rounded charge: $rounded_charge"

#-----
#calls final3a.py to change inp file
python final3a.py "$rounded_charge"

echo "File updated successfully."
orca_mm -convff -AMBER complex_dry.prmtop

#last is submitting the job with $ sbatch mm_orca_job.sh