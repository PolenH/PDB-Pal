#!/bin/bash

#Instructions:
#Drop scripts and ONLY pdb files in a seperate directory, this script will convert all pdb files in directory. Make sure you only have ligand pdb files in your directory

#This part of script is taken from ADFR_prep_lig.sh written by Geoff
export BABEL_DATADIR="/mnt/c/Program Files (x86)/ADFRsuite-1.0/OpenBabel-2.4.1/data"
export PATH="/mnt/c/Program Files (x86)/ADFRsuite-1.0/bin:/mnt/c/Program Files (x86)/ADFRsuite-1.0/OpenBabel-2.4.1:$PATH"

# Get the name of the convScript
script_name=$(basename "$0")

# Directory containing the files and script
script_dir="$(dirname "$0")"

for file in "$script_dir"/*
do
    # Get the basename of the current file, file we are currently indexing within the folder
    filename=$(basename "$file")

    if [[ "$filename" != "$script_name" && -f "$file" && "$file" == *.pdb && "$file" != "$optnolig.pdb" ]]; then    #checks if we are processing the scripts file and if the file exists and is a regular file
        echo "Processing file: $filename" #debugging command, you can comment this out if you want
        # edited from Geoff's original script, replace with receptor part of command if mass compilling receptors
        "/mnt/c/Program Files (x86)/ADFRsuite-1.0/python.exe" "C:\Program Files (x86)\ADFRsuite-1.0\Lib\site-packages\AutoDockTools\Utilities24\prepare_ligand4.py" -l $file -U -o test
    fi
done
