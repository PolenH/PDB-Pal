#!/bin/bash
#we need to copy all of the up-to-date INI files and compile them in tleap

# Navigate to the INI-Gen/Out directory
cd INI-Gen/Out

# Iterate over each folder in the current directory
for folder in */; do
    # Create the same folder structure in the Out directory
    mkdir -p "../../Out/${folder}"
    
    # Check if INI-rot.pdb exists in the folder
    if [[ -f "${folder}INI-rot.pdb" ]]; then
        # If INI-rot.pdb exists, copy it to the Out directory and rename it to INI.pdb
        cp "${folder}INI-rot.pdb" "../../Out/${folder}INI.pdb"
    elif [[ -f "${folder}INI.pdb" ]]; then
        # If INI-rot.pdb does not exist but INI.pdb does, copy INI.pdb
        cp "${folder}INI.pdb" "../../Out/${folder}INI.pdb"
    fi
done

# Navigate back to the original directory
cd - > /dev/null

for folder in Out/*/; do
    # Copy all content from the script/ directory to the current folder in Out/
    cp -r script/* "${folder}"
done