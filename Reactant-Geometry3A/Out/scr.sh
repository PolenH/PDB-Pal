#!/bin/bash

# Assuming this script is located in the parent directory of script0, script1, ..., script19
parent_dir=$(dirname "$0")
scripts_dir="$parent_dir"  # Adjust this if scripts are in a specific subdirectory like "scripts"

# Loop through script directories script0 to script19
for i in {0..19}; do
    script_dir="script$i"
    if [ -d "$scripts_dir/$script_dir" ]; then
        echo "Entering directory: $script_dir"
        cd "$scripts_dir/$script_dir" || exit 1
        
        # Check conditions for i == 4, 11, 16
        if [ "$i" -eq 4 ] || [ "$i" -eq 11 ] || [ "$i" -eq 16 ]; then
            # Execute script.sh with specific options
            bash script.sh -nc -2 && sbatch mm_orca_job.sh
        else
            # Execute default command
            bash script.sh && sbatch mm_orca_job.sh
        fi
        
        # Optionally, you can add a delay or any other logic between script executions
        echo "Execution in $script_dir completed."
        echo
        
        # Move back to parent directory
        cd .. || exit 1
    else
        echo "Directory $script_dir not found."
    fi
done
