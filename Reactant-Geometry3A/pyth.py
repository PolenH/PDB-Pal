import os
import shutil

# Paths to the source and destination directories
script_template_path = 'script-template'
pdbs_nonaligned_path = 'PDBSnonaligned'
destination_root = 'Out'  # This is where the script0-19 folders will be created

# Ensure the destination directory exists
os.makedirs(destination_root, exist_ok=True)

# Loop to create script0 through script19 folders
for i in range(20):
    # Define new script folder name
    new_script_folder = os.path.join(destination_root, f'script{i}')
    
    # Copy the template folder to the new script folder
    shutil.copytree(script_template_path, new_script_folder)
    
    # Define the path to the Input folder inside the new script folder
    input_folder = os.path.join(new_script_folder, 'Input')
    
    # Ensure the Input folder exists (it should if it's in the template)
    os.makedirs(input_folder, exist_ok=True)
    
    # Define the path to the corresponding PDB file
    pdb_file = os.path.join(pdbs_nonaligned_path, f'{i}.pdb')
    
    # Define the destination path for the PDB file
    destination_pdb_file = os.path.join(input_folder, f'{i}.pdb')
    
    # Copy the PDB file to the Input folder
    shutil.copy(pdb_file, destination_pdb_file)

print("Folders and files have been created and copied successfully.")
