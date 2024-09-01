import workflowhead
from workflowhead import *
import os
import glob

# Get the input filename from user
filepath = "Input/" + next((file for file in os.listdir("Input") if file.endswith('.pdb')), None)

# Now we open the subdirectory within 'Input' and convert files

# Find all 'LIG.pdb' files in the specific subdirectory within 'Input'
with open(filepath, 'r') as original_file, open('LIG.pdb', 'w') as temp_file:
    for line in original_file:
        # Replace 'UNL' with 'LIG' in the current line
        new_line = line.replace('UNL', 'LIG')
        # Write the modified line to the temporary file
        temp_file.write(new_line)
#Now we have made LIG.pdb in correct format

#Now we have to run find head and align the protein

DFTopt = np.matrix([[-32.727,-37.694,21.911],[-34.009,-37.578,21.807],[-34.841,-38.437,22.123],[-34.526,-36.280,21.093],[-33.728,-35.523,20.564]])
DFTopID = ['O10','C16','O11','C5','O7']


cordlist=[]
headlist=[]
with open('Temp/output.txt', 'w') as file:
    head, cords = FindHead("LIG.pdb")
    
    print("filename: LIG.pdb, head: "+ str(head)+ ", cords: "+ str(cords),file=file)


#---------------------------------------------------------
# R,t,rmsd = kabsch_numpy(cordlist[i],DFTopt)
# print(cordlist)
R,t = kabsch_algorithm(cords,DFTopt)
atoms = extract_all_atoms_manual('LIG.pdb')

newcords = []
for atom in atoms:
    new_coords = np.dot(R, atom[2]) + t
    newcords.append(new_coords)

with open('LIG.pdb', 'r') as file:
    lines = file.readlines()

atom_lines = [line for line in lines if line.startswith('ATOM') or line.startswith('HETATM')]

with open('LIG.pdb', 'w') as file:
    for line, new_coord in zip(atom_lines, newcords):
        x, y, z = new_coord  # change this line
        line = line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[54:]
        file.write(line)

    last_hetatm_index = max(idx for idx, line in enumerate(lines) if line.startswith('ATOM') or line.startswith('HETATM'))

    for line in lines[last_hetatm_index+1:]:
        file.write(line)
    
#---------------------------------------------------------
#Now we have an aligned fully prepared ligand file
print("Ligand has been formatted and aligned, saved as LIG.pdb")
#Now we are done with python and need to run bash files to convert to amber format