from Bio.PDB import PDBParser, PDBIO, Select
import workflowhead
from workflowhead import *
import os
class LigandSelect(Select):
    def __init__(self, residue_names):
        self.residue_names = residue_names

    def accept_residue(self, residue):
        return residue.get_resname() in self.residue_names

# Define the input and output file paths
input_pdb_file = "opt.pdb"  # Replace with your input PDB file path
output_pdb_file = "INI-pre.pdb"  # Replace with your desired output PDB file path

# Initialize the PDB parser
parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", input_pdb_file)

# Define the residue names to filter
residue_names = {"LIG", "TPP"}

# Initialize the PDBIO writer
io = PDBIO()
io.set_structure(structure)
io.save(output_pdb_file, LigandSelect(residue_names))

# print(f"Filtered residues saved to {output_pdb_file}")

id,crds,lastatom=FindHeadProtein(output_pdb_file,0.1)
# print(id)
with open('INI-pre.pdb', 'r') as file:
    lines = file.readlines()

with open('INI.pdb', 'w') as file:
    for line in lines:
        if (line.startswith('ATOM') or line.startswith('HETATM')) and line.__contains__('LIG'):
            line_parts = line.split()
            atom_id = int(line_parts[1])  # Assuming the ID is always in the second position
            if atom_id not in id[0]:
                continue
        file.write(line)

#now every non tail atom is in INI.pdb 

#now align the tail and insert it

head, cords = FindHead('LIG.pdb')
atoms = extract_all_atoms_manual('LIG.pdb')
cordlist = []
newlist=[]
# headlist.append(head)
for atom in atoms:
    if atom[0] not in head:
        newlist.append(atom[0])
        cordlist.append(atom[2])

cc_arrayINI=[[-35.005,-35.981,22.098],[-36.426,-35.820,21.476]]
cc_array=[]
# print(atom_indexes)#all indexes of atoms in initialsysdry.pdb that are lig not head
atom=FindHead('LIG.pdb')
#now we need to find that carbon carbon bonds
fourthhead = np.array(atom[1][3])[0]
#distance is 1.517A
cc_array.append(fourthhead)

for cordfs in cordlist:
    if abs(get_dist(cordfs,fourthhead)) <= 1.526:
        cc_array.append(cordfs)
current_directory = os.path.dirname(os.path.abspath(__file__))

# Print the directory
print(f"The script is running in: {current_directory}")
# print('cc array: ' + str(cc_array)+ ' cc arrayINI: ' + str(cc_arrayINI))
r,t=kabsch_algorithm(cc_array,cc_arrayINI)

head,cords=FindHead('LIG.pdb')
atoms = extract_all_atoms_manual('LIG.pdb')
# headlist.append(head)
tail=list()
for atom in atoms:
    # print(atom[0],head)
    if atom[0] not in head:
        tail.append(atom)
# print(tail)
crdlist=[]
for cords in tail:
    cordn=np.dot(r, cords[2]) + t
    crdlist.append(cordn)
print(cordfs)
def find_last_atom_id(file_path):
    last_atom_id = 0
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                try:
                    atom_id = int(line[6:11].strip())
                    if atom_id > last_atom_id:
                        last_atom_id = atom_id
                except ValueError:
                    continue
    return last_atom_id

new_atom_id = find_last_atom_id('INI.pdb') + 1
with open('LIG.pdb','r') as file:
    lines = file.readlines()
# Open INI.pdb for appending
with open('INI.pdb', 'a') as ini_file:
    for i, (atom_index, *_) in enumerate(tail):
        # Get the corresponding line from LIG.pdb
        line = lines[atom_index - 1]
        
        # Split the line into parts
        line_parts = line.split()
        
        # Update the coordinates
        new_x = f"{crdlist[i][0]:8.3f}"
        new_y = f"{crdlist[i][1]:8.3f}"
        new_z = f"{crdlist[i][2]:8.3f}"
        number=1116
        line_parts[5] = str(number).rjust(4)
        
        # Construct the new line with updated coordinates and residue sequence number
        new_line = f"{line[:6]}{str(new_atom_id).rjust(5)}{line[11:22]}{line_parts[5]}{line[26:30]}{new_x:>8}{new_y:>8}{new_z:>8}{line[54:]}"
        
        # Write the modified line to INI.pdb
        ini_file.write(new_line)
        
        # Increment the atom ID for the next line
        new_atom_id += 1

with open('INI.pdb', 'r') as file:
    lines = file.readlines()
with open('INI.pdb', 'w') as file:
    ter=''
    end=''
    for line in lines:
        if line.startswith('TER'):
            ter=line
            continue
        if line.startswith('END'):
            end=line
            continue
        file.write(line)
    file.write(ter)
    file.write(end)

with open('INI.pdb', 'r') as file:
    lines = file.readlines()
with open('INI.pdb', 'w') as file:
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            atom_type = line[16:20].strip()  # Extract and strip the type field
            
            # if atom_type == 'TPP' or atom_type == 'LIG':
            #     atom_type = 'INI'
            
            # Construct the new line with preserved formatting
            new_line = f"{line[:17]}{atom_type.ljust(4)}{line[21:]}"
            file.write(new_line)
        elif line.startswith('TER'):
            ter_type_start = 17
            ter_type_end = 20
            new_line = line[:ter_type_start] + 'INI'.ljust(ter_type_end - ter_type_start) + line[ter_type_end:]
            
            file.write(new_line)
        else:
            file.write(line)

        