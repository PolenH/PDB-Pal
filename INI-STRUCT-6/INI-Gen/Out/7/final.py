#run split.py, tppmak.py, then this py file
import workflowhead
from workflowhead import *
#open

#A keto acid head location
#ATOM     10  O   LIG  1116     -33.461 -37.829  22.530  1.00  1.00           O  
#ATOM      1  C   LIG  1116     -34.712 -37.546  22.439  1.00  1.00           C  
#ATOM     12  O2  LIG  1116     -35.657 -38.330  22.541  1.00  1.00           O  
#ATOM      2  C1  LIG  1116     -35.005 -35.981  22.098  1.00  1.00           C  
#ATOM     13  O3  LIG  1116     -34.054 -35.508  21.246  1.00  1.00           O  

#tail first carbon
# ATOM      3  C2  LIG  1116     -36.426 -35.820  21.476  1.00  1.00           C  

def get_vector(point1, point2):
    return np.array(point2) - np.array(point1)

keto_acid_head_location = np.array([-35.005,-35.981,22.098])  #4
tail_first_carbon = np.array([-36.426, -35.820, 21.476])  # tail

vector = get_vector(keto_acid_head_location, tail_first_carbon)

#   algorithm should follow like this:
#   get all non head atoms in ligand, then get the kabsch to align the 4th carbon and tail first carbon
#   then get the rotation matrix and translation vector and add to rest of file
#   Output a new file of INI know now as our INI structure


headlist=[]
newlist=[]
cordlist = []  # create a new cordlist
with open('Temp/output.txt', 'w') as file:
    head, cords = FindHead('LIG.pdb')
    atoms = extract_all_atoms_manual('LIG.pdb')
    
    headlist.append(head)
    for atom in atoms:
        if atom[0] not in head:
            newlist.append(atom[0])
            cordlist.append(atom[2])
    
   
    
    grouped_numbers = group_consecutive_numbers(newlist)
    range_strs = [f"{group[0] - 1}:{group[-1] - 1}" for group in grouped_numbers]
    print(range_strs, file=file)
    print(cordlist, file=file)
    print("filename: LIG.pdb, head: "+ str(head)+ ", cords: "+ str(cords), file=file)
    # Convert each sublist in cordlist to a string and join them with a space
atom_indexes = []
with open('initial_system_dry.pdb', 'r') as file:
    for line in file:
        for cords in cordlist:
            # Convert each coordinate to a string with exactly 3 digits after the decimal point
            cords_str = [f"{cord:.3f}" for cord in cords]
            # Check if all coordinates are in the line
            if all(cord_str in line for cord_str in cords_str):
                atom_index = int(line[6:11])
                atom_indexes.append(atom_index)
cc_arrayINI=[[-35.005,-35.981,22.098],[-36.426,-35.820,21.476]]
cc_array=[]
# print(atom_indexes)#all indexes of atoms in initialsysdry.pdb that are lig not head
atom=FindHead('LIG.pdb')
#now we need to find that carbon carbon bonds
fourthhead = np.array(atom[1][3])[0]
#distance is 1.517A
cc_array.append(fourthhead)

for cordfs in cordlist:
    if abs(get_dist(cordfs,fourthhead)) <= 1.517:
        cc_array.append(cordfs)

r,t=kabsch_algorithm(cc_array,cc_arrayINI)

# print(head)

# print(cordlist)
# print('crdlist: ')
# print(crdlist)
#now we have a list of new cords to add to the INI file and there correct numbers with them in LIG

with open('../Input/opt.pdb', 'r') as file:
    lines = file.readlines()
with open('opt.pdb', 'w') as file:
    head,cords,lastatom=FindHeadProtein('../Input/opt.pdb',0.1)
    for line in lines:
        if (line.startswith("ATOM") or line.startswith("HETATM")) and line.__contains__("LIG"):
            x, y, z = map(float, [line[30:38], line[38:46], line[46:54]])
            original_cord = np.array([x, y, z])
            id = int(line[6:11])
            if id in head[0]:
                file.write(line)
                continue
            if np.allclose(original_cord,cordlist,atol=0.5):
                file.write(line)
                continue
        else:
            file.write(line)
head,cords,lastatom=FindHeadProtein('LIG.pdb',0.1)
atoms = extract_all_atoms_manual('LIG.pdb')
headlist.append(head)
tail=list()
for atom in atoms:
    # print(atom[0],head)
    if atom[0] not in head[0]:
        tail.append(atom)
print(tail)
crdlist=[]
for cords in tail:
    cordn=np.dot(r, cords[2]) + t
    crdlist.append(cordn)
# print('cordlist: '+str(tail))
# print('crdlist: '+str(crdlist))

for i in range(len(tail)):
    # Create a new tuple with the updated third element
    coordinates_as_list = crdlist[i].tolist()
    # Create a new tuple with the updated third element as a list
    new_tuple = (tail[i][0], tail[i][1], coordinates_as_list)
    # Replace the old tuple with the new tuple in the list
    tail[i] = new_tuple

#__________________
# Assuming tail contains tuples with (atom_id, atom_name, [x, y, z]) and is already prepared
# Assuming tail contains tuples with (atom_id, atom_name, [x, y, z]) and is already prepared
# Function to get the next atom ID based on the highest ID in the file
def get_next_atom_id(lines):
    highest_id = 0
    for line in lines:
        if line.startswith("ATOM"):
            atom_id = int(line.split()[1])
            highest_id = max(highest_id, atom_id)
    return highest_id + 1

# Read the opt.pdb file
with open('opt.pdb', 'r') as file:
    lines = file.readlines()

# Get the next atom ID
next_atom_id = get_next_atom_id(lines)

# Generate new LIG atoms lines with updated atom IDs
new_lig_atoms_lines = []
for _, atom_name, coords in tail:
    element_symbol = atom_name[0] if atom_name[0].isalpha() else 'C'
    occupancy = "1.00"
    temp_factor = "0.00"
    pdb_line = f"ATOM  {next_atom_id:5d} {atom_name:<4} LIG  1116    {coords[0]:8.3f}{coords[1]:8.3f}{coords[2]:8.3f}  {occupancy:>4}{temp_factor:>6}           {element_symbol}\n"
    new_lig_atoms_lines.append(pdb_line)
    next_atom_id += 1

# Find the index of the last line containing a LIG residue
last_lig_index = 0
for i, line in enumerate(lines):
    if "LIG" in line:
        last_lig_index = i

# Insert the new LIG atoms lines after the last LIG residue
lines = lines[:last_lig_index+1] + new_lig_atoms_lines + lines[last_lig_index+1:]

# Write the updated content back to opt.pdb
with open('opt.pdb', 'w') as file:
    file.writelines(lines)
# # Write to opt.pdb
# with open('opt.pdb', 'w') as file:
#     file.writelines(lines)
#plan so far we are rewriting any the coords. we have the shifted cords in a list, just need to find old cords
#and replace them with new cords

#do this by checking where that carbon is that is in the opt.pdb then replacing all LIG cords with the new ones
#besides ones with head ids and that carbon
print("opt.pdb has been updated with new coordinates.")
