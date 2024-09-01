from Bio import PDB

import re
#this is your code you should understand it

def replace_frags_section(input_file, output_file, new_data):
    # Read the content of the file
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Find the line with "FRAGS"
    for i, line in enumerate(lines):
        if line.startswith("    Frags"):
            break

    # Replace the line after the "FRAGS" line with new_data
    lines[i+1] = new_data + '\n'

    # Write the modified content to the output file
    with open(output_file, 'w') as file:
        file.writelines(lines)

def get_atoms_of_residues(pdb_file, residue_numbers):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)

    atom_indexes = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[1] in residue_numbers:
                    for atom in residue:
                        atom_indexes.append(atom.get_serial_number()-1) # adjust for orca indexing with 0

    return atom_indexes

def simplify_integer_list(int_list):
    if not int_list:
        return ""

    # Sort the list to make sure the integers are in ascending order
    sorted_list = sorted(int_list)
    ranges = []
    start = sorted_list[0]
    end = start

    for number in sorted_list[1:]:
        if number == end + 1:
            end = number
        else:
            if start == end:
                ranges.append(f"{start}")
            else:
                ranges.append(f"{start}:{end}")
            start = number
            end = start

    # Add the last range or number
    if start == end:
        ranges.append(f"{start}")
    else:
        ranges.append(f"{start}:{end}")

    return ' '.join(ranges)

def get_residue_numbers(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)

    residue_numbers = set()
    for model in structure:
        for chain in model:
            for residue in chain:
                residue_numbers.add(residue.id[1])

    return sorted(residue_numbers)


# 02 ext_2
residue_numbers = [31, 32, 78, 117, 118, 946, 947, 951, 968, 969, 970, 971, 972, 1031, 1177] # yi, QM 

#print(residue_numbers)

full_pdb_file = "initial_system_dry.pdb"
atom_indexes = get_atoms_of_residues(full_pdb_file, residue_numbers)
orca_list = simplify_integer_list(atom_indexes)

print(orca_list)

#_____________________________________________________________________________________

import workflowhead
from workflowhead import *


# O10 -32.727 -37.694 21.911
# C16 -34.009 -37.578 21.807
# O11 -34.841 -38.437 22.123
# C5 -34.526 -36.280 21.093
# O7 -33.728 -35.523 20.564

DFTopt = np.matrix([[-32.727,-37.694,21.911],[-34.009,-37.578,21.807],[-34.841,-38.437,22.123],[-34.526,-36.280,21.093],[-33.728,-35.523,20.564]])
DFTopID = ['O10','C16','O11','C5','O7']

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
grouped_numbers = group_consecutive_numbers(atom_indexes)
range_strs = [f"{group[0] - 1}:{group[-1] - 1}" for group in grouped_numbers]
print(' '.join(range_strs))

newlist="       2 {"+orca_list + ' ' + ' '.join(range_strs)+'} end'
input_file = 'mm_opt.inp'
output_file = 'mm_opt.inp'
# print(newlist)

replace_frags_section(input_file, output_file, newlist)
