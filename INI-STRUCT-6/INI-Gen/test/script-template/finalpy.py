import workflowhead
from workflowhead import *
import os
import glob

#this finds the atom indexes of the ligand in the initial system
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
    # print(range_strs, file=file)
    # print(cordlist, file=file)
    # print(cordlist.__len__())
    print("filename: LIG.pdb, head: "+ str(head)+ ", cords: "+ str(cords), file=file)
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
# print(atom_indexes)
grouped_numbers = group_consecutive_numbers(atom_indexes)
range_strs = [f"{group[0] - 1}:{group[-1] - 1}" for group in grouped_numbers]
print(range_strs)
    