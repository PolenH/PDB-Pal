import workflowhead
from workflowhead import *
from collections import OrderedDict
import sys
import os

def euclidean_distance(point1, point2):
    point1 = np.array(point1)  # Convert to NumPy array, excluding the last element
    point2 = np.array(point2)  # Convert to NumPy array, excluding the last element
    return np.linalg.norm(point1 - point2)

def findnewidx(list):
    newidx=[]
    with open('initial_system_dry.pdb', 'r') as file:
        for line in file:
            for cords in list:
                # Convert each coordinate to a string with exactly 3 digits after the decimal point
                cords_str = [f"{cord:.3f}" for cord in cords]
                # Check if all coordinates are in the line
                if all(cord_str in line for cord_str in cords_str):
                    atom_index = int(line[6:11])
                    cords = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                    newidx.append([atom_index,cords])
    return newidx

def findR(atom_indexes):
    atoms=[]
    with open('initial_system_dry.pdb', 'r') as file:
        for line in file:
            if line.startswith('HETATM') or line.startswith('ATOM'):
                resid = int(line[22:26].strip())
                if resid in atom_indexes:
                    atom_id = int(line[6:11])
                    # atom_name = line[12:16].strip()
                    # x = float(line[30:38])
                    # y = float(line[38:46])
                    # z = float(line[46:54])
                    atoms.append(atom_id)

    return atoms
#extracts head atoms and cords. Also extracts all atoms and cords. Makes cords a list for consistency
#These functions are stored in workflowhead.py
head, cords = FindHead('LIG.pdb')
atoms = extract_all_atoms_manual('LIG.pdb')
cords=cords.tolist()

acords=[]
hhead=[]
#seperates cords and atom ids
for atom in atoms:
    acords.append(atom[2])
    hhead.append(atom[0])

#now we need to find all the new atom ids for whole ligand and just head, functions from workflowhead.py
newhead = findnewidx(cords)
newall = findnewidx(acords)

#now we need to find the atoms closest to head within three angstroms but make sure we keep head frozen
atomid = extract_all_atoms_manualW('initial_system_dry.pdb')
seenids=[]
resids=[]
for ids in atomid:
    if ids[2] in newhead:
        continue
    if ids[0] in seenids:
        continue
    for i in newall:
        distance=euclidean_distance(ids[2],i[1])
        if (distance <= 3):  # Corrected syntax error here
            resids.append(ids[1])
            seenids.append(ids[0])
idfr= findR(resids)
seenids.extend(idfr)

#now we need to group these atoms using a set python method which will get rid of copies as well as sort it
seenids=sorted(OrderedDict.fromkeys(seenids))
grouped_numbers = group_consecutive_numbers(seenids)
#then we group it so it can be used to replace in the inp file
range_strs = [f"{group[0] - 1}:{group[-1] - 1}" for group in grouped_numbers]
string=' '.join(range_strs)
#IT FRIKEN WORKS

#now we edit inp file with the new charge and atomids
rounded_ch = sys.argv[1]
with open('mm_opt.inp','r') as file:
    lines = file.readlines()

with open('tmp.inp', 'w') as file:
    found_pdbfile = False
    found_frags = False

    for line in lines:
        if line.startswith('*pdbfile') and not found_pdbfile:
            parts = line.split()
            parts[1] = str(rounded_ch)  # Update the charge
            line = ' '.join(parts) + '\n'
            found_pdbfile = True
        
        if line.strip() == "Frags":
            found_frags = True
            file.write(line)
            continue

        if found_frags and line.strip().startswith('2 {'):
            line = f"       2 {{{string}}} end\n"
            found_frags = False  # Assuming only one Frags section to modify

        file.write(line)
os.replace('tmp.inp','mm_opt.inp')
