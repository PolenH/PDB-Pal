import workflowhead
from workflowhead import *
from collections import OrderedDict
import sys
import os

def euclidean_distance(point1, point2):
    point1 = np.array(point1)  # Convert to NumPy array, excluding the last element
    point2 = np.array(point2)  # Convert to NumPy array, excluding the last element
    return np.linalg.norm(point1 - point2)
#basically this function will take a list of cords and find the matching cords in the initial_system_dry file
#this will give us a list of atom indexes after they are run through orca and amber when combined into our one file
#for MM.
def findnewidx(list):
    newidx=[]
    with open('initial_system_dry.pdb', 'r') as file:
        for line in file:
            for cords in list:
                if line.startswith('HETATM') or line.startswith('ATOM'):
                    # Convert each coordinate to a string with exactly 3 digits after the decimal point
                    cords_str = [f"{cord:.3f}" for cord in cords]
                    # Check if all coordinates are in the line
                    if all(cord_str in line for cord_str in cords_str):
                        atom_index = int(line[6:11].strip())    #find matching cord and extract atom_idx
                        cords = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                        newidx.append([atom_index,cords])#add new list of new indexes
    return newidx

#find all atom residues given atom indexes in the initial_system_dry
def findR(atom_indexes):
    atoms=[]
    with open('initial_system_dry.pdb', 'r') as file:
        for line in file:
            if line.startswith('HETATM') or line.startswith('ATOM'):
                resid = int(line[22:26].strip())
                if resid in atom_indexes:
                    atom_id = int(line[6:11])
                    atoms.append(atom_id)

    return atoms
#here we find the atoms and cords of our intermediate structure. We also get the all atoms in INI
head, cords = FindHead('INI.pdb')
atoms = extract_all_atoms_manual('INI.pdb')
cords=cords.tolist()

acords=[]
hhead=[]
for atom in atoms:
    acords.append(atom[2])
    hhead.append(atom[0])

#now we need to find all the new atom ids for whole ligand and just head
newhead = findnewidx(cords)
newall = findnewidx(acords)

#now we need to find the atoms closest to head within three angstroms but make sure we keep head frozen
atomid = extract_all_atoms_manualW('initial_system_dry.pdb')
seenids=[]
resids=[]       #now we go through all atoms and only relax non head atoms that are within 3A of our Intermediate
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

#now we get rid of any duplicates then group them so we can write it to our inp file
seenids=sorted(OrderedDict.fromkeys(seenids))
grouped_numbers = group_consecutive_numbers(seenids)
range_strs = [f"{group[0] - 1}:{group[-1] - 1}" for group in grouped_numbers]
string=' '.join(range_strs)
#IT FRIKEN WORKS

#write to inp file so we know which atoms to be relaxed
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
