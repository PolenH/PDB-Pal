import numpy as np
import re

def group_consecutive_numbers(numbers):
    ranges = []
    for number in numbers:
        if not ranges or number > ranges[-1][-1] + 1:
            ranges.append([number])
        else:
            ranges[-1].append(number)
    return ranges

def get_dist(cord1,cord2):
    d=np.sqrt((cord1[0]-cord2[0])**2+(cord1[1]-cord2[1])**2+(cord1[2]-cord2[2])**2)
    return d

def extract_all_atoms_manual(pdb_file):
    atoms = []

    with open(pdb_file, 'r') as file:
        for line in file:
            # print(line)
            if line.startswith('HETATM') or line.startswith('ATOM'):
                # print(line)
                atom_id = int(line[6:11])
                atom_name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append((atom_id, atom_name, [x, y, z]))
                # print(line)
                # print(atom_id, atom_name, [x, y, z])
                # print(atoms)
    return atoms

def extract_all_atoms_manualW(pdb_file):
    atoms = []

    with open(pdb_file, 'r') as file:
        for line in file:
            # print(line)
            if line.startswith('HETATM') or line.startswith('ATOM'):
                # print(line)
                atom_id = int(line[6:11])
                atom_name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                resid = int(line[22:26].strip())
                atoms.append((atom_id, resid, [x, y, z]))
                # print(line)
                # print(atom_id, atom_name, [x, y, z])
                # print(atoms)
    return atoms

def rep_atoms_pdb(input_structure, rep_atoms):
    o = PDBIO()
    o.set_structure(input_structure)

    my_atoms = []
    for model in input_structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name() in rep_atoms:
                        my_atoms.append(atom)

    return my_atoms

def extract_atom_type(pdb_file,atom_name):
    atomret=[]
    atoms = extract_all_atoms_manual(pdb_file)
    pattern = re.compile(f'^{re.escape(atom_name)}\\d*$')

    # Append matching atoms to atomret
    for atom in atoms:
        if pattern.match(atom[1]):
            atomret.append(atom)
    return atomret

def find_matching_pairs(atomsnC, atomsnO, target_distance, tolerance):
    matching_pairs = []

    for atomC in atomsnC:
        for atomO in atomsnO:
            distance = get_dist(atomC[2], atomO[2])
            if abs(distance - target_distance) <= tolerance:
                matching_pairs.append((atomC, atomO))

    return matching_pairs
def find_matching_pairs_no_duplicates(atoms1, atoms2, target_distance, tolerance):
    matching_pairs = []

    for i, atom1 in enumerate(atoms1):
        for j, atom2 in enumerate(atoms2):
            if i < j:
                distance = get_dist(atom1[2], atom2[2])
                if abs(distance - target_distance) <= tolerance:
                    matching_pairs.append((atom1, atom2))

    return matching_pairs
def FindHead(filename):
    small_pdb = 'Head/head.pdb'
    large_pdb = filename

    atoms=extract_all_atoms_manual(small_pdb)
    # print(atoms)
    ato=[]
    ato.append(atoms[2])
    ato.append(atoms[0])
    ato.append(atoms[4])
    ato.append(atoms[1])
    ato.append(atoms[5])
    # print(ato)#I append the atoms in the order described in excel allignment file, chimera one.
    # print(ato[0][2])
    disthead=[]
    disthead.append(get_dist(ato[0][2],ato[1][2]))
    disthead.append(get_dist(ato[1][2],ato[2][2]))
    disthead.append(get_dist(ato[1][2],ato[3][2]))
    disthead.append(get_dist(ato[3][2],ato[4][2]))
    # print(disthead)#checked on chimera this works well for head pdb
    #check 2rd distance for large pdb C-O connections
    # print(disthead)


    atomsC=extract_atom_type(large_pdb,'C')
    # print(atomsC)
    atomsO=extract_atom_type(large_pdb,'O')

    matching_pairs = find_matching_pairs(atomsC, atomsO, disthead[1],0.1)

    #it properly found both pairs of C=O atoms in tolerance
    # print(matching_pairs)
    unique_pairs = find_matching_pairs(atomsC, atomsO, disthead[0],0.1)#change this to carbon single oxygen distance of 0.1A
    # unique_pairs = [pair for pair in matching if pair not in matching_pairs] #C-O pairs

    # Print the unique pairs
    # print(unique_pairs)
    # Find pairs of carbon atoms that are within a certain distance of each other
    matching_pairs_C_C = find_matching_pairs(atomsC, atomsC, disthead[2], 0.1) #all C-C pairs

    # Find all pairs in unique_pairs that share a carbon atom with a pair in matching_pairs_C_C
    matching_unique_pairs = []

    for pair1 in unique_pairs:
        for pair2 in matching_pairs:
            if pair1[0][0] == pair2[0][0] or pair1[0][0] == pair2[1][0]:
                matching_unique_pairs.append((pair1,pair2))
    # print(matching_unique_pairs)
    new_unique_pairs = []

    #adding last matching pairs
    for pair in matching_unique_pairs:
        for pairC_C in matching_pairs_C_C:
            # print("pairC: "+ str(pairC_C)+"/n")
            if pair[0][0][0] == pairC_C[0][0] or pair[0][0][0] == pairC_C[1][0]:
                new_unique_pairs.append((pair,pairC_C))
    # print(new_unique_pairs)
    # Find all pairs in matching_unique_pairs where pair[0] shares a carbon-carbon bond with a carbon that also has a C=O bond
    # use matching pairs
    # print(new_unique_pairs)
    newest_pair = []

    for pair in new_unique_pairs:
        for matching_pair in matching_pairs:
        
            if pair[-1][-1] == matching_pair[0]:
            
            
                newest_pair.append((pair,matching_pair))
    # print(newest_pair)
    final_pair = []

    for pair in newest_pair:
    
        if pair[0][0][0][0][0] != pair[-1][0][0]:
            final_pair.append(pair)

    # print(final_pair)
    atom_ids = []
    for pair in final_pair:
        atom_ids.append(pair[0][0][0][1][0])
        atom_ids.append(pair[0][0][1][0][0])
        atom_ids.append(pair[0][0][1][1][0])
        atom_ids.append(pair[1][0][0])
        atom_ids.append(pair[1][1][0])

    atoms_by_id = tuple(atom_ids)
    print(atoms_by_id)
    atom_cords = []
    for pair in final_pair:
        atom_cords.append(pair[0][0][0][1][2])
        atom_cords.append(pair[0][0][1][0][2])
        atom_cords.append(pair[0][0][1][1][2])
        atom_cords.append(pair[1][0][2])
        atom_cords.append(pair[1][1][2])

    atoms_by_cord = np.matrix(atom_cords)
    # print(atoms_by_cord)
    atoms_by_id = atoms_by_id[:5]
    atoms_by_cord = atoms_by_cord[:5]
    return atoms_by_id, atoms_by_cord

def FindHeadProtein(filename, error):
    small_pdb = 'Head/head.pdb'
    large_pdb = filename

    disthead=[]
    disthead.append(1.286)
    disthead.append(1.232)
    disthead.append(1.628)
    disthead.append(1.362)

    atomsC=extract_atom_type(large_pdb,'C')
    # print(atomsC)
    atomsO=extract_atom_type(large_pdb,'O')

    matching_pairs = find_matching_pairs(atomsC, atomsO, disthead[1],error)
    # print(matching_pairs)
    #it properly found both pairs of C=O atoms in tolerance
    
    unique_pairs = find_matching_pairs(atomsC, atomsO, disthead[0],error)#change this to carbon single oxygen distance of 0.1A
    # unique_pairs = [pair for pair in matching if pair not in matching_pairs] #C-O pairs
    # print(unique_pairs)
    # Print the unique pairs
    
    # Find pairs of carbon atoms that are within a certain distance of each other
    matching_pairs_C_C = find_matching_pairs(atomsC, atomsC, disthead[2],error) #all C-C pairs
    # print(matching_pairs_C_C)
    
    matching_unique_pairs = []
    used_oxygen_atoms = set()  # Corrected from used_carbon_atoms to used_oxygen_atoms

    for pair1 in unique_pairs:
        for pair2 in matching_pairs:
            if pair1==pair2:
                continue
            carbon_atom_pair1 = pair1[0][0]  # Assuming this is the carbon atom index
            carbon_atom_pair2 = pair2[0][0] if pair1[0][0] == pair2[0][0] else pair2[1][0]

            # Check if the carbon atoms match and the oxygen atom hasn't been used yet
            oxygen_atom_pair1 = pair1[1][0]  # Assuming this is the oxygen atom index in pair1
            if carbon_atom_pair1 == carbon_atom_pair2 and oxygen_atom_pair1 not in used_oxygen_atoms:
                # Determine which atom from pair2 matches the carbon atom and append the correct structure
                if pair1[0][0] == pair2[0][0]:
                    matching_unique_pairs.append((pair1, pair2))
                    used_oxygen_atoms.add(oxygen_atom_pair1)  # Track the used oxygen atom

    # if matching_unique_pairs:
    #     for pair in matching_unique_pairs:
    #         print(pair)
    #         print("-" * 20)  # Separator line
    # else:
    #     print("No matching unique pairs found.")

    riddup=[]
    seen = set()
    for pair in matching_unique_pairs:
        if pair[0][0][0] not in seen:
            riddup.append((pair[0][0],pair[0][1],pair[1][1]))
            seen.add(pair[0][0][0])
    # for pair in riddup:
    #     print(pair)
    #     print("-" * 20)  # Separator line
    fouratom=[]
    for pair in matching_pairs_C_C:
        # print('pair0: '+ str(pair[0])+ '\n riddup0: ' +str(riddup[0]))
        if pair[0] == riddup[0][0]:
            fouratom.append((riddup[0][0],riddup[0][1],riddup[0][2],pair[1]))
    # for pair in fouratom:
    #     print(pair)
    #     print("-" * 20)
    lastatom=[]
    for pair in unique_pairs:
        # print('pair0: '+ str(pair[0])+ '\n riddup0: ' +str(fouratom[0][3]))
        if pair[0][0] == fouratom[0][3][0]:
            lastatom.append((fouratom[0][1],fouratom[0][0],fouratom[0][2],fouratom[0][3],pair[1]))
    # print(lastatom[0][0])
    # print(lastatom[0][1])
    # print(lastatom[0][2])
    # print(lastatom[0][3])
    # print(lastatom[0][4])
    atom_ids = []
    atom_ids.append((lastatom[0][0][0],lastatom[0][1][0],lastatom[0][2][0],lastatom[0][3][0],lastatom[0][4][0]))
    atom_crds = []
    atom_crds.append((lastatom[0][0][2],lastatom[0][1][2],lastatom[0][2][2],lastatom[0][3][2],lastatom[0][4][2]))
    return atom_ids, atom_crds, lastatom

def kabsch_numpy(P, Q):
    """
    Credit: Hunter Heidenreich https://hunterheidenreich.com/posts/kabsch_algorithm/
    Computes the optimal rotation and translation to align two sets of points (P -> Q),
    and their RMSD. Dont work anymore

    :param P: A Nx3 matrix of points
    :param Q: A Nx3 matrix of points
    :return: A tuple containing the optimal rotation matrix, the optimal
             translation vector, and the RMSD.
    """
    assert P.shape == Q.shape, "Matrix dimensions must match"
    # print(P)
    # print('-----')
    # print(Q)
    # Compute centroids
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)

    # Optimal translation
    t = centroid_Q - centroid_P

    # Center the points
    p = P - centroid_P
    q = Q - centroid_Q

    # Compute the covariance matrix
    H = np.dot(p.T, q)

    # SVD
    U, S, Vt = np.linalg.svd(H)

    # Validate right-handed coordinate system
    if np.linalg.det(np.dot(Vt.T, U.T)) < 0.0:
        Vt[-1, :] *= -1.0

    # Optimal rotation
    R = np.dot(Vt.T, U.T)

    # RMSD
    rmsd = np.sqrt(np.sum(np.square(np.dot(p, R.T) - q)) / P.shape[0])

    return R, t, rmsd

def kabsch_algorithm(P, Q):
    # Ensure the points are numpy arrays
    P = np.array(P)
    Q = np.array(Q)

    # Step 1: Compute the centroids of P and Q
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)

    # Step 2: Center the points
    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q

    # Step 3: Compute the covariance matrix
    H = np.dot(P_centered.T, Q_centered)

    # Step 4: Compute the optimal rotation matrix using SVD
    U, S, Vt = np.linalg.svd(H)
    R = np.dot(Vt.T, U.T)

    # Special reflection case handling
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = np.dot(Vt.T, U.T)

    # Step 5: Compute the translation vector
    t = centroid_Q - np.dot(R, centroid_P)

    return R, t

