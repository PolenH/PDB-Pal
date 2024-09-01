import numpy as np

def group_consecutive_numbers(numbers):
    #this function groups all of our consecutive numbers in a list this can then be converted to string
    #and used in the inp file
    ranges = []
    for number in numbers:
        if not ranges or number > ranges[-1][-1] + 1:
            ranges.append([number])
        else:
            ranges[-1].append(number)
    return ranges

def get_dist(cord1,cord2):
    #gets distance between two cords in 3d ONLY
    if cord1[2] == None or cord2[2] == None:
        print('Z cord is missing')
        return None
    d=np.sqrt((cord1[0]-cord2[0])**2+(cord1[1]-cord2[1])**2+(cord1[2]-cord2[2])**2)
    return d

def extract_all_atoms_manual(pdb_file):

    atoms = []
    #extracts all atoms basically stored as a tuple of atom_id, atom_name, and cords 
    #just pass pdb_file name
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
    #extracts all atoms but instead a tuple of atom_id, resid, and cords
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

def extract_atom_type(pdb_file,atom_name):
    #extracts all atoms then filters by atom type. This is more robust in future version in INI-STRUCT-5/ workflowhead.py
    #look there if your having issues with this function
    #the issue in this older version is if you using a pdb file with residue names that have numerics
    #like C10 instead of just C, it was updated in future versions to handle this.
    atoms = extract_all_atoms_manual(pdb_file)
    # print('pdb: '+pdb_file+' atoms: '+str(atoms))
    atomret=[]
    for atom in atoms:
        # print(atom[1])
        if atom[1]==atom_name:
            atomret.append(atom)
    return atomret

def find_matching_pairs(atomsnC, atomsnO, target_distance, tolerance):
    #finds matching pairs of atoms in two lists of atoms C and O basically but works with any two types
    #in our purposes we pass the C first since we are looking for C=O bonds not O=C bonds
    matching_pairs = []

    for atomC in atomsnC:
        for atomO in atomsnO:
            distance = get_dist(atomC[2], atomO[2])
            if abs(distance - target_distance) <= tolerance:
                matching_pairs.append((atomC, atomO))

    return matching_pairs
def find_matching_pairs_no_duplicates(atoms1, atoms2, target_distance, tolerance):
    #find matching pairs but gets rid of duplicates
    matching_pairs = []

    for i, atom1 in enumerate(atoms1):
        for j, atom2 in enumerate(atoms2):
            if i < j:
                distance = get_dist(atom1[2], atom2[2])
                if abs(distance - target_distance) <= tolerance:
                    matching_pairs.append((atom1, atom2))

    return matching_pairs
def FindHead(filename):
    #this function finds the head atoms and their cords in the pdb file
    #my suggestion in understanding how it works is to use print statements to see what is happening
    #basically it uses distances to find matching bonds in the pdb file
    #the series of for loops is the process in which I check distances and similar atoms to find those bonds.
    #It works even for new geometry seen in INI-STRUCT-5/, it has been updated to a protein standard
    #in workflowhead.py in INI-STRUCT-5/ That version might be easier to read.
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
    # print(atoms_by_id)
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
    #beautiful algorithm so keep it!
    return R, t, rmsd

def kabsch_algorithm(P, Q):
    #kabsch algorithm I wrote, heavily based on Heidenreich python implimentation 
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

