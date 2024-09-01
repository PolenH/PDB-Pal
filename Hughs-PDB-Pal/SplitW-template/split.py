import Bio
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, Select
#go into opt.pdb
parser = PDB.PDBParser(QUIET=True)
io = PDB.PDBIO()
#open it up and read it. Write 5 files, protein.pdb, ligand.pdb, TPP.pdb, MG.pdb, and water.pdb.
#find the residues by name, LIG, TPP, MG, and HOH then seperate them while preserving correct connections
def separatepdb(input_pdb, output_pdb,resname):
    structure = parser.get_structure("structure", input_pdb)
    
    ligands = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname().startswith(resname):
                    ligands.append(residue)
    
    class LigandSelect(PDB.Select):
        def accept_residue(self, residue):
            return residue in ligands

    io.set_structure(structure)
    io.save(output_pdb, LigandSelect())

input_pdb = "Input/opt.pdb"

separatepdb(input_pdb,"Output/ligands.pdb","LIG")
separatepdb(input_pdb,"Output/TPP.pdb","TPP")
separatepdb(input_pdb,"Output/MG.pdb","MG")
separatepdb(input_pdb,"Output/solvent.pdb","WAT")

#now we take all these files and output a protein file that doesn't have a TPP, MG, LIG, or WAT residue

excluded_residues = ["TPP","MG","LIG","WAT"]
class ResidueSelect(Select):
    def __init__(self, excluded_residues):
        self.excluded_residues = excluded_residues

    def accept_residue(self, residue):
        return residue.get_resname().strip() not in self.excluded_residues

# Assuming 'structure' is your loaded PDB structure
structure = parser.get_structure("ID", "Input/opt.pdb")

# Define the residues to exclude
excluded_residues = ["TPP", "MG", "LIG", "WAT"]

io.set_structure(structure)

# Use the custom Select class to filter out unwanted residues
io.save("Output/protein.pdb", ResidueSelect(excluded_residues))
#now test with tleap

# Open the original PDB file for reading
with open('Output/protein.pdb', 'r') as file:
    lines = file.readlines()

# Prepare a list to hold the modified lines
modified_lines = []

for index, line in enumerate(lines):
    # Add the current line to the modified lines list
    modified_lines.append(line)
    # Check if the line contains an OXT atom
    if " OXT " in line:
        # Check if this is not the last line and if the next line contains "TER"
        if index + 1 < len(lines) and "TER" not in lines[index + 1]:
            # Extract necessary information to construct the TER line
            fields = line.split()
            atom_id = fields[1]
            res_name = fields[3]
            res_id = fields[4]
            # Construct and append the TER line
            ter_line = f"TER   {atom_id:>5}      {res_name:<3} {res_id}\n"
            modified_lines.append(ter_line)
        elif index + 1 == len(lines):  # If it's the last line, just add the TER line
            fields = line.split()
            atom_id = fields[1]
            res_name = fields[3]
            res_id = fields[4]
            ter_line = f"TER   {atom_id:>5}      {res_name:<3} {res_id}\n"
            modified_lines.append(ter_line)

# Write the modified lines back to a new file or overwrite the original
with open('Output/protein_modified.pdb', 'w') as file:
    file.writelines(modified_lines)