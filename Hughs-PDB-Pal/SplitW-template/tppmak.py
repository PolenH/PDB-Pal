#now we combine the tpp and lig files to one to make INI
from split import *
from Bio.PDB import PDBParser, PDBIO, Structure

def combinepdb(ligand, ligand2, output_pdb):
    # Initialize the PDB parser
    parser = PDBParser()
    residue_id = ''

    # Parse the first PDB file
    structure1 = parser.get_structure("structure1", ligand)

    # Parse the second PDB file
    structure2 = parser.get_structure("structure2", ligand2)

    # Create a new structure to hold both
    combined_structure = Structure.Structure("combined")

    # Add models from the first structure
    for model in structure1:
        combined_structure.add(model.copy())

    # Add models from the second structure
    # Ensure model IDs do not clash
    model_id_offset = len(combined_structure)
    for model in structure2:
        # Copy the model to avoid altering the original structure
        new_model = model.copy()
        # Update the model ID to ensure uniqueness
        new_model.id += model_id_offset
        combined_structure.add(new_model)

    # Get the residue ID of each atom in the model
    for chain in new_model:
        for i in chain.get_residues():
            residue_id = i.get_resname()
            break
    # Initialize PDBIO and set the structure
    io = PDBIO()
    io.set_structure(combined_structure)
    # Save the combined structure to a new PDB file
    io.save(output_pdb)
    with open(output_pdb, 'r') as file:
        lines = file.readlines()
    with open(output_pdb, 'w') as file:
        for line in lines:
            if (line.startswith("TER") and not line.__contains__(residue_id)) or line.startswith("ENDMDL") or line.startswith("MODEL      0"):
                continue
            file.write(line)
    #reorder lines
    residue_counter = 1

    # Read the file
    with open(output_pdb, 'r') as file:
        lines = file.readlines()

    # Open the file to write the changes
    with open(output_pdb, 'w') as file:
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Split the line into components
                parts = line.split()
                parts[1] = str(residue_counter)  # Update the integer part
                new_line = "{:<6}{:>5}  {:<4}{:<3} {:>5}    {:>8}{:>8}{:>8}  {:>4}  {:>4}           {:<2}\n".format(
                    parts[0], parts[1], parts[2], parts[3], parts[4],
                    parts[5], parts[6], parts[7], parts[8], parts[9], parts[10])
                file.write(new_line)
                residue_counter += 1
            elif line.startswith("TER"):
                parts = line.split()
                parts[1] = str(residue_counter)  # Update the integer part in TER line
                new_line = "{:<6}{:>5}      {:<4}{:>5}\n".format(parts[0], parts[1], parts[2], parts[3])
                file.write(new_line)
            else:
                file.write(line)


combinepdb("Output/ligands.pdb","Output/TPP.pdb","Output/ini.pdb")
#now I have to remove
# TER      21      LIG  1116                                                       
# ENDMDL
# MODEL      0
# MODEL 0 at top of file, all ENDMDL, reorder atomid
