import os

if os.path.exists('mm_opt.pdb'):
    filename = 'mm_opt.pdb'
else:
    print("Both mm_opt.pdb does not exist.")
    exit()
# Rest of the code...
with open(filename, 'r') as file:
    lines = file.readlines()
with open(filename, 'w') as file:
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            atom_type = line[16:20].strip()  # Extract and strip the type field
            
            if atom_type == 'TPP' or atom_type == 'LIG':
                atom_type = 'INI'
            
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