import math

# Given center point coordinates for the line parallel to the z-axis
x0 = 0.616461
y0 = 0.348754
z0 = 4.676463  # Center z-coordinate

# Distance thresholds
max_distance = 20.0    # Main cutoff distance from the line
secondary_distance = 15.0  # Secondary cutoff distance for index listing
buffer_zone = 5.0  # Buffer zone for z-coordinate
cylinder_length = 72.950196  # Half-length of the cylinder

# Input and output file names
input_file = '/Users/uqtduign/Dropbox/Research/KcsA/data/input/md/waterbox/K_24-1k4c/kcsa.xyz'
output_file = '/Users/uqtduign/Dropbox/Research/KcsA/kcsa-mod-filtered-cylinder-4K.xyz'

# Atom type replacements
replacements = {
    'CO': 'C',
    'CA': 'C',
    'OH': 'O',
    'HO': 'H',
    'HN': 'H',
    'NH': 'N',
    'HS': 'H',
    'SH': 'S',
    'M+': 'K',
    'K+': 'K'
}

# Atom types to delete
delete_types = {'Cl-', 'Na+'}

# Lists to store atom data
all_atoms = []  # All atoms with their distances from the line

# Open and read the input file
with open(input_file, 'r') as file:
    lines = file.readlines()

# Initialize atom index
atom_index = -1  # Start from -1 since we increment after parsing each atom

# Process each line starting from the third line
i = 2
while i < len(lines):
    line = lines[i]
    if not line.strip():
        i += 1
        continue  # Skip empty lines
    parts = line.strip().split()
    if len(parts) < 5:
        i += 1
        continue  # Skip lines that don't contain atom data
    try:
        # Parse atom data
        atom_symbol = parts[1]
        x = float(parts[2])
        y = float(parts[3])
        z = float(parts[4])

        # Apply atom type replacements
        atom_symbol = replacements.get(atom_symbol, atom_symbol)

        # Skip unwanted atom types
        if atom_symbol in delete_types:
            i += 1
            continue

        # Increment atom index after successfully parsing an atom
        atom_index += 1

        # Calculate the distance from the line parallel to the z-axis
        distance = math.sqrt((x - x0)**2 + (y - y0)**2)

        atom_data = {
            'original_index': atom_index,
            'symbol': atom_symbol,
            'x': x,
            'y': y,
            'z': z,
            'distance': distance,
            'included': False  # Flag to indicate if the atom is included in the final output
        }

        all_atoms.append(atom_data)
    except ValueError:
        # Line might be a header or comment; skip it
        pass

    i += 1

# Lists to store final atoms and indexes to delete if cutoff was 15 Angstroms
final_atoms = []            # Atoms to include in the output file
indexes_to_delete_if_15 = []  # Indexes in the new XYZ file

# Initialize new atom index for the filtered list
new_index = -1  # Start from -1 to increment before adding to final_atoms

# Include atoms within max_distance directly
i = 0
while i < len(all_atoms):
    atom = all_atoms[i]
    if atom['original_index'] >= 5353 and atom['symbol'] == 'O' and i + 2 < len(all_atoms):
        # Check if this is a water molecule (O followed by two Hs)
        h1 = all_atoms[i + 1]
        h2 = all_atoms[i + 2]
        if h1['symbol'] == 'H' and h2['symbol'] == 'H':
            # Check the distance for the oxygen atom
            if atom['distance'] <= max_distance:
                # Increment new index for the water molecule
                new_index += 1
                atom['new_index'] = new_index
                atom['included'] = True
                final_atoms.append(atom)

                # Increment new index for the hydrogens
                for h in [h1, h2]:
                    new_index += 1
                    h['new_index'] = new_index
                    h['included'] = True
                    final_atoms.append(h)

                # Check if the molecule would be deleted if cutoff was 15 Angstroms
                if atom['distance'] > secondary_distance or \
                   atom['z'] < z0 - cylinder_length + buffer_zone or \
                   atom['z'] > z0 + cylinder_length - buffer_zone:
                    indexes_to_delete_if_15.extend([atom['new_index'], h1['new_index'], h2['new_index']])

            # Skip the entire water molecule
            i += 3
            continue

    # Process non-water atoms
    if atom['distance'] <= max_distance:
        # Increment new index
        new_index += 1
        atom['new_index'] = new_index
        atom['included'] = True
        final_atoms.append(atom)

        # Check if atom would be deleted if cutoff was 15 Angstroms
        if atom['distance'] > secondary_distance or \
           atom['z'] < z0 - cylinder_length + buffer_zone or \
           atom['z'] > z0 + cylinder_length - buffer_zone:
            indexes_to_delete_if_15.append(new_index)

    i += 1

# Update the atom count
num_atoms = len(final_atoms)

# Write the filtered atoms to the output XYZ file
with open(output_file, 'w') as file:
    file.write(f"{num_atoms}\n")
    file.write(f"Filtered atoms within {max_distance} Angstroms from the line parallel to the z-axis\n")
    for atom in final_atoms:
        file.write(f"{atom['symbol']} {atom['x']} {atom['y']} {atom['z']}\n")

# Print the list of new indexes of atoms that would be deleted if cutoff was 15 Angstroms
print("Indexes of atoms that would be deleted if cutoff was 15 Angstroms (in the new XYZ file):")
print(" ".join(map(str, indexes_to_delete_if_15)))
print(f"Total number of indices: {len(indexes_to_delete_if_15)}")

print(f"\nFiltered XYZ file has been saved as '{output_file}'.")