import matplotlib.pyplot as plt
from ase.io import read
import re
import os
import numpy as np

# Define the set_plot_style function before it is used
def set_plot_style(ax=None):
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 32
    plt.tight_layout()
    if ax is None:
        ax = plt.gca()
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(3)
    ax.spines['left'].set_linewidth(3)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.tick_params(axis='both', which='major', labelsize=24, length=10, width=3)
# List of your trajectory files in order
trajectory_files = [
    '/Users/uqtduign/Desktop/temp/Trajectories/KcsA/KcsA_4K_1.00e-01F.xyz',
    '/Users/uqtduign/Desktop/temp/Trajectories/KcsA/KcsA_4K_1.00e-01F-r1.xyz',
    '/Users/uqtduign/Desktop/temp/Trajectories/KcsA/KcsA_4K_1.00e-01F-r2.xyz',
    '/Users/uqtduign/Desktop/temp/Trajectories/KcsA/KcsA_4K_1.00e-01F-r3.xyz',
    '/Users/uqtduign/Desktop/temp/Trajectories/KcsA/KcsA_4K_1.00e-01F-r4.xyz',
    '/Users/uqtduign/Desktop/temp/Trajectories/KcsA/KcsA_4K_1.00e-01F-r5.xyz',
    '/Users/uqtduign/Desktop/temp/Trajectories/KcsA/KcsA_4K_1.00e-01F-r6.xyz'
]

# Limit the number of frames to read from each trajectory
max_frames_per_traj = 12500  # Process only the first 100 frames of each trajectory

# Initialize lists to store time and z-positions
times = []
z_positions = {}

# Regular expression to extract time from comment line
time_regex = re.compile(r'time:\s*([-\d\.eE]+)\s*fs')


# Total time accumulated from previous trajectories
cumulative_time = 0.0

# Flag to check if K atom indices have been initialized
k_indices = None

# List of water molecule indices
water_indices = [4462, 4594, 4216, 4123]

# Initialize z_positions for water molecules
water_z_positions = {i: [] for i in water_indices}

# Pairs of atom indices for RDF calculation
atom_pairs = [
    (1609, 1630),
    (2662, 2683),
    (555, 576),
    (3714, 3735)
]

# Initialize a list to store distances for RDF calculation
distances = {pair: [] for pair in atom_pairs}

# Open a text file to write positions
with open('z_positions_all.txt', 'w') as f_out:
    # Header for the output file
    header_written = False

    # Initialize a list to store RDF results for all trajectories
    all_rdf_results = []

    for traj_file in trajectory_files:
        if not os.path.isfile(traj_file):
            print(f"File not found: {traj_file}")
            continue

        # Read the trajectory file lines directly to extract times
        with open(traj_file, 'r') as f:
            lines = f.readlines()

        # Calculate the number of atoms from the first line
        natoms = int(lines[0].strip())
        frame_length = natoms + 2  # Number of lines per frame (atoms + 2 header lines)

        # Calculate the total number of frames in the file
        total_frames = len(lines) // frame_length

        # Adjust 'max_frames_per_traj' if it's larger than the total frames available
        num_frames = min(max_frames_per_traj, total_frames)

        # Read frames using ASE
        traj = read(traj_file, index=f'0:{num_frames}')

        # Initialize K atom indices and z_positions dict
        if k_indices is None:
            # Get indices of K atoms from the first frame
            first_frame = traj[0]
            k_indices = [atom.index for atom in first_frame if atom.symbol == 'K']
            # Initialize z_positions dict
            z_positions = {i: [] for i in range(len(k_indices))}

            # Write header to output file
            f_out.write('Time(fs)\t' + '\t'.join([f'K{i+1}_Z(Å)' for i in range(len(k_indices))]) + '\n')
            header_written = True

        # Loop over each frame
        for frame_number in range(num_frames):
            # Index of the comment line in the file
            comment_line_idx = frame_number * frame_length + 1
            comment_line = lines[comment_line_idx].strip()

            # Extract time from the comment line
            time_match = time_regex.search(comment_line)
            if time_match:
                time = float(time_match.group(1))
            else:
                time = None  # Handle missing time information if necessary

            # Adjust time with cumulative time
            if time is not None:
                adjusted_time = time + cumulative_time
            else:
                # Estimate time if missing
                if times:
                    # Assume constant time step
                    time_step = times[-1] - times[-2] if len(times) > 1 else 0.0
                    adjusted_time = times[-1] + time_step
                else:
                    adjusted_time = cumulative_time  # Start from cumulative_time

            times.append(adjusted_time)

            # Get the corresponding frame from the trajectory
            frame = traj[frame_number]

            # Extract z-positions of K atoms
            z_pos = []
            for i, k_index in enumerate(k_indices):
                z = frame.positions[k_index][2]  # z-coordinate
                z_positions[i].append(z)
                z_pos.append(z)

            # Extract z-positions of specified water molecules
            water_z_pos = []
            for water_index in water_indices:
                z = frame.positions[water_index][2]  # z-coordinate
                water_z_positions[water_index].append(z)
                water_z_pos.append(z)

            # Calculate distances for each pair of atoms
            for pair in atom_pairs:
                atom1, atom2 = pair
                pos1 = frame.positions[atom1]
                pos2 = frame.positions[atom2]
                distance = np.linalg.norm(pos1 - pos2)
                distances[pair].append(distance)

            # Write to output file
            if header_written:
                line = f"{adjusted_time}\t" + '\t'.join([f"{z}" for z in z_pos]) + '\t' + '\t'.join([f"{z}" for z in water_z_pos]) + '\n'
                f_out.write(line)

        # Update cumulative time to the last adjusted time of the current trajectory
        if times:
            cumulative_time = times[-1]

        # Compute RDF from distances for the current trajectory
        # Define bins for RDF
        bin_edges = np.linspace(0, 10, 100)  # Adjust range and number of bins as needed
        rdf_results = []

        for pair, dists in distances.items():
            hist, _ = np.histogram(dists, bins=bin_edges, density=True)
            rdf_results.append(hist)

        # Append the RDF results of the current trajectory to the list
        all_rdf_results.append(rdf_results)

        # Print a message indicating the trajectory has been processed
        print(f"Trajectory {traj_file} processed.")

    # Calculate the average RDF across all trajectories
    average_rdf = np.mean(np.concatenate(all_rdf_results), axis=0)

    # Save RDF data to a text file
    with open('average_rdf_data.txt', 'w') as rdf_file:
        rdf_file.write('Distance(Å)\tAverage RDF\n')
        for distance, rdf_value in zip(bin_edges[:-1], average_rdf):
            rdf_file.write(f'{distance}\t{rdf_value}\n')

    # Plot the average RDF for all trajectories
    plt.figure(figsize=(10, 8))
    plt.plot(bin_edges[:-1], average_rdf, label='O---H', linewidth=5, color='blue')  # Use a bold line for the average

    plt.xlabel('Distance (Å)', fontsize=32)
    plt.ylabel('RDF', fontsize=32)
    plt.legend(loc='upper right', fontsize=24, frameon=False)
    set_plot_style()

    # Save the plot as a PDF
    plt.savefig('average_rdf_plot.pdf')

    plt.show()