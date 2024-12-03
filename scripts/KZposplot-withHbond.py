import matplotlib.pyplot as plt
from ase.io import read
import re
import os
import numpy as np

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

# Initialize cumulative time
cumulative_time = 0.0

# Flag to check if K atom indices have been initialized
k_indices = None

# List of water molecule indices
water_indices = [4462, 4594, 4216, 4213]

# Initialize z_positions for water molecules
water_z_positions = {i: [] for i in water_indices}

# Indices of hydrogen atoms to check proximity
h_indices = [576, 1630, 2683, 3735]

# Distance threshold in Å
distance_threshold = 2.0

# Initialize a list to store proximity status for each water molecule
water_proximity_status = {i: [] for i in water_indices}

# Open a text file to write positions
with open('z_positions_all.txt', 'w') as f_out:
    # Header for the output file
    header_written = False

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
                time = 0.0  # Default to 0 if time is missing

            # Adjust time with cumulative time
            adjusted_time = (time + cumulative_time) / 1000.0
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

            # Check proximity of water molecules to hydrogen atoms
            for water_index in water_indices:
                is_close = False
                water_pos = frame.positions[water_index]

                for h_index in h_indices:
                    h_pos = frame.positions[h_index]
                    distance = np.linalg.norm(water_pos - h_pos)
                    if distance <= distance_threshold:
                        is_close = True
                        break

                water_proximity_status[water_index].append(is_close)

            # Write to output file
            if header_written:
                line = f"{adjusted_time}\t" + '\t'.join([f"{z}" for z in z_pos]) + '\t' + '\t'.join([f"{z}" for z in water_z_pos]) + '\n'
                f_out.write(line)

        # Update cumulative time after processing each trajectory
        if times:
            cumulative_time = times[-1] * 1000.0  # Convert back to fs for the next trajectory

        # Print completion message for each trajectory
        print(f"Trajectory {traj_file} completed")

bright_colors = ['xkcd:dark green', 'xkcd:light blue']

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

# Check if data was collected
if not times or not any(z_positions.values()) or not any(water_z_positions.values()):
    print("No data was collected. Please check the trajectory files and script.")
else:
    plt.figure(figsize=(10, 8))
    
    # Plot z-position vs. time for all K⁺ ions
    for i in range(len(k_indices)):
        label = 'K⁺' if i == 0 else None  # Label only the first K⁺
        plt.plot(times, z_positions[i], label=label, linewidth=5, color=bright_colors[0])

    # Plot z-position vs. time for all specified water molecules
    for j, water_index in enumerate(water_indices):
        proximity = water_proximity_status[water_index]
        start_idx = 0

        for k in range(1, len(times)):
            # Check if proximity status changes or if it's the last point
            if proximity[k] != proximity[k - 1] or k == len(times) - 1:
                # Determine the color based on proximity
                color = 'xkcd:dark blue' if proximity[k - 1] else bright_colors[1]
                
                # Plot the segment
                plt.plot(times[start_idx:k], water_z_positions[water_index][start_idx:k], 
                         color=color, label='Water' if j == 0 and start_idx == 0 else "", linewidth=5)
                
                # Update the start index for the next segment
                start_idx = k

    plt.xlabel('Time (ps)', fontsize=32)
    plt.ylabel('Z position (Å)', fontsize=32)
    plt.ylim(-15, 25)  # Set y-axis limit to 15
    plt.legend(loc='upper right', fontsize=24, frameon=False)
    set_plot_style()
    
    # Save the plot as a PDF
    plt.savefig('z_positions_plot.pdf', format='pdf')
    
    plt.show()