import numpy as np
import matplotlib.pyplot as plt

# Define the simulation box length and half its size for periodic boundary condition corrections
BOX_LENGTH = 500.0  # Length of the simulation box (assumed cubic)
HALF_BOX_LENGTH = BOX_LENGTH / 2  # Half-box length for PBC handling

def parse_pdb_frames(file_path):
    """
    Reads a PDB trajectory file and extracts atomic positions for each frame.
    
    Parameters:
    - file_path (str): Path to the PDB file.
    
    Returns:
    - frames (list of np.ndarray): List of frames, where each frame is an array of atomic positions (N x 3).
    """
    frames = []
    current_frame = []

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):  # Extract atomic coordinates from 'ATOM' lines
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                current_frame.append([x, y, z])
            elif line.startswith("END"):  # Identify end of a frame and store it
                frames.append(np.array(current_frame))
                current_frame = []
    
    return frames

def compute_msd(atom_positions):
    """
    Computes the mean squared displacement (MSD) between two atomic positions,
    applying periodic boundary corrections if necessary.
    
    Parameters:
    - atom_positions (np.ndarray): Array containing two atomic positions (2 x 3).
    
    Returns:
    - msd (float): The squared displacement between the two positions.
    """
    corrected_positions = np.copy(atom_positions)
    
    # Apply periodic boundary corrections
    for j in range(3):  # Loop over x, y, z coordinates
        delta = atom_positions[0, j] - atom_positions[1, j]
        if delta > HALF_BOX_LENGTH:
            corrected_positions[0, j] -= BOX_LENGTH
        elif delta < -HALF_BOX_LENGTH:
            corrected_positions[0, j] += BOX_LENGTH
    
    # Compute squared displacement
    msd_value = np.sum((corrected_positions[0] - corrected_positions[1]) ** 2)
    
    return msd_value

def calculate_msd_rolling_average(frames):
    """
    Computes the MSD over increasing time intervals and calculates propagated errors.
    
    Parameters:
    - frames (list of np.ndarray): List of atomic position frames.
    
    Returns:
    - msd_all_frames (np.ndarray): MSD values for each atom at different time intervals.
    - avg_msd (np.ndarray): Average MSD over all atoms at each time interval.
    - propagated_error_per_time (np.ndarray): Error in MSD calculations.
    """
    n_frames = len(frames)
    n_atoms = len(frames[0])

    msd_all_frames = []  # MSD for each atom over time intervals
    avg_msd_per_time = []  # Average MSD per time interval
    propagated_error_per_time = []  # Error propagation per time interval

    for t in range(1, n_frames):  # Time intervals from 1 to n_frames-1
        msd_per_time_per_atom = []

        for atom_index in range(n_atoms):
            msd_intervals = []
            
            for start in range(n_frames - t):  # Compute MSD at intervals of length `t`
                ref = start
                target = start + t
                atom_positions = np.array([frames[target][atom_index], frames[ref][atom_index]])
                msd_value = compute_msd(atom_positions)
                msd_intervals.append(msd_value)

            # Compute mean and standard error for the MSD at this atom and time interval
            mean_msd_atom = np.mean(msd_intervals)
            error_msd_atom = np.std(msd_intervals) / np.sqrt(len(msd_intervals))
            msd_per_time_per_atom.append((mean_msd_atom, error_msd_atom))

        # Compute average MSD across all atoms at this time interval
        mean_msds = [msd[0] for msd in msd_per_time_per_atom]
        errors = [msd[1] for msd in msd_per_time_per_atom]
        avg_msd_per_time.append(np.mean(mean_msds))

        # Calculate propagated error
        squared_diffs = [(mean_msd - np.mean(mean_msds)) ** 2 for mean_msd in mean_msds]
        variance = np.mean(squared_diffs)
        propagated_error = np.sqrt(variance / n_atoms)

        propagated_error_per_time.append(propagated_error)

        msd_all_frames.append(mean_msds)

        print(f"Completed computation for time interval {t}")

    return np.array(msd_all_frames).T, np.array(avg_msd_per_time), np.array(propagated_error_per_time)

def plot_msd(msd_all_frames, avg_msd, error_per_time):
    """
    Plots MSD over time, including individual atom MSDs and the average MSD with error bars.
    
    Parameters:
    - msd_all_frames (np.ndarray): MSD values for each atom over time.
    - avg_msd (np.ndarray): Average MSD across all atoms.
    - error_per_time (np.ndarray): Propagated error in MSD values.
    """
    time = np.arange(1, avg_msd.size + 1)

    # Plot individual msds for each atom
    for atom_msd in msd_all_frames:
        plt.plot(time, atom_msd, color='gray', alpha=0.5, linewidth=0.5)

    # Plot the average msd with error bars
    plt.errorbar(time, avg_msd, yerr=error_per_time, color='blue', label='Average msd', linewidth=2, capsize=3)

    # Save average msd to file
    output_file_path = 'msd_calvados_sc.txt'
    with open(output_file_path, 'w') as output_file:
        for i in range(len(time)):
            output_file.write(f"{time[i]}\t{avg_msd[i]}\t{error_per_time[i]}\n")

    # Customize the plot appearance
    plt.xlabel("Time (frames)")
    plt.ylabel("msd (Å)")
    plt.title("msd for Each Atom and Average msd Over Time with Error")
    plt.legend()
    plt.show()

# File path to your PDB file
pdb_file_path = 'traj_cm.pdb'

# Parse PDB file and calculate msd
frames = parse_pdb_frames(pdb_file_path)
msd_all_frames, avg_msd, error_per_time = calculate_msd_rolling_average(frames)

# Plot msd results
plot_msd(msd_all_frames, avg_msd, error_per_time)
