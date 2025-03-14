import numpy as np
import matplotlib.pyplot as plt

def parse_pdb_frames(file_path):
    """
    Parse a PDB file containing multiple frames separated by 'END' lines.
    Returns a list of NumPy arrays, each containing atomic coordinates for a frame.
    """
    frames = []
    current_coords = []
    
    with open(file_path, 'r') as file:
        for line_num, line in enumerate(file, 1):
            if line.startswith("ATOM"):
                parts = line.split()
                try:
                    x, y, z = float(parts[6]), float(parts[7]), float(parts[8])
                    current_coords.append([x, y, z])
            elif line.startswith("END"):
                if current_coords:
                    frames.append(np.array(current_coords))
                    current_coords = []
    
    if current_coords:
        frames.append(np.array(current_coords))
    
    return frames

def compute_distances(coords1, coords2):
    """
    Compute pairwise Euclidean distances between two sets of atomic coordinates.
    """
    diff = coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :]
    distances = np.sqrt(np.sum(diff**2, axis=2))
    return distances

def update_contact_matrix(coords, num_atoms_per_protein, cutoff):
    """
    Update the contact matrix based on distances less than the cutoff.
    """
    num_proteins = len(coords) // num_atoms_per_protein
    matrix = np.zeros((num_atoms_per_protein, num_atoms_per_protein))

    for i in range(num_proteins):
        for j in range(num_proteins):
            if i != j:
                start_i, end_i = i * num_atoms_per_protein, (i + 1) * num_atoms_per_protein
                start_j, end_j = j * num_atoms_per_protein, (j + 1) * num_atoms_per_protein
                coords_i = coords[start_i:end_i]
                coords_j = coords[start_j:end_j]
                distances = compute_distances(coords_i, coords_j)
                contacts = distances < cutoff
                
                # Update the matrix using a weighted contribution:
                # This approach sums the number of contacts for each atom in i and each atom in j,
                # then multiplies these sums to get a combined interaction score.
                matrix += contacts.sum(axis=1).reshape(-1, 1) * contacts.sum(axis=0).reshape(1, -1)

                # Alternative approach:
                # Instead of weighting contacts, update the matrix directly where contacts exist.
                # This simply adds 1 to matrix[i, j] whenever atom i and atom j are in contact.
                
                #matrix += contacts

    return matrix

def normalize_matrix(matrix):
    """
    Normalize the matrix to get probabilities.
    """
    total_contacts = np.sum(matrix)
    if total_contacts > 0:
        return matrix / total_contacts
    else:
        return matrix


def plot_contact_matrix(matrix, title):
    """
    Plot the contact matrix as a heatmap.
    """
    plt.figure(figsize=(8, 6))
    plt.imshow(matrix, cmap='viridis', interpolation='nearest', vmin=0, vmax=0.0018)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.title(title, fontsize=16)
    plt.colorbar(label='Contact Probability', orientation='vertical')
    plt.show()

# Main script
if __name__ == "__main__":
    pdb_file_path = 'traj_20.pdb'       # Path to the PDB file containing multiple frames
    num_atoms_per_protein = 229           # Number of atoms per protein in each frame
    cutoff_distance = 8.2                 # Distance threshold in angstroms

    # Parse all frames from the PDB file
    frames = parse_pdb_frames(pdb_file_path)
    num_frames = len(frames)
    print(f"Total number of frames found: {num_frames}\n")

    # Initialize a matrix to accumulate contact matrices
    accumulated_matrix = np.zeros((num_atoms_per_protein, num_atoms_per_protein))

    for frame_index, coords in enumerate(frames, 1):
        # Update the contact matrix for the current frame
        contact_matrix = update_contact_matrix(coords, num_atoms_per_protein, cutoff_distance)

        # Accumulate the contact matrices
        accumulated_matrix += contact_matrix

    # Compute the average contact matrix
    average_matrix = accumulated_matrix / num_frames

    # Normalize the averaged contact matrix to get probabilities
    average_prob_matrix = normalize_matrix(average_matrix)

    # Plot the averaged contact probability matrix
    plot_contact_matrix(average_prob_matrix, 'Averaged Contact Probability Matrix')



 
