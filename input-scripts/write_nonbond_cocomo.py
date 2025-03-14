import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import csv
import math
    
# Define the amino acid order for indexing
DESIRED_ORDER = [
    'TRP', 'TYR', 'PHE', 'MET', 'LEU', 'ILE', 'VAL', 'ALA', 'PRO', 'GLY',
    'CYS', 'GLN', 'ASN', 'THR', 'SER', 'GLU', 'ASP', 'LYS', 'HIS', 'ARG'
]
    
# Load data from COCOMO repulsion file
def load_repulsion_data(filename):
    amino_acids, a_values, a0_values, radi_values = [], [], [], []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            amino_acids.append(parts[0])
            a_values.append(float(parts[2]))
            a0_values.append(float(parts[3]))
            radi_values.append(float(parts[4]))
    return amino_acids, a_values, a0_values, radi_values

# Load data from COCOMO attraction file
def load_attraction_data(filename):
    amino_acids, en_values = [], []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            amino_acids.append(parts[0])
            en_values.append(float(parts[1]))
    return amino_acids, en_values
    
# Load data
amino_acids_rep, a_values, a0_values, radi_values = load_repulsion_data('cocomo_rep.txt')
amino_acids_att, en_values = load_attraction_data('cocomo_att.txt')

# Create a mapping for amino acid indices
amino_acid_indices = {amino: i for i, amino in enumerate(DESIRED_ORDER)}

# Initialize interaction matrices
size = len(DESIRED_ORDER)
D_matrix = np.zeros((size, size))
alpha_matrix = np.zeros((size, size))
r0_matrix = np.zeros((size, size))
A_matrix = np.zeros((size, size))

# Define special amino acid pairs for D modification
SPECIAL_PAIRS = {('ARG', 'PHE'), ('ARG', 'TYR'), ('ARG', 'TRP'),
                 ('LYS', 'PHE'), ('LYS', 'TYR'), ('LYS', 'TRP')}

# Compute interaction parameters
for i in range(size):
    for j in range(size):
        aa_i, aa_j = DESIRED_ORDER[i], DESIRED_ORDER[size - 1 - j]
        
        # Retrieve properties from loaded data
        idx_i, idx_j = amino_acids_rep.index(aa_i), amino_acids_rep.index(aa_j)
        r_i, r_j = radi_values[idx_i], radi_values[idx_j]
        Ai_i, Ai_j = a_values[idx_i], a_values[idx_j]
        A0_i, A0_j = a0_values[idx_i], a0_values[idx_j]
        en_i, en_j = en_values[idx_i], en_values[idx_j]

        # Compute interaction parameters
        D = np.sqrt(en_i * en_j) * 0.239
        if (aa_i, aa_j) in SPECIAL_PAIRS or (aa_j, aa_i) in SPECIAL_PAIRS:
            D += 0.3 * 0.239
        
        r0 = (2 * r_i * 2**(1/5) + 2 * r_j * 2**(1/5)) / 2
        alpha = 5 / r0
        A = (Ai_i * Ai_j + (A0_i + A0_j)) * 1390 * 0.2388 / (4 * math.pi * 80)

        # Store values in matrices
        D_matrix[i, j] = D
        alpha_matrix[i, j] = alpha
        r0_matrix[i, j] = r0
        A_matrix[i, j] = A

# Save interaction parameters to file
with open('nonbond_COCOMO.txt', 'w') as f:
    for i in range(size):
        for j in range(size - i):
            f.write(f"pair_coeff\t{i+1}\t{20-j}\tmorse\t{D_matrix[i,j]:.2f}\t{alpha_matrix[i,j]:.2f}\t{r0_matrix[i,j]:.2f}\t25\n")
            f.write(f"pair_coeff\t{i+1}\t{20-j}\tyukawa\t{A_matrix[i,j]:.4f}\t40\n")


# Define function to plot a triangular matrix
def plot_triangular_matrix(matrix, vmin, vmax, label, cmap=cm.jet):
    mask = np.zeros_like(matrix, dtype=bool)
    for i in range(size):
        for j in range(size):
            if i + j <= size - 1:
                mask[i, j] = True
    masked_matrix = np.ma.array(matrix, mask=~mask)
    
    plt.imshow(masked_matrix, cmap=cmap, vmin=vmin, vmax=vmax)
    cbar = plt.colorbar()
    cbar.set_label(label, fontsize=18)
    cbar.ax.tick_params(labelsize=16)
    
    # Map amino acid names to single-letter codes
    aa_mapping = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                  'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                  'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                  'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
    single_letter_order = [aa_mapping[aa] for aa in DESIRED_ORDER]
    
    plt.xticks(np.arange(size), single_letter_order[::-1], fontsize=16)
    plt.yticks(np.arange(size), single_letter_order, fontsize=16)
    plt.show()

# Plot matrices
plot_triangular_matrix(D_matrix, vmin=0, vmax=0.2, label='kcal/mol')
plot_triangular_matrix(alpha_matrix, vmin=0.3, vmax=1.3, label='1/$\AA$')
plot_triangular_matrix(r0_matrix, vmin=5.5, vmax=16, label='$\AA$')
