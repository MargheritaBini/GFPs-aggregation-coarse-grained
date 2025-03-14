import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import csv
import math
    
# Define the desired order of amino acids
desired_order = ['TRP', 'TYR', 'PHE', 'MET', 'LEU', 'ILE', 'VAL', 'ALA', 'PRO', 'GLY', 'CYS', 'GLN', 'ASN', 'THR', 'SER', 'GLU', 'ASP', 'LYS', 'HIS', 'ARG']

# Load radii values from a text file
with open('calvados_radius.txt', 'r') as f:
    data = f.readlines()

amino_acids_radii = []
radii_values = []

for line in data:
    parts = line.strip().split()
    amino_acids_radii.append(parts[0])
    radii_values.append(float(parts[1]))

# Load lambda and charge values from CSV file
csv_file_path = 'calvados.csv'
amino_acids_lambda = []
lambda_values = []
charge_values = []

with open(csv_file_path, mode='r') as csvfile:
    csvreader = csv.DictReader(csvfile)
    for row in csvreader:
        amino_acids_lambda.append(row['three'])
        lambda_values.append(float(row['CALVADOS2']))
        charge_values.append(float(row['q']))


# Create a dictionary for amino acids
amino_acid_indices = {amino: i for i, amino in enumerate(desired_order)}

# Initialize matrices for different parameters
D_matrix = np.zeros((len(desired_order), len(desired_order)))
alpha_matrix = np.zeros((len(desired_order), len(desired_order)))
r0_matrix = np.zeros((len(desired_order), len(desired_order)))
A_matrix = np.zeros((len(desired_order), len(desired_order)))

# Compute interaction parameters and store in matrices
for i in range(len(desired_order)):
    for j in range(len(desired_order)):
        lambda_i = lambda_values[amino_acids_lambda.index(desired_order[i])]
        sigma_i = radii_values[amino_acids_radii.index(desired_order[i])]
        charge_i = charge_values[amino_acids_lambda.index(desired_order[i])]

        lambda_j = lambda_values[amino_acids_lambda.index(desired_order[19 - j])]
        sigma_j = radii_values[amino_acids_radii.index(desired_order[19 - j])]
        charge_j = charge_values[amino_acids_lambda.index(desired_order[19 - j])]
        
        rcut = 24
        D = (lambda_i + lambda_j) / 2 * 0.8368 / 4.184 + (lambda_i + lambda_j) / 2 * 0.8368 / 4.184 * (((sigma_i + sigma_j) / 2 / rcut) ** 12 - ((sigma_i + sigma_j) / 2 / rcut) ** 6)
        alpha = (9 - 3 * (lambda_i + lambda_j) / 2) / (((sigma_i + sigma_j) / 2) * 2 ** (1 / 6))
        r0 = ((sigma_i + sigma_j) / 2) * 2 ** (1 / 6)
        A = charge_i * charge_j * 1390 * 0.2388 / (4 * math.pi * 80)
        
        D_matrix[i, j] = D
        alpha_matrix[i, j] = alpha
        r0_matrix[i, j] = r0
        A_matrix[i, j] = A

# Save interaction parameters to a text file
with open('nonbond_CALVADOS.txt', 'w') as txtfile:
    for i in range(len(desired_order)):
        for j in range(len(desired_order) - i):
            txtfile.write(f"pair_coeff\t{i+1}\t{20-j}\tmorse\t{D_matrix[i,j]:.2f}\t{alpha_matrix[i,j]:.2f}\t{r0_matrix[i,j]:.2f}\t24\n")
            txtfile.write(f"pair_coeff\t{i+1}\t{20-j}\tyukawa\t{A_matrix[i,j]:.4f}\t40\n")

# Function to plot a matrix

def plot_matrix(matrix, vmin, vmax, label):
    n = matrix.shape[0]
    mask = np.zeros_like(matrix, dtype=bool)
    for i in range(n):
        for j in range(n):
            if i + j <= n - 1:
                mask[i, j] = True
    masked_matrix = np.ma.array(matrix, mask=~mask)
    cmap = cm.jet
    
    plt.imshow(masked_matrix, cmap=cmap, vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(label=label)
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label(label, fontsize=18)
    
    single_letter_order = [aa_mapping[aa] for aa in desired_order]
    plt.xticks(np.arange(len(single_letter_order)), single_letter_order[::-1], fontsize=16)
    plt.yticks(np.arange(len(single_letter_order)), single_letter_order, fontsize=16)
    plt.show()

# Amino acid mapping
aa_mapping = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

# Plot all matrices
plot_matrix(D_matrix, vmin=0, vmax=0.2, label='kcal/mol')
plot_matrix(alpha_matrix, vmin=0.3, vmax=1.3, label='1/$\AA$')
plot_matrix(r0_matrix, vmin=5.5, vmax=16, label='$\AA$')
plot_matrix(A_matrix, vmin=-0.3, vmax=0.3, label='kcal/mol')
