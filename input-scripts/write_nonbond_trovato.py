import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import csv
import math
    
# Define the desired order of amino acids
desired_order = ['TRP', 'TYR', 'PHE', 'MET', 'LEU', 'ILE', 'VAL', 'ALA', 'PRO', 'GLY',
                 'CYS', 'GLN', 'ASN', 'THR', 'SER', 'GLU', 'ASP', 'LYS', 'HIS', 'ARG']
# Initialize 20x20 matrices (for energy, alpha, r, and A) with zeros
matrix_size = 20
energy_matrix = np.zeros((matrix_size, matrix_size))
alpha_matrix = np.zeros((matrix_size, matrix_size))
r_matrix = np.zeros((matrix_size, matrix_size))
A_matrix = np.zeros((matrix_size, matrix_size))

# h1-h1 interaction parameters and update matrices for h1-h1 ((positions 4 to 7 and position 10))
e2 = 0.398 - 0.11
a1 = 0.761
a2 = 0.217
r2 = 10.96
r1 = 9.36

for i in range(4, 8):
    for j in range(4, 8):
        energy_matrix[i][19 - j] = e2
        alpha_matrix[i][19 - j] = (9 * a2 + a1) / 10
        r_matrix[i][19 - j] = (r2 + 2 * r1) / 3
        A_matrix[i][19 - j] = 0
        
        energy_matrix[10][19 - j] = e2
        alpha_matrix[10][19 - j] = (9 * a2 + a1) / 10
        r_matrix[10][19 - j] = (r2 + 2 * r1) / 3
        A_matrix[10][19 - j] = 0

    energy_matrix[i][19 - 10] = e2
    alpha_matrix[i][19 - 10] = (9 * a2 + a1) / 10
    r_matrix[i][19 - 10] = (r2 + 2 * r1) / 3
    A_matrix[i][19 - 10] = 0
    
    energy_matrix[10][19 - 10] = e2
    alpha_matrix[10][19 - 10] = (9 * a2 + a1) / 10
    r_matrix[10][19 - 10] = (r2 + 2 * r1) / 3
    A_matrix[10][19 - 10] = 0

# h1-h2 interaction parameters. Update matrices for h1-h2 interaction (positions 4 to 7 and positions 0 to 3)
e2 = 0.184
a1 = 0.719
a2 = 0.278
r2 = 11.35
r1 = 9.50

for i in range(4, 8):
    for j in range(0, 4):
        energy_matrix[i][19 - j] = e2
        energy_matrix[10][19 - j] = e2
        alpha_matrix[i][19 - j] = (9 * a2 + a1) / 10
        alpha_matrix[10][19 - j] = (9 * a2 + a1) / 10
        r_matrix[i][19 - j] = (r2 + 2 * r1) / 3
        r_matrix[10][19 - j] = (r2 + 2 * r1) / 3
        A_matrix[i][19 - j] = 0
        A_matrix[10][19 - j] = 0

    energy_matrix[i][19 - 8] = e2
    energy_matrix[i][19 - 9] = e2
    energy_matrix[10][19 - 8] = e2
    energy_matrix[10][19 - 9] = e2
    
    alpha_matrix[i][19 - 8] = (9 * a2 + a1) / 10
    alpha_matrix[i][19 - 9] = (9 * a2 + a1) / 10
    alpha_matrix[10][19 - 8] = (9 * a2 + a1) / 10
    alpha_matrix[10][19 - 9] = (9 * a2 + a1) / 10
    
    r_matrix[i][19 - 8] = (r2 + 2 * r1) / 3
    r_matrix[i][19 - 9] = (r2 + 2 * r1) / 3
    r_matrix[10][19 - 8] = (r2 + 2 * r1) / 3
    r_matrix[10][19 - 9] = (r2 + 2 * r1) / 3
    
    A_matrix[i][19 - 8] = 0
    A_matrix[i][19 - 9] = 0
    A_matrix[10][19 - 8] = 0
    A_matrix[10][19 - 9] = 0

 
# h1-p1 interaction parameters
e1 = 0.9824
e2 = 0.115
a1 = 0.357
a2 = 0.260
r2 = 12.24
r1 = 6.54

# Update matrices for h1-p1 interaction (positions 4 to 7 and positions 12 to 14)
for i in range(4, 8):
    for j in range(12, 15):
        energy_matrix[i][19 - j] = e2
        energy_matrix[10][19 - j] = e2
        alpha_matrix[i][19 - j] = a2
        alpha_matrix[10][19 - j] = a2
        r_matrix[i][19 - j] = (19 * r2 + r1) / 20
        r_matrix[10][19 - j] = (19 * r2 + r1) / 20
        A_matrix[i][19 - j] = 0
        A_matrix[10][19 - j] = 0

    energy_matrix[i][19 - 18] = e2
    energy_matrix[10][19 - 18] = e2
    alpha_matrix[i][19 - 18] = a2
    alpha_matrix[10][19 - 18] = a2
    r_matrix[i][19 - 18] = (19 * r2 + r1) / 20
    r_matrix[10][19 - 18] = (19 * r2 + r1) / 20
    A_matrix[i][19 - 18] = 0
    A_matrix[10][19 - 18] = 0

    
# h1-p2 interaction parameters
e2 = 0.082
a2 = 0.288
a1 = 0.277
r2 = 12.42
r1 = 6.55
e1 = 2.11

# Update matrices for h1-p2 interaction (positions 4 to 7 and positions 15 to 17)
for i in range(4, 8):
    for j in range(15, 18):
        energy_matrix[i][19 - j] = e2
        energy_matrix[10][19 - j] = e2
        alpha_matrix[i][19 - j] = a2
        alpha_matrix[10][19 - j] = a2
        r_matrix[i][19 - j] = (19 * r2 + r1) / 20
        r_matrix[10][19 - j] = (19 * r2 + r1) / 20
        A_matrix[i][19 - j] = 0
        A_matrix[10][19 - j] = 0

    energy_matrix[i][19 - 11] = e2
    energy_matrix[i][19 - 19] = e2
    energy_matrix[10][19 - 11] = e2
    energy_matrix[10][19 - 19] = e2
    alpha_matrix[i][19 - 11] = a2
    alpha_matrix[i][19 - 19] = a2
    alpha_matrix[10][19 - 11] = a2
    alpha_matrix[10][19 - 19] = a2
    r_matrix[i][19 - 11] = (19 * r2 + r1) / 20
    r_matrix[i][19 - 19] = (19 * r2 + r1) / 20
    r_matrix[10][19 - 11] = (19 * r2 + r1) / 20
    r_matrix[10][19 - 19] = (19 * r2 + r1) / 20
    A_matrix[i][19 - 11] = 0
    A_matrix[i][19 - 19] = 0
    A_matrix[10][19 - 11] = 0
    A_matrix[10][19 - 19] = 0

    
# h1-h2 interaction parameters
e2 = 0.184
a1 = 0.719
a2 = 0.278
r2 = 11.35
r1 = 9.50

# Update matrices for h1-h2 interaction (positions 0 to 3 and 4 to 7)
for i in range(0, 4):
    for j in range(4, 8):
        energy_matrix[i][19 - j] = e2
        energy_matrix[8][19 - j] = e2
        energy_matrix[9][19 - j] = e2
        alpha_matrix[i][19 - j] = (9 * a2 + a1) / 10
        alpha_matrix[8][19 - j] = (9 * a2 + a1) / 10
        alpha_matrix[9][19 - j] = (9 * a2 + a1) / 10
        r_matrix[i][19 - j] = (r2 + 2 * r1) / 3
        r_matrix[8][19 - j] = (r2 + 2 * r1) / 3
        r_matrix[9][19 - j] = (r2 + 2 * r1) / 3
        A_matrix[i][19 - j] = 0
        A_matrix[8][19 - j] = 0
        A_matrix[9][19 - j] = 0

    # Update for position 10
    energy_matrix[i][19 - 10] = e2
    energy_matrix[8][19 - 10] = e2
    energy_matrix[9][19 - 10] = e2
    alpha_matrix[i][19 - 10] = (9 * a2 + a1) / 10
    alpha_matrix[8][19 - 10] = (9 * a2 + a1) / 10
    alpha_matrix[9][19 - 10] = (9 * a2 + a1) / 10
    r_matrix[i][19 - 10] = (r2 + 2 * r1) / 3
    r_matrix[8][19 - 10] = (r2 + 2 * r1) / 3
    r_matrix[9][19 - 10] = (r2 + 2 * r1) / 3
    A_matrix[i][19 - 10] = 0
    A_matrix[8][19 - 10] = 0
    A_matrix[9][19 - 10] = 0

# h2-h2 interaction parameters
e2 = 0.109
a1 = 0.539
a2 = 0.285
r2 = 11.56
r1 = 6.73

# Update matrices for h2-h2 interaction (positions 0 to 3 and 0 to 3)
for i in range(0, 4):
    for j in range(0, 4):
        energy_matrix[i][19 - j] = e2
        energy_matrix[8][19 - j] = e2
        energy_matrix[9][19 - j] = e2
        alpha_matrix[i][19 - j] = a2
        alpha_matrix[8][19 - j] = a2
        alpha_matrix[9][19 - j] = a2
        r_matrix[i][19 - j] = (r2 * 19 + r1) / 20
        r_matrix[8][19 - j] = (r2 * 19 + r1) / 20
        r_matrix[9][19 - j] = (r2 * 19 + r1) / 20
        A_matrix[i][19 - j] = 0
        A_matrix[8][19 - j] = 0
        A_matrix[9][19 - j] = 0

    # Specific positions 8 and 9
    energy_matrix[i][19 - 8] = e2
    energy_matrix[i][19 - 9] = e2
    energy_matrix[8][19 - 8] = e2
    energy_matrix[9][19 - 8] = e2
    energy_matrix[8][19 - 9] = e2
    energy_matrix[9][19 - 9] = e2
    alpha_matrix[i][19 - 8] = a2
    alpha_matrix[i][19 - 9] = a2
    alpha_matrix[8][19 - 8] = a2
    alpha_matrix[9][19 - 8] = a2
    alpha_matrix[8][19 - 9] = a2
    alpha_matrix[9][19 - 9] = a2
    r_matrix[i][19 - 8] = (r2 * 19 + r1) / 20
    r_matrix[i][19 - 9] = (r2 * 19 + r1) / 20
    r_matrix[8][19 - 8] = (r2 * 19 + r1) / 20
    r_matrix[9][19 - 8] = (r2 * 19 + r1) / 20
    r_matrix[8][19 - 9] = (r2 * 19 + r1) / 20
    r_matrix[9][19 - 9] = (r2 * 19 + r1) / 20
    A_matrix[i][19 - 8] = 0
    A_matrix[i][19 - 9] = 0
    A_matrix[8][19 - 8] = 0
    A_matrix[9][19 - 8] = 0
    A_matrix[8][19 - 9] = 0
    A_matrix[9][19 - 9] = 0

# h2p1 interaction parameters
e2 = 0.048
a2 = 0.314
a1 = 0.321
r2 = 12.4
r1 = 6.24
e1 = 1.5615

# Update matrices for h2p1 interaction (positions 0 to 3 and 12 to 15)
for i in range(0, 4):
    for j in range(12, 15):
        energy_matrix[i][19 - j] = e2
        energy_matrix[8][19 - j] = e2
        energy_matrix[9][19 - j] = e2
        alpha_matrix[i][19 - j] = a2
        alpha_matrix[8][19 - j] = a2
        alpha_matrix[9][19 - j] = a2
        r_matrix[i][19 - j] = (19 * r2 + r1) / 20
        r_matrix[8][19 - j] = (19 * r2 + r1) / 20
        r_matrix[9][19 - j] = (19 * r2 + r1) / 20
        A_matrix[i][19 - j] = 0
        A_matrix[8][19 - j] = 0
        A_matrix[9][19 - j] = 0

    # Specific position 18
    energy_matrix[i][19 - 18] = e2
    energy_matrix[8][19 - 18] = e2
    energy_matrix[9][19 - 18] = e2
    alpha_matrix[i][19 - 18] = a2
    alpha_matrix[8][19 - 18] = a2
    alpha_matrix[9][19 - 18] = a2
    r_matrix[i][19 - 18] = (19 * r2 + r1) / 20
    r_matrix[8][19 - 18] = (19 * r2 + r1) / 20
    r_matrix[9][19 - 18] = (19 * r2 + r1) / 20
    A_matrix[i][19 - 18] = 0
    A_matrix[8][19 - 18] = 0
    A_matrix[9][19 - 18] = 0

# h2p2 interaction parameters
e2 = 0.02
a2 = 0.378
a1 = 0.444
r2 = 12.8
r1 = 6.41
e1 = 0.5198

# Update matrices for h2p2 interaction (positions 0 to 3 and 15 to 18)
for i in range(0, 4):
    for j in range(15, 18):
        energy_matrix[i][19 - j] = e2
        energy_matrix[8][19 - j] = e2
        energy_matrix[9][19 - j] = e2
        alpha_matrix[i][19 - j] = a2
        alpha_matrix[8][19 - j] = a2
        alpha_matrix[9][19 - j] = a2
        r_matrix[i][19 - j] = (4 * r2 + r1) / 5
        r_matrix[8][19 - j] = (4 * r2 + r1) / 5
        r_matrix[9][19 - j] = (4 * r2 + r1) / 5
        A_matrix[i][19 - j] = 0
        A_matrix[8][19 - j] = 0
        A_matrix[9][19 - j] = 0

    # Specific positions 11 and 19
    for pos in [11, 19]:
        energy_matrix[i][19 - pos] = e2
        energy_matrix[8][19 - pos] = e2
        energy_matrix[9][19 - pos] = e2
        alpha_matrix[i][19 - pos] = a2
        alpha_matrix[8][19 - pos] = a2
        alpha_matrix[9][19 - pos] = a2
        r_matrix[i][19 - pos] = (4 * r2 + r1) / 5
        r_matrix[8][19 - pos] = (4 * r2 + r1) / 5
        r_matrix[9][19 - pos] = (4 * r2 + r1) / 5
        A_matrix[i][19 - pos] = 0
        A_matrix[8][19 - pos] = 0
        A_matrix[9][19 - pos] = 0

# p1h1 interaction parameters
e2 = 0.115
a1 = 0.357
a2 = 0.260
r2 = 12.24
r1 = 6.54
e1 = 0.9824

# Update matrices for p1h1 interaction (positions 12 to 14 and 4 to 7)
for i in range(12, 15):
    for j in range(4, 8):
        energy_matrix[i][19 - j] = e2
        energy_matrix[18][19 - j] = e2
        alpha_matrix[i][19 - j] = a2
        alpha_matrix[18][19 - j] = a2
        r_matrix[i][19 - j] = (19 * r2 + r1) / 20
        r_matrix[18][19 - j] = (19 * r2 + r1) / 20
        A_matrix[i][19 - j] = 0
        A_matrix[18][19 - j] = 0

    # Specific position 10
    energy_matrix[i][19 - 10] = e2
    energy_matrix[18][19 - 10] = e2
    alpha_matrix[i][19 - 10] = a2
    alpha_matrix[18][19 - 10] = a2
    r_matrix[i][19 - 10] = (19 * r2 + r1) / 20
    r_matrix[18][19 - 10] = (19 * r2 + r1) / 20
    A_matrix[i][19 - 10] = 0
    A_matrix[18][19 - 10] = 0

# p1h2 interaction parameters
e2 = 0.048
a2 = 0.314
a1 = 0.321
r2 = 12.4
r1 = 6.24
e1 = 1.5615

# Update matrices for p1h2 interaction (positions 12 to 14 and 0 to 3)
for i in range(12, 15):
    for j in range(0, 4):
        energy_matrix[i][19 - j] = e2
        energy_matrix[18][19 - j] = e2
        alpha_matrix[i][19 - j] = a2
        alpha_matrix[18][19 - j] = a2
        r_matrix[i][19 - j] = (19 * r2 + r1) / 20
        r_matrix[18][19 - j] = (19 * r2 + r1) / 20
        A_matrix[i][19 - j] = 0
        A_matrix[18][19 - j] = 0

    # Specific positions 8 and 9
    for j in [8, 9]:
        energy_matrix[i][19 - j] = e2
        energy_matrix[18][19 - j] = e2
        alpha_matrix[i][19 - j] = a2
        alpha_matrix[18][19 - j] = a2
        r_matrix[i][19 - j] = (19 * r2 + r1) / 20
        r_matrix[18][19 - j] = (19 * r2 + r1) / 20
        A_matrix[i][19 - j] = 0
        A_matrix[18][19 - j] = 0

# p1-p1 interaction parameters
e2 = 0.011
a2 = 0.46
a1 = 0.364
r2 = 12.42
r1 = 6.26
e1 = 1.045
A = 0

# Update matrices for p1-p1 interaction (positions 12 to 14)
for i in range(12, 15):
    for j in range(12, 15):
        energy_matrix[i][19 - j] = e2
        energy_matrix[18][19 - j] = e2
        r_matrix[i][19 - j] = (14 * r2 + r1) / 15
        r_matrix[18][19 - j] = (14 * r2 + r1) / 15
        alpha_matrix[i][19 - j] = a1
        alpha_matrix[18][19 - j] = a1
        A_matrix[i][19 - j] = 0
        A_matrix[18][19 - j] = 0

    # Specific position for j = 18
    energy_matrix[i][19 - 18] = e2
    energy_matrix[18][19 - 18] = e2
    r_matrix[i][19 - 18] = (14 * r2 + r1) / 15
    r_matrix[18][19 - 18] = (14 * r2 + r1) / 15
    alpha_matrix[i][19 - 18] = a1
    alpha_matrix[18][19 - 18] = a1
    A_matrix[i][19 - 18] = 0
    A_matrix[18][19 - 18] = A

# p1-p2 interaction parameters
e2 = 3.71e-5
a1 = 0.662
a2 = 0.479
r2 = 17.80
r1 = 13.1
A_neg = 0
A_pos = 0

# Update matrices for p1-p2 interaction (positions 12 to 14 for i and 15 to 17 for j)
for i in range(12, 15):
    for j in range(15, 18):
        energy_matrix[i][19 - j] = e2
        alpha_matrix[i][19 - j] = a2
        r_matrix[i][19 - j] = (3 * r2 + r1) / 4
        A_matrix[i][19 - j] = 0

    # Specific positions for j = 15, 16, 17, 11, and 19
    for j in [15, 16, 17, 11, 19]:
        energy_matrix[18][19 - j] = e2
        alpha_matrix[18][19 - j] = a2
        r_matrix[18][19 - j] = (3 * r2 + r1) / 4
        A_matrix[18][19 - j] = A_neg if j in [15, 16] else A_pos

    # Specific positions for i = 12, 13, 14 and j = 11, 19
    for j in [11, 19]:
        energy_matrix[i][19 - j] = e2
        alpha_matrix[i][19 - j] = a2
        r_matrix[i][19 - j] = (3 * r2 + r1) / 4
        A_matrix[i][19 - j] = 0

    # Specific positions for i = 18 and j = 11, 19
    for j in [11, 19]:
        A_matrix[18][19 - j] = A_pos

# p2-h1 interaction parameters
e2 = 0.082
a2 = 0.288
a1 = 0.277
r2 = 12.42
r1 = 6.55
e1 = 2.11

# Update matrices for p2-h1 interaction (positions 15 to 17 for i and 4 to 7 for j)
for i in range(15, 18):
    for j in range(4, 8):
        energy_matrix[i][19 - j] = e2
        energy_matrix[11][19 - j] = e2
        energy_matrix[19][19 - j] = e2
        alpha_matrix[i][19 - j] = a2
        alpha_matrix[11][19 - j] = a2
        alpha_matrix[19][19 - j] = a2
        r_matrix[i][19 - j] = (r1 + 19 * r2) / 20
        r_matrix[11][19 - j] = (r1 + 19 * r2) / 20
        r_matrix[19][19 - j] = (r1 + 19 * r2) / 20
        A_matrix[i][19 - j] = 0
        A_matrix[11][19 - j] = 0
        A_matrix[19][19 - j] = 0

    # Specific position for j = 10
    for j in [10]:
        energy_matrix[i][19 - j] = e2
        energy_matrix[11][19 - j] = e2
        energy_matrix[19][19 - j] = e2
        alpha_matrix[i][19 - j] = a2
        alpha_matrix[11][19 - j] = a2
        alpha_matrix[19][19 - j] = a2
        r_matrix[i][19 - j] = (r1 + 19 * r2) / 20
        r_matrix[11][19 - j] = (r1 + 19 * r2) / 20
        r_matrix[19][19 - j] = (r1 + 19 * r2) / 20
        A_matrix[i][19 - j] = 0
        A_matrix[11][19 - j] = 0
        A_matrix[19][19 - j] = 0

# p2-h2 interaction parameters
e2 = 0.02
a2 = 0.378
a1 = 0.444
r2 = 12.8
r1 = 6.41
e1 = 0.5198

# Update matrices for p2-h2 interaction (positions 15 to 17 for i and 0 to 3 for j)
for i in range(15, 18):
    for j in range(0, 4):
        energy_matrix[i][19 - j] = e2
        energy_matrix[11][19 - j] = e2
        energy_matrix[19][19 - j] = e2
        alpha_matrix[i][19 - j] = a2
        alpha_matrix[11][19 - j] = a2
        alpha_matrix[19][19 - j] = a2
        r_matrix[i][19 - j] = (4 * r2 + r1) / 5
        r_matrix[11][19 - j] = (4 * r2 + r1) / 5
        r_matrix[19][19 - j] = (4 * r2 + r1) / 5
        A_matrix[i][19 - j] = 0
        A_matrix[11][19 - j] = 0
        A_matrix[19][19 - j] = 0

    # Specific positions for j = 8 and 9
    for j in [8, 9]:
        energy_matrix[i][19 - j] = e2
        energy_matrix[11][19 - j] = e2
        energy_matrix[19][19 - j] = e2
        alpha_matrix[i][19 - j] = a2
        alpha_matrix[11][19 - j] = a2
        alpha_matrix[19][19 - j] = a2
        r_matrix[i][19 - j] = (4 * r2 + r1) / 5
        r_matrix[11][19 - j] = (4 * r2 + r1) / 5
        r_matrix[19][19 - j] = (4 * r2 + r1) / 5
        A_matrix[i][19 - j] = 0
        A_matrix[11][19 - j] = 0
        A_matrix[19][19 - j] = 0

# p2-p1 interaction parameters
e2 = 3.71e-5
a1 = 0.662
a2 = 0.479
r2 = 17.80
r1 = 13.1
A_neg = 0
A_pos = 0

# Update matrices for p2-p1 interaction (positions 15 to 17 for i and 12 to 14 for j)
for i in range(15, 18):
    for j in range(12, 15):
        energy_matrix[i][19 - j] = e2
        energy_matrix[11][19 - j] = e2
        energy_matrix[19][19 - j] = e2
        alpha_matrix[i][19 - j] = a2
        alpha_matrix[11][19 - j] = a2
        alpha_matrix[19][19 - j] = a2
        r_matrix[i][19 - j] = (r1 + 3 * r2) / 4
        r_matrix[11][19 - j] = (r1 + 3 * r2) / 4
        r_matrix[19][19 - j] = (r1 + 3 * r2) / 4
        A_matrix[i][19 - j] = 0
        A_matrix[11][19 - j] = 0
        A_matrix[19][19 - j] = 0

# Specific position for j = 18
for i in range(15, 18):
    energy_matrix[i][19 - 18] = e2
    alpha_matrix[i][19 - 18] = a2
    r_matrix[i][19 - 18] = (r1 + 3 * r2) / 4
    A_matrix[i][19 - 18] = A_neg if i != 17 else A_pos

energy_matrix[11][19 - 18] = e2
alpha_matrix[11][19 - 18] = a2
r_matrix[11][19 - 18] = (r1 + 3 * r2) / 4
A_matrix[11][19 - 18] = 0

energy_matrix[19][19 - 18] = e2
alpha_matrix[19][19 - 18] = a2
r_matrix[19][19 - 18] = (r1 + 3 * r2) / 4
A_matrix[19][19 - 18] = A_pos


# p2-p2 interaction parameteris for ++, +-, ++
e2_1 = 5.11e-5
a2_1 = 0.393 -0.03
a1_1 = 0.419
r2_1 = 20.04
r1_1 = 7.07
A_1 =  1 * (1)* 1390 * 0.2388/(4 * math.pi * 80)

e2_2 = 3.0e-4
rcut = 30
a2_2 = 0.585 -0.05
a1_2 = 0.714
r2_2 = 14.39
r1_2 = 12.1
A_2 =  1 * (-1)* 1390 * 0.2388/(4 * math.pi * 80)

e2_3 = 2.5e-4
a2_3 = 0.360
a1_3 = 0.643
r2_3 = 18.68
r1_3 = 9.50

energy_matrix[15][19-15] = e2_1 + 0.011
alpha_matrix[15][19-15] = a2_1
r_matrix[15][19-15] = (r1_1 +  r2_1)/2
A_matrix[15][19-15] = A_1

energy_matrix[16][19-15] = e2_1 + 0.011
alpha_matrix[16][19-15] = a2_1
r_matrix[16][19-15] = (r1_1 + r2_1)/2
A_matrix[16][19-15] = A_1

energy_matrix[17][19-15] = e2_2
alpha_matrix[17][19-15] = a2_2
r_matrix[17][19-15] = (r1_2 + r2_2)/2
A_matrix[17][19-15] = A_2

energy_matrix[15][19-16] = e2_1 + 0.011
alpha_matrix[15][19-16] = a2_1
r_matrix[15][19-16] = (r1_1 + r2_1)/2
A_matrix[15][19-16] = A_1

energy_matrix[16][19-16] = e2_1 + 0.011
alpha_matrix[16][19-16] = a2_1
r_matrix[16][19-16] = (r1_1 + r2_1)/2
A_matrix[16][19-16] = A_1

energy_matrix[17][19-16] = e2_2
alpha_matrix[17][19-16] = a2_2
r_matrix[17][19-16] =(r1_2 + r2_2)/2
A_matrix[17][19-16] = A_2

energy_matrix[15][19-17] = e2_2
alpha_matrix[15][19-17] = a2_2
r_matrix[15][19-17] = (r1_2 + r2_2)/2
A_matrix[15][19-17] = A_2

energy_matrix[16][19-17] = e2_2
alpha_matrix[16][19-17] = a2_2
r_matrix[16][19-17] = (r1_2 + r2_2)/2
A_matrix[16][19-17] = A_2

energy_matrix[17][19-17] = e2_1 + 0.011
alpha_matrix[17][19-17] = a2_1
r_matrix[17][19-17] = (r1_1 + r2_1)/2
A_matrix[17][19-17] = A_1

energy_matrix[11][19-15] = e2_3
alpha_matrix[11][19-15] = a2_3
r_matrix[11][19-15] = (9* r2_3 + r1_3)/10
A_matrix[11][19-15] = 0

energy_matrix[11][19-16] = e2_3
alpha_matrix[11][19-16] = a2_3
r_matrix[11][19-16] = (9* r2_3 + r1_3)/10
A_matrix[11][19-16] = 0

energy_matrix[11][19-17] = e2_3
alpha_matrix[11][19-17] = a2_3
r_matrix[11][19-17] = (9* r2_3 + r1_3)/10
A_matrix[11][19-17] = 0

energy_matrix[17][19-11] = e2_3
alpha_matrix[17][19-11] = a2_3
r_matrix[17][19-11] = (9* r2_3 + r1_3)/10
A_matrix[17][19-11] = 0

energy_matrix[19][19-15] = e2_2
alpha_matrix[19][19-15] = a2_2
r_matrix[19][19-15] = (r1_2 + r2_2)/2
A_matrix[19][19-15] = A_2

energy_matrix[19][19-16] = e2_2
alpha_matrix[19][19-16] = a2_2
r_matrix[19][19-16] = r1_2
A_matrix[19][19-16] = A_2

energy_matrix[19][19-17] = e2_1 + 0.011
alpha_matrix[19][19-17] = a2_1
r_matrix[19][19-17] = (r1_1 + r2_1)/2
A_matrix[19][19-17] = A_1

energy_matrix[15][19-11] = e2_3
alpha_matrix[15][19-11] = a2_3
r_matrix[15][19-11] = (9* r2_3 + r1_3)/10
A_matrix[15][19-11] = 0

energy_matrix[16][19-11] = e2_3
alpha_matrix[16][19-11] = a2_3
r_matrix[16][19-11] = (9* r2_3 + r1_3)/10
A_matrix[16][19-11] = 0

energy_matrix[15][19-19] =e2_2
alpha_matrix[15][19-19] =a2_2
r_matrix[15][19-19] = (r1_2 + r2_2)/2
A_matrix[15][19-19] = A_2

energy_matrix[16][19-19] =e2_2
alpha_matrix[16][19-19] =a2_2
r_matrix[16][19-19] =(r1_2 + r2_2)/2
A_matrix[16][19-19] =A_2

energy_matrix[17][19-19] =e2_1 + 0.011
alpha_matrix[17][19-19] =a2_1
r_matrix[17][19-19] =(r1_1 + r2_1)/2
A_matrix[17][19-19] =A_1

energy_matrix[11][19-11] = e2_3
alpha_matrix[11][19-11] = a2_3
r_matrix[11][19-11] = (9* r2_3 + r1_3)/10
A_matrix[11][19-11] = 0

energy_matrix[11][19-19] =e2_3
alpha_matrix[11][19-19] =a2_3
r_matrix[11][19-19] =(9* r2_3 + r1_3)/10
A_matrix[11][19-19] = 0

energy_matrix[19][19-11] = e2_3
alpha_matrix[19][19-11] = a2_3
r_matrix[19][19-11] = (9* r2_3 + r1_3)/10
A_matrix[19][19-11] = 0

energy_matrix[19][19-19] =e2_1 + 0.011
alpha_matrix[19][19-19] =a2_1
r_matrix[19][19-19] =(r1_1 + r2_1)/2
A_matrix[19][19-19] =A_1

# Convert lists to numpy arrays for matrix operations
energy_matrix = np.array(energy_matrix)
alpha_matrix = np.array(alpha_matrix)
r_matrix = np.array(r_matrix)
A_matrix = np.array(A_matrix)

# Write the nonbonding parameters to a text file for use in simulations
with open('nonbond_Trovato.txt', 'w') as txtfile:
    for i in range(len(desired_order)):
        for j in range(len(desired_order)-i):
            txtfile.write(f"pair_coeff\t{i+1}\t{20-j}\tmorse\t{energy_matrix[i,j]:.2f}\t{alpha_matrix[i,j]:.2f}\t{r_matrix[i,j]:.2f}\t24\n")
            txtfile.write(f"pair_coeff\t{i+1}\t{20-j}\tyukawa\t{A_matrix[i,j]:.2f}\t40\n")

# Function to create and display masked heatmaps for the matrices
def plot_masked_heatmap(matrix, cmap, vmin, vmax, colorbar_label, axis_labels):
    n = matrix.shape[0]
    mask = np.zeros_like(matrix, dtype=bool)
    # Mask elements above the diagonal for visualization
    for i in range(n):
        for j in range(n):
            if i + j <= n - 1:
                mask[i, j] = True
    masked_matrix = np.ma.array(matrix, mask=~mask)
    
    # Plot the masked matrix
    plt.imshow(masked_matrix, cmap=cmap, vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(label=colorbar_label)
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label(colorbar_label, fontsize=18)
    
    # Map amino acids to single-letter codes for axis labels
    aa_mapping = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    single_letter_order = [aa_mapping[aa] for aa in desired_order]
    
    plt.xticks(np.arange(len(single_letter_order)), single_letter_order[::-1], fontsize=16)
    plt.yticks(np.arange(len(single_letter_order)), single_letter_order, fontsize=16)

    plt.show()

plot_masked_heatmap(energy_matrix, cm.jet, 0, 0.2, 'kcal/mol', 'Energy')
plot_masked_heatmap(alpha_matrix, cm.jet, 0.3, 1.3, '1/$\AA$', 'Alpha')
plot_masked_heatmap(r_matrix, cm.jet, 5.5, 16, '$\AA$', 'Radius')
plot_masked_heatmap(A_matrix, cm.jet, -0.3, 0.3, 'kcal/mol', 'A')
