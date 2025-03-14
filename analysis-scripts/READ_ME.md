# Trajectory Analysis for Coarse-Grained Protein Simulations

This set of Python scripts provides post-processing tools for analyzing generated trajectories of coarse-grained protein systems. The tools allow for:

- Mean Squared Displacement (MSD) analysis to study protein diffusion.
- Protein contact dynamics to track interactions over time.
- Contact probability matrix computation to assess interaction patterns.

## 1. compute_msd.py

### Description:
Computes the Mean Squared Displacement (MSD) for a system of proteins, averaging over all proteins in the simulation.

### Workflow:
1. Reads a PDB trajectory file containing the center-of-mass (CM) positions of proteins.
2. Computes the MSD while accounting for periodic boundary conditions (PBCs).
3. Uses a rolling average for improved statistical accuracy.
4. Outputs a text file with MSD values and associated errors, along with a plot.

---

## 2. contacts_evolution.py

### Description:
Analyzes the evolution of protein-protein contacts over time based on a distance cutoff.

### Workflow:
1. Reads a PDB trajectory file containing multiple frames.
2. Computes pairwise distances between proteins.
3. Identifies unique protein-protein contacts in each frame.
4. Outputs a `contacts_over_time.txt` file storing the time and number of contacts.

---

## 3. prob_contact_matrix.py

### Description:
Computes the contact probability matrix for proteins in a trajectory, providing insight into interaction frequency.

### Workflow:
1. Reads a PDB trajectory file containing multiple frames.
2. Calculates pairwise distances between residues or beads across different proteins.
3. Identifies contacts based on a distance cutoff and updates the contact probability matrix.
4. Normalizes the matrix to reflect contact probabilities over the entire trajectory.
5. Outputs the contact probability matrix file and generates a visual plot.
