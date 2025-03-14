# Non-Bonded Interaction Parameters: Generation & Visualization

This suite of Python scripts generates non-bonded interaction parameters for the **CALVADOS, COCOMO, and Trovato** force fields in **LAMMPS** format. The scripts:

- Read input files containing interaction parameters.
- Format and export the data into LAMMPS-compatible output files.
- Generate visual maps of the interaction matrices for analysis.

---

## 1. write_nonbonded_calvados.py

### Description:
Generates non-bonded interaction parameters for the **CALVADOS** force field.

### Workflow:
1. Reads two input files:
   - `calvados_radius.txt` → Contains radius data.
   - `calvados.csv` → Contains interaction parameters.
2. Writes the formatted interaction parameters to `nonbond_calvados.txt`.
3. Generates visual interaction maps for CALVADOS.

---

## 2. write_nonbonded_cocomo.py

### Description:
Generates non-bonded interaction parameters for the **COCOMO** force field.

### Workflow:
1. Reads two input files:
   - `cocomo_att.txt` → Contains attractive interaction parameters.
   - `cocomo_rep.txt` → Contains repulsive interaction parameters.
2. Writes the formatted interaction parameters to `nonbond_COCOMO.txt`.
3. Generates visual interaction maps for COCOMO.

---

## 3. write_nonbonded_Trovato.py

### Description:
Generates non-bonded interaction parameters for the **Trovato** force field.

### Workflow:
1. Processes built-in data (**no external input files required**).
2. Writes the formatted interaction parameters to `nonbond_Trovato.txt`.
3. Generates visual interaction maps for Trovato.
