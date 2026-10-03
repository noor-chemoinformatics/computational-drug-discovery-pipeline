 # Computational-drug-discovery-pipeline
Computational drug discovery pipeline for structure-based virtual screening. Combines deep docking (ML-accelerated screening of ultra-large libraries), multi-criteria hit selection (consensus docking, ADMET, diversity, interaction fingerprints), relative binding free energy calculations (OpenFE), and molecular dynamics.
<br>
A modular, reproducible pipeline for **structure-based virtual screening** from ultra-large chemical libraries down to a prioritized shortlist of compounds for experimental testing.

Combines **deep docking** (machine learning-accelerated screening), **multi-criteria hit selection** (consensus docking, ADMET, chemical diversity, interaction fingerprints), **relative binding free energy (RBFE) calculations** with OpenFE, and **molecular dynamics** for stability and mechanism analysis.

Designed as a **generalizable toolkit** applicable to any protein target with a known or predicted structure.
# Computational Drug Discovery Pipeline

Built with Python, RDKit, AutoDock Vina, and Meeko. Validated on a test system (lysozyme + NAG) and demonstrated on a 500-compound chemical library.

**Status:** Working prototype — validated on a test system, not yet applied to a real drug target.

 **What This Pipeline Does**

This project implements the core computational workflow used in early-stage drug discovery:

- Single docking — dock a known ligand to a known protein and validate against experimental data
- Batch docking — screen a chemical library against a target and rank by binding affinity
- Deep docking — use machine learning to screen ultra-large libraries efficiently
- Hit selection — filter docking hits by drug-likeness and physicochemical properties
- Visualization — plot score distributions and identify top candidates

The pipeline is generalizable — any protein with a known or predicted structure can be used as a target.

**Pipeline at a Glance**

1. PROTEIN PREPARATION — PDB -> clean -> add hydrogens -> PDBQT
2. LIGAND PREPARATION — SMILES -> 3D conformer -> hydrogens -> PDBQT
3. MOLECULAR DOCKING — AutoDock Vina -> ranked poses -> binding affinities
4. DEEP DOCKING (ML-ACCELERATED) — sample -> train Random Forest -> predict -> filter -> repeat
5. HIT SELECTION — Lipinski filter -> property calculation -> top N shortlist

**Results**

### 1. Single Docking — Validation on a Test System

| Metric | Value |
|--------|-------|
| Protein | Lysozyme (PDB: 1HEW) |
| Ligand | NAG (tri-N-acetylglucosamine) |
| Best docking affinity | -6.5 kcal/mol |
| Experimental range | -5 to -7 kcal/mol |
| Result | Matches experimental range |

### 2. Batch Docking — 98 ZINC Compounds

| Metric | Value |
|--------|-------|
| Compounds docked | 98 |
| Best affinity | -7.55 kcal/mol |
| Compounds better than NAG | 10+ |
| Mean affinity | ~-5.8 kcal/mol |

### 3. Deep Docking — ML-Accelerated Screening

| Iteration | Library size | Model R2 | RMSE (kcal/mol) |
|-----------|--------------|----------|-----------------|
| 1 | 500 -> 100 | 0.353 | 0.68 |
| 2 | 100 -> 100 | 0.354 | 0.67 |
| 3 | 100 -> 100 | 0.246 | 0.73 |

The Random Forest model predicted docking scores with RMSE < 0.8 kcal/mol, accurate enough for ranking. Deep docking compressed a 500-compound library to 100 top candidates while only docking 300 compounds.

### 4. Hit Selection — Drug-Like Filtering

Compounds were filtered using Lipinski's Rule of Five (MW <= 500, logP <= 5, HBD <= 5, HBA <= 10) and physicochemical properties calculated with RDKit.

Top 20 drug-like hits: results/tables/top_hits_shortlist.csv
### 5. β2AR Validation (Real GPCR Target)

To validate the pipeline on a real drug target, we applied it to the beta-2 adrenergic receptor (β2AR), a well-characterized GPCR.

| Metric | Value |
|--------|-------|
| Protein | β2AR (PDB: 2RH1) |
| Ligand | Carazolol (inverse agonist, Ki ~1 nM) |
| Best docking affinity | -10.04 kcal/mol |
| RMSD vs crystal pose | 1.74 Å |
| Result | Validation PASSED |

This confirms the pipeline can recover correct binding poses for a real GPCR target.
## Installation

Requires Linux (Ubuntu 20.04+ or WSL2 on Windows) and Miniconda.

    git clone https://github.com/noor-chemoinformatics/computational-drug-discovery-pipeline.git
    cd computational-drug-discovery-pipeline
    conda env create -f environment.yml
    conda activate cddp
    python -c "import numpy, pandas, rdkit, sklearn, prolif, meeko; print('All packages OK')"
    vina --version

## 5. Project Structure

    computational-drug-discovery-pipeline/
    |-- README.md
    |-- LICENSE
    |-- environment.yml
    |-- requirements.txt
    |-- data/
    |   |-- raw/                       # Raw downloaded data
    |   |-- processed/                 # Cleaned, prepared PDBQT files
    |-- scripts/
    |   |-- common/                    # Reusable utility scripts
    |   |-- pipeline/                  # Main pipeline scripts
    |   |-- targets/                   # Target-specific scripts
    |-- results/
    |   |-- tables/                    # CSV results
    |   |-- figures/                   # Plots
    |   |-- deep_docking_test/         # Deep docking outputs
    |-- notebooks/
    |-- docs/
    |-- tests/

## 6. Methods

| Stage | Tool | Purpose |
|-------|------|---------|
| Protein preparation | Biopython, Meeko | Clean PDB, add H, convert to PDBQT |
| Ligand preparation | RDKit, Meeko | 3D conformer, hydrogens, PDBQT |
| Docking | AutoDock Vina 1.2.7 | Binding pose and affinity prediction |
| Fingerprinting | RDKit (Morgan, r=2, 2048 bits) | Molecular representation for ML |
| ML model | Scikit-learn (Random Forest) | Predict docking scores |
| Drug-likeness | RDKit (Lipinski descriptors) | Filter for drug-like compounds |
| Visualization | Matplotlib | Score distributions, metrics |

## 7. Limitations

1. Validated on a test system only — lysozyme, not a real drug target
2. Small library — 500 compounds for deep docking; real runs use millions
3. No experimental validation — all results are computational predictions
4. No RBFE or MD yet — these are planned next
5. ML applicability domain — Random Forest is limited to training chemical space

## 8. Roadmap

- [x] Single docking validated against known complex
- [x] Batch docking of chemical library
- [x] Deep docking (ML-accelerated screening)
- [x] Drug-likeness filtering and hit selection
- [ ] Apply pipeline to a real GPCR target (B2AR, GPR35)
- [ ] Relative binding free energy calculations (OpenFE)
- [ ] Molecular dynamics simulations (GROMACS)
- [ ] Interaction fingerprints (ProLIF)
- [ ] Scale to ultra-large libraries (millions of compounds)

## 9. Author

**Manahil Noor**
Master 2 Student
manahilnoor339@gmail.com
GitHub: https://github.com/noor-chemoinformatics

## License

This project is licensed under the MIT License — see LICENSE for details.
