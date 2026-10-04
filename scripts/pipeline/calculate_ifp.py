"""
Calculate protein-ligand interaction fingerprints with ProLIF.

Usage:
    python calculate_ifp.py <protein.pdb> <ligand.pdbqt> <output.csv>
"""
import sys
import warnings
warnings.filterwarnings('ignore')

import MDAnalysis as mda
import prolif as plf
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def pdbqt_to_mol(pdbqt_file):
    """Convert PDBQT to RDKit mol via PDB."""
    # ProLIF can read PDBQT directly via MDAnalysis
    u = mda.Universe(pdbqt_file)
    return u


def calculate_ifp(protein_pdb, ligand_pdbqt, output_csv):
    """Calculate interaction fingerprint."""
    print(f"Loading protein: {protein_pdb}")
    print(f"Loading ligand: {ligand_pdbqt}")

    # Load protein
    protein_u = mda.Universe(protein_pdb)
    protein = plf.Molecule.from_mda(protein_u)

    # Load ligand
    ligand_u = mda.Universe(ligand_pdbqt)
    ligand = plf.Molecule.from_mda(ligand_u)

    # Calculate IFP
    fp = plf.Fingerprint()
    fp.run_from_iterable([ligand], protein)

    # Convert to DataFrame
    df = fp.to_dataframe()

    # Save
    df.to_csv(output_csv)
    print(f"Saved interaction fingerprint to {output_csv}")

    # Print summary
    print(f"\nInteractions found: {df.shape[1]}")
    print(f"Interaction types: {df.columns.tolist()[:10]}...")

    return df


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python calculate_ifp.py <protein.pdb> <ligand.pdbqt> <output.csv>")
        sys.exit(1)
    calculate_ifp(sys.argv[1], sys.argv[2], sys.argv[3])
