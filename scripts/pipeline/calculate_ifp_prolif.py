"""
Calculate interaction fingerprint with ProLIF.

Usage:
    python calculate_ifp_prolif.py <protein.pdb> <ligand.pdb> <output.csv>
"""
import sys
import warnings
warnings.filterwarnings('ignore')

import MDAnalysis as mda
import prolif as plf
import pandas as pd
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def calculate_ifp(protein_pdb, ligand_pdb, output_csv):
    print(f"Loading protein: {protein_pdb}")
    print(f"Loading ligand: {ligand_pdb}")
    
    # Load with MDAnalysis
    protein_u = mda.Universe(protein_pdb)
    protein = plf.Molecule.from_mda(protein_u)
    
    ligand_u = mda.Universe(ligand_pdb)
    ligand = plf.Molecule.from_mda(ligand_u)
    
    # Generate fingerprint
    fp = plf.Fingerprint()
    fp.run_from_iterable([ligand], protein)
    
    # Convert to DataFrame
    df = fp.to_dataframe()
    df.to_csv(output_csv)
    
    print(f"Saved to {output_csv}")
    print(f"Total interaction types detected: {df.shape[1]}")
    
    if len(df) > 0:
        row = df.iloc[0]
        active = row[row > 0]
        print(f"\nDetected interactions ({len(active)}):")
        for col, val in active.items():
            print(f"  {col}")
    
    return df


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python calculate_ifp_prolif.py <protein.pdb> <ligand.pdb> <output.csv>")
        sys.exit(1)
    calculate_ifp(sys.argv[1], sys.argv[2], sys.argv[3])
