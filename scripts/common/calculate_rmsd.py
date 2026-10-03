"""
Calculate RMSD between docked pose and crystal ligand.

Usage:
    python calculate_rmsd.py <crystal.pdb> <docked.pdbqt>
"""
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def load_coords_from_pdb(pdb_file):
    """Load heavy-atom coordinates from PDB."""
    coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                element = line[76:78].strip() if len(line) > 76 else line[12:16].strip()[0]
                if element != 'H':
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
    return np.array(coords)


def calculate_rmsd(crystal_pdb, docked_pdbqt):
    crystal = load_coords_from_pdb(crystal_pdb)
    docked = load_coords_from_pdb(docked_pdbqt)

    if len(crystal) != len(docked):
        print(f"Atom count mismatch: crystal {len(crystal)}, docked {len(docked)}")
        print("Using minimum atom count")
        n = min(len(crystal), len(docked))
        crystal = crystal[:n]
        docked = docked[:n]

    # Simple RMSD (not aligned — same coordinate frame)
    diff = crystal - docked
    rmsd = np.sqrt((diff ** 2).sum() / len(crystal))

    print(f"RMSD: {rmsd:.2f} Å")
    if rmsd < 2.0:
        print("Pose matches crystal structure (RMSD < 2 Å) — validation PASSED")
    else:
        print("Pose differs from crystal structure (RMSD > 2 Å) — check binding site")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python calculate_rmsd.py <crystal.pdb> <docked.pdbqt>")
        sys.exit(1)
    calculate_rmsd(sys.argv[1], sys.argv[2])
