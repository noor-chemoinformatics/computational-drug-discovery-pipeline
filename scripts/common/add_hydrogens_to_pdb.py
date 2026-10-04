"""
Add explicit hydrogens to a PDB file using RDKit.

Usage:
    python add_hydrogens_to_pdb.py <input.pdb> <output.pdb>
"""
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def add_hydrogens(input_pdb, output_pdb):
    """Read PDB, add hydrogens, write PDB."""
    mol = Chem.MolFromPDBFile(input_pdb, removeHs=False)
    if mol is None:
        print(f"Failed to load {input_pdb}")
        return False

    # Add hydrogens (with coordinates)
    mol = Chem.AddHs(mol, addCoords=True)

    # Write
    Chem.MolToPDBFile(mol, output_pdb)
    print(f"Added hydrogens and saved to {output_pdb}")
    return True


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python add_hydrogens_to_pdb.py <input.pdb> <output.pdb>")
        sys.exit(1)
    add_hydrogens(sys.argv[1], sys.argv[2])
