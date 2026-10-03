"""
Find the center of a ligand in a PDB file to define the docking box.

Usage:
    python find_binding_site.py <pdb_file> <ligand_resname>
"""
import sys
import numpy as np
from Bio.PDB import PDBParser


def find_ligand_center(pdb_file, ligand_resname):
    """Calculate the center of mass of a ligand."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    coords = []
    for residue in structure.get_residues():
        if residue.get_resname().strip() == ligand_resname:
            for atom in residue.get_atoms():
                coords.append(atom.get_coord())

    if not coords:
        print(f"No atoms found for ligand {ligand_resname}")
        return None

    coords = np.array(coords)
    center = coords.mean(axis=0)

    print(f"Ligand: {ligand_resname}")
    print(f"Number of atoms: {len(coords)}")
    print(f"Center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
    print(f"Box size suggestion: 25-30 Å")

    # Also print the range for reference
    print(f"\nCoordinate ranges:")
    print(f"  X: {coords[:, 0].min():.2f} to {coords[:, 0].max():.2f}")
    print(f"  Y: {coords[:, 1].min():.2f} to {coords[:, 1].max():.2f}")
    print(f"  Z: {coords[:, 2].min():.2f} to {coords[:, 2].max():.2f}")

    return center


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python find_binding_site.py <pdb_file> <ligand_resname>")
        sys.exit(1)
    find_ligand_center(sys.argv[1], sys.argv[2])
