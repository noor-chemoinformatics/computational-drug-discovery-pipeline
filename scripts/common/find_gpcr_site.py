"""
Estimate GPCR binding site from protein structure.

For class A GPCRs, the orthosteric site is roughly at the center
of the transmembrane helical bundle.

Usage:
    python find_gpcr_site.py <protein.pdb> <output.txt>
"""
import sys
import numpy as np
from Bio.PDB import PDBParser


def estimate_center(pdb_file, output_txt):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    ca_coords = []
    for residue in structure.get_residues():
        if 'CA' in residue:
            ca_coords.append(residue['CA'].get_coord())

    ca_coords = np.array(ca_coords)
    center = ca_coords.mean(axis=0)

    output = f"Estimated binding site center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})"
    print(output)
    print("Box size suggestion: 40 x 40 x 40 Angstrom (blind docking)")

    with open(output_txt, 'w') as f:
        f.write(output + "\n")
        f.write(f"Center: {center[0]:.2f} {center[1]:.2f} {center[2]:.2f}\n")
        f.write("Box size: 40 40 40\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python find_gpcr_site.py <protein.pdb> <output.txt>")
        sys.exit(1)
    estimate_center(sys.argv[1], sys.argv[2])
