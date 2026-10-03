"""
Prepare a protein structure for docking.

Removes waters and heteroatoms, keeps only standard amino acids.

Usage:
    python prepare_protein.py <input_pdb> <output_pdb>
"""
import sys
from Bio.PDB import PDBParser, PDBIO, Select


class ProteinSelect(Select):
    """Select only protein atoms (standard amino acids)."""

    def accept_residue(self, residue):
        return residue.id[0] == ' '

    def accept_atom(self, atom):
        return atom.element != 'H'


def prepare_protein(input_pdb, output_pdb):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, ProteinSelect())

    print(f"Prepared protein saved to {output_pdb}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python prepare_protein.py <input_pdb> <output_pdb>")
        sys.exit(1)

    prepare_protein(sys.argv[1], sys.argv[2])
