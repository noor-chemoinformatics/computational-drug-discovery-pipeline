"""
Extract a specific ligand from a PDB file.

Usage:
    python extract_ligand.py <pdb_file> <ligand_resname> <output_file>
"""
import sys
from Bio.PDB import PDBParser, PDBIO, Select


class LigandSelect(Select):
    """Select only atoms belonging to a specific ligand."""

    def __init__(self, ligand_resname):
        self.ligand_resname = ligand_resname

    def accept_residue(self, residue):
        return residue.get_resname().strip() == self.ligand_resname


def extract_ligand(pdb_file, ligand_resname, output_file):
    """Extract ligand from PDB and save to a new file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file, LigandSelect(ligand_resname))

    print(f"Extracted {ligand_resname} to {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_ligand.py <pdb_file> <ligand_resname> <output_file>")
        sys.exit(1)

    extract_ligand(sys.argv[1], sys.argv[2], sys.argv[3])
