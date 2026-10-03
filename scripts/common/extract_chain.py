"""
Extract a specific chain from a PDB file.

Usage:
    python extract_chain.py <input.pdb> <chain_id> <output.pdb>
"""
import sys
from Bio.PDB import PDBParser, PDBIO, Select


class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id


def extract_chain(input_pdb, chain_id, output_pdb):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, ChainSelect(chain_id))
    print(f"Extracted chain {chain_id} to {output_pdb}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_chain.py <input.pdb> <chain_id> <output.pdb>")
        sys.exit(1)
    extract_chain(sys.argv[1], sys.argv[2], sys.argv[3])
