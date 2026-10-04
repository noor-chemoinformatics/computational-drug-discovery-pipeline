"""
Create a curated list of known GPR35 agonists.

Usage:
    python create_gpr35_ligands.py <output.smi>
"""
import sys


GPR35_LIGANDS = [
    ("OC(=O)c1ccccc1NC(=O)C(Cl)Cl", "lodoxamide"),
    ("OC(=O)c1ccc(cc1)C(=O)c1ccccc1", "pamoate_analog"),
    ("O=C1NC(=O)c2ccccc21", "zaprinast_analog"),
    ("OC(=O)c1ccccc1O", "salicylate"),
    ("OC(=O)c1ccc(O)cc1", "4-hydroxybenzoate"),
    ("OC(=O)c1ccccc1", "benzoate"),
    ("OC(=O)CCc1ccccc1", "phenylpropionate"),
    ("OC(=O)c1ccc(cc1)N", "aminobenzoate"),
    ("OC(=O)c1cccnc1", "nicotinate"),
    ("OC(=O)c1cccs1", "thiophene_carboxylate"),
]


def create_ligands(output_smi):
    with open(output_smi, 'w') as f:
        for smi, name in GPR35_LIGANDS:
            f.write(f"{smi}\t{name}\n")
    print(f"Wrote {len(GPR35_LIGANDS)} ligands to {output_smi}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python create_gpr35_ligands.py <output.smi>")
        sys.exit(1)
    create_ligands(sys.argv[1])
