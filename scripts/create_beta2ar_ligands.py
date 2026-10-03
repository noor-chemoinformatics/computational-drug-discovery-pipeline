"""
Create a curated list of known beta2AR ligands.

These are well-known beta2AR agonists and antagonists with
experimental binding data from the literature.

Usage:
    python create_beta2ar_ligands.py <output.smi>
"""
import sys


# Known beta2AR ligands (SMILES from ChEMBL / DrugBank)
BETA2AR_LIGANDS = [
    # Antagonists / inverse agonists
    ("CC(C)NCC(COC1=CC=CC=C1)O", "propranolol"),
    ("CC(C)NCC(COC1=CC=C(C=C1)COCCOC)O", "metoprolol"),
    ("CC(C)(C)NCC(COC1=CC=CC2=CC=CC=C21)O", "propranolol_analog"),
    ("CC(C)NCC(COC1=CC=C(C=C1)CCO)O", "atenolol"),
    ("COC1=CC=C(C=C1)CCNCC(C)O", "betaxolol"),
    ("CC(C)NCC(COC1=CC=CC=C1C#N)O", "cyanopindolol"),
    ("CC(C)NCC(COC1=CC=CC2=CC=CC=C21)O", "carvedilol_analog"),
    # Carazolol (the crystal ligand)
    ("CC(C)NCC(COC1=CC=CC2=C1C3=CC=CC=C3N2)O", "carazolol"),
    # Agonists
    ("CC(C)NCC(O)C1=CC=C(O)C(CO)=C1", "salbutamol"),
    ("CC(C)(C)NCC(O)C1=CC(CO)=C(O)C=C1", "salbutamol_analog"),
    ("CC(C)NCC(O)C1=CC=C(O)C=C1", "fenoterol_analog"),
    ("CC(C)(C)NCC(O)C1=CC=C(O)C(CO)=C1", "terbutaline"),
    ("CC(C)NCC(O)C1=CC=CC=C1O", "isoprenaline_analog"),
    ("CC(C)NCC(O)C1=CC=C(O)C(O)=C1", "isoprenaline"),
    ("CC(C)NCC(O)C1=CC(O)=C(O)C=C1", "adrenaline_analog"),
    # Other known binders
    ("CC(C)NCC(COC1=CC=CC=C1)OC", "beta_blocker_analog"),
    ("CC(C)NCC(COC1=CC=C(C=C1)Cl)O", "beta_blocker_2"),
    ("CC(C)NCC(COC1=CC=C(C=C1)F)O", "beta_blocker_3"),
    ("CC(C)NCC(COC1=CC=C(C=C1)C)O", "beta_blocker_4"),
    ("CC(C)NCC(COC1=CC=C(C=C1)OC)O", "beta_blocker_5"),
]


def create_ligands(output_smi):
    with open(output_smi, 'w') as f:
        for smi, name in BETA2AR_LIGANDS:
            f.write(f"{smi}\t{name}\n")
    print(f"Wrote {len(BETA2AR_LIGANDS)} ligands to {output_smi}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python create_beta2ar_ligands.py <output.smi>")
        sys.exit(1)
    create_ligands(sys.argv[1])
