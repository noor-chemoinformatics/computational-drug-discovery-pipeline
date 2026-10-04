"""
Simple interaction fingerprint: count contacts within distance cutoffs.

Usage:
    python calculate_ifp_simple.py <protein.pdb> <ligand.pdb> <output.csv>
"""
import sys
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser


def load_atoms(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("mol", pdb_file)
    atoms = []
    for atom in structure.get_atoms():
        atoms.append({
            'name': atom.get_name(),
            'resname': atom.get_parent().get_resname(),
            'resid': atom.get_parent().get_id()[1],
            'element': atom.element,
            'coord': atom.get_coord(),
        })
    return atoms


def calculate_contacts(protein_pdb, ligand_pdb, output_csv, cutoff=4.0):
    protein_atoms = load_atoms(protein_pdb)
    ligand_atoms = load_atoms(ligand_pdb)

    print(f"Protein atoms: {len(protein_atoms)}")
    print(f"Ligand atoms: {len(ligand_atoms)}")

    contacts = []
    for p_atom in protein_atoms:
        if p_atom['element'] == 'H':
            continue
        for l_atom in ligand_atoms:
            if l_atom['element'] == 'H':
                continue
            dist = np.linalg.norm(p_atom['coord'] - l_atom['coord'])
            if dist <= cutoff:
                contacts.append({
                    'protein_residue': f"{p_atom['resname']}{p_atom['resid']}",
                    'protein_atom': p_atom['name'],
                    'ligand_atom': l_atom['name'],
                    'distance': round(dist, 2),
                })

    df = pd.DataFrame(contacts)
    df.to_csv(output_csv, index=False)
    print(f"\nFound {len(df)} contacts within {cutoff} A")

    if len(df) > 0:
        print("\nTop interacting residues:")
        print(df['protein_residue'].value_counts().head(10))
    else:
        print("No contacts found — try increasing the cutoff")

    return df


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python calculate_ifp_simple.py <protein.pdb> <ligand.pdb> <output.csv>")
        sys.exit(1)
    calculate_contacts(sys.argv[1], sys.argv[2], sys.argv[3])
