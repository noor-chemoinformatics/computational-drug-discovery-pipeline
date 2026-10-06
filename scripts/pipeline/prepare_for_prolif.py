"""
Prepare protein and ligand for ProLIF analysis.

Usage:
    python prepare_for_prolif.py <protein.pdb> <ligand.pdb> <output_dir>
"""
import sys
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def prepare_ligand(ligand_pdb, output_pdb):
    """Add hydrogens and assign bond orders to ligand."""
    mol = Chem.MolFromPDBFile(ligand_pdb, removeHs=False)
    if mol is None:
        print(f"Failed to load ligand: {ligand_pdb}")
        return None
    
    mol = Chem.AddHs(mol, addCoords=True)
    Chem.MolToPDBFile(mol, output_pdb)
    print(f"Ligand prepared: {output_pdb}")
    return mol


def prepare_protein(protein_pdb, output_pdb):
    """Prepare protein by adding hydrogens."""
    mol = Chem.MolFromPDBFile(protein_pdb, removeHs=False)
    if mol is None:
        print(f"Failed to load protein: {protein_pdb}")
        return None
    
    mol = Chem.AddHs(mol, addCoords=True)
    Chem.MolToPDBFile(mol, output_pdb)
    print(f"Protein prepared: {output_pdb}")
    return mol


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python prepare_for_prolif.py <protein.pdb> <ligand.pdb> <output_dir>")
        sys.exit(1)
    
    protein_pdb = sys.argv[1]
    ligand_pdb = sys.argv[2]
    output_dir = sys.argv[3]
    
    os.makedirs(output_dir, exist_ok=True)
    
    prepare_protein(protein_pdb, os.path.join(output_dir, "protein_h.pdb"))
    prepare_ligand(ligand_pdb, os.path.join(output_dir, "ligand_h.pdb"))
