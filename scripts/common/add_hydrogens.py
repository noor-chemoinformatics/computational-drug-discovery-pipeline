"""
Add explicit hydrogens to a ligand SDF file.

Usage:
    python add_hydrogens.py <input.sdf> <output.sdf>
"""
import sys
from rdkit import Chem
from rdkit.Chem import AllChem


def add_hydrogens(input_sdf, output_sdf):
    """Read SDF, add hydrogens, and write new SDF."""
    supplier = Chem.SDMolSupplier(input_sdf, removeHs=False)
    
    writer = Chem.SDWriter(output_sdf)
    count = 0
    
    for mol in supplier:
        if mol is None:
            continue
        
        # Add explicit hydrogens
        mol = Chem.AddHs(mol, addCoords=True)
        
        # Write
        writer.write(mol)
        count += 1
    
    writer.close()
    print(f"Processed {count} molecule(s)")
    print(f"Saved to {output_sdf}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python add_hydrogens.py <input.sdf> <output.sdf>")
        sys.exit(1)
    add_hydrogens(sys.argv[1], sys.argv[2])
