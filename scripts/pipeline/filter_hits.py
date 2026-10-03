"""
Filter docking hits by drug-likeness and physicochemical properties.

Usage:
    python filter_hits.py <docking_results.csv> <ligand_id_map.csv> <output_csv>
"""
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski


def calc_properties(smiles):
    """Calculate physicochemical properties."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    return {
        "MW": Descriptors.MolWt(mol),
        "logP": Descriptors.MolLogP(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "TPSA": Descriptors.TPSA(mol),
        "RotBonds": Lipinski.NumRotatableBonds(mol),
        "AromaticRings": Descriptors.NumAromaticRings(mol),
    }


def passes_lipinski(props):
    """Lipinski's Rule of Five."""
    if props is None:
        return False

    violations = 0
    if props["MW"] > 500: violations += 1
    if props["logP"] > 5: violations += 1
    if props["HBD"] > 5: violations += 1
    if props["HBA"] > 10: violations += 1

    return violations <= 1


def main(docking_csv, id_map_csv, output_csv):
    # Load docking results
    dock_df = pd.read_csv(docking_csv)
    dock_df["ligand_id"] = dock_df["ligand"].str.replace(".pdbqt", "", regex=False)

    # Load ID map
    id_map = pd.read_csv(id_map_csv)

    # Merge
    df = dock_df.merge(id_map, on="ligand_id", how="left")

    # Calculate properties
    props_list = []
    for smi in df["smiles"]:
        props_list.append(calc_properties(smi))

    props_df = pd.DataFrame(props_list)
    df = pd.concat([df, props_df], axis=1)

    # Apply Lipinski filter
    df["passes_lipinski"] = df.apply(
        lambda row: passes_lipinski({
            "MW": row["MW"],
            "logP": row["logP"],
            "HBD": row["HBD"],
            "HBA": row["HBA"],
        }), axis=1
    )

    # Sort by affinity
    df = df.sort_values("affinity_kcal_mol")

    # Save full results
    df.to_csv(output_csv, index=False)

    # Show summary
    print(f"\n=== Filtering Summary ===")
    print(f"Total docked: {len(df)}")
    print(f"Pass Lipinski: {df['passes_lipinski'].sum()}")

    # Top hits that pass Lipinski
    top_hits = df[df["passes_lipinski"]].head(10)
    print(f"\n=== Top 10 Drug-Like Hits ===")
    cols = ["ligand_id", "affinity_kcal_mol", "MW", "logP", "HBD", "HBA"]
    print(top_hits[cols].to_string(index=False))

    return df


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python filter_hits.py <docking_results.csv> <ligand_id_map.csv> <output_csv>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3])
