"""
Select top drug-like hits.

Usage:
    python select_top_hits.py <filtered_hits.csv> <output_csv> [n_top]
"""
import sys
import pandas as pd


def select_top_hits(input_csv, output_csv, n_top=20):
    df = pd.read_csv(input_csv)
    df = df[df["passes_lipinski"] == True]
    df = df.sort_values("affinity_kcal_mol")
    top = df.head(n_top)
    cols = ["ligand_id", "smiles", "affinity_kcal_mol",
            "MW", "logP", "HBD", "HBA", "TPSA", "RotBonds"]
    top = top[cols]
    top.to_csv(output_csv, index=False)
    print(f"Selected top {len(top)} hits")
    print(top[["ligand_id", "affinity_kcal_mol", "MW", "logP"]].to_string(index=False))
    return top


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python select_top_hits.py <filtered_hits.csv> <output_csv> [n_top]")
        sys.exit(1)
    n_top = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    select_top_hits(sys.argv[1], sys.argv[2], n_top)
