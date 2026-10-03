"""
Analyze AutoDock Vina output.

Usage:
    python analyze_docking.py <vina_output.pdbqt> [output.csv]
"""
import sys
import pandas as pd


def parse_vina_output(output_pdbqt):
    """Parse Vina output file for binding affinities."""
    results = []
    with open(output_pdbqt, 'r') as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT"):
                parts = line.split()
                results.append({
                    "mode": len(results) + 1,
                    "affinity_kcal_mol": float(parts[3]),
                    "rmsd_lb": float(parts[4]),
                    "rmsd_ub": float(parts[5]),
                })
    return pd.DataFrame(results)


def main(output_pdbqt, output_csv=None):
    df = parse_vina_output(output_pdbqt)
    if df.empty:
        print("No results found.")
        return
    print("\n=== Docking Results ===")
    print(df.to_string(index=False))
    print(f"\nBest affinity: {df['affinity_kcal_mol'].min():.2f} kcal/mol")
    if output_csv:
        df.to_csv(output_csv, index=False)
        print(f"Saved to {output_csv}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python analyze_docking.py <vina_output.pdbqt> [output.csv]")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
