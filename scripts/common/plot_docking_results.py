"""
Plot docking results distribution.

Usage:
    python plot_docking_results.py <results.csv> <output.png>
"""
import sys
import pandas as pd
import matplotlib.pyplot as plt


def plot_results(csv_file, output_png):
    df = pd.read_csv(csv_file)
    df = df.dropna(subset=["affinity_kcal_mol"])

    # Handle both 'ligand' and 'ligand_id' column names
    if "ligand_id" not in df.columns and "ligand" in df.columns:
        df["ligand_id"] = df["ligand"].str.replace(".pdbqt", "", regex=False)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Histogram of affinities
    axes[0].hist(df["affinity_kcal_mol"], bins=30, edgecolor='black', alpha=0.7)
    axes[0].axvline(df["affinity_kcal_mol"].min(), color='red',
                    linestyle='--', label=f'Best: {df["affinity_kcal_mol"].min():.2f}')
    axes[0].set_xlabel("Binding Affinity (kcal/mol)")
    axes[0].set_ylabel("Number of Compounds")
    axes[0].set_title("Distribution of Docking Scores")
    axes[0].legend()

    # Top 20 bar chart
    top20 = df.nsmallest(20, "affinity_kcal_mol")
    axes[1].barh(range(len(top20)), top20["affinity_kcal_mol"])
    axes[1].set_yticks(range(len(top20)))
    axes[1].set_yticklabels(top20["ligand_id"], fontsize=8)
    axes[1].set_xlabel("Binding Affinity (kcal/mol)")
    axes[1].set_title("Top 20 Compounds")
    axes[1].invert_yaxis()

    plt.tight_layout()
    plt.savefig(output_png, dpi=150)
    print(f"Plot saved to {output_png}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python plot_docking_results.py <results.csv> <output.png>")
        sys.exit(1)
    plot_results(sys.argv[1], sys.argv[2])
