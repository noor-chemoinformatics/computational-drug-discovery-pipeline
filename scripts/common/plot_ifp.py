"""
Plot interaction fingerprint heatmap.

Usage:
    python plot_ifp.py <ifp.csv> <output.png>
"""
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def plot_ifp(csv_file, output_png):
    df = pd.read_csv(csv_file, index_col=0)

    # Take first row (single ligand)
    if len(df) > 0:
        row = df.iloc[0]
        # Keep only interactions that are present
        active = row[row > 0]

        if len(active) == 0:
            print("No interactions found")
            return

        fig, ax = plt.subplots(figsize=(10, max(4, len(active) * 0.3)))
        active.plot(kind='barh', ax=ax, color='steelblue')
        ax.set_xlabel("Interaction strength")
        ax.set_title("Protein-Ligand Interactions (ProLIF)")
        plt.tight_layout()
        plt.savefig(output_png, dpi=150)
        print(f"Plot saved to {output_png}")
    else:
        print("Empty DataFrame")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python plot_ifp.py <ifp.csv> <output.png>")
        sys.exit(1)
    plot_ifp(sys.argv[1], sys.argv[2])
