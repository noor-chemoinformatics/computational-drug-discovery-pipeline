"""
Plot interaction contacts from the simple IFP analysis.

Usage:
    python plot_interaction_contacts.py <contacts.csv> <output.png>
"""
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def plot_contacts(csv_file, output_png):
    df = pd.read_csv(csv_file)

    if len(df) == 0:
        print("No contacts found")
        return

    # Count contacts per residue
    residue_counts = df['protein_residue'].value_counts().head(20)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Bar chart of top residues
    residue_counts.plot(kind='barh', ax=axes[0], color='steelblue')
    axes[0].set_xlabel("Number of contacts")
    axes[0].set_ylabel("Residue")
    axes[0].set_title("Top Interacting Residues (within 4 A)")
    axes[0].invert_yaxis()

    # Histogram of distances
    axes[1].hist(df['distance'], bins=20, color='coral', edgecolor='black', alpha=0.7)
    axes[1].axvline(4.0, color='red', linestyle='--', label='Cutoff (4 A)')
    axes[1].set_xlabel("Distance (A)")
    axes[1].set_ylabel("Number of contacts")
    axes[1].set_title("Contact Distance Distribution")
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(output_png, dpi=150)
    print(f"Plot saved to {output_png}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python plot_interaction_contacts.py <contacts.csv> <output.png>")
        sys.exit(1)
    plot_contacts(sys.argv[1], sys.argv[2])
