"""
Plot MD analysis results.

Usage:
    python plot_md_analysis.py <output_dir>
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt


def read_xvg(filepath):
    data = []
    with open(filepath) as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                data.append([float(parts[0]), float(parts[1])])
    return np.array(data)


def plot_md(output_dir):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    files = {
        "rmsd.xvg": ("Time (ns)", "RMSD (nm)", "Backbone RMSD", 'blue', axes[0, 0]),
        "rmsf.xvg": ("Residue", "RMSF (nm)", "Per-Residue RMSF", 'red', axes[0, 1]),
        "gyrate.xvg": ("Time (ps)", "Rg (nm)", "Radius of Gyration", 'green', axes[1, 0]),
        "hbonds.xvg": ("Time (ps)", "H-bonds", "Hydrogen Bonds", 'purple', axes[1, 1]),
    }

    for filename, (xlabel, ylabel, title, color, ax) in files.items():
        filepath = os.path.join(output_dir, filename)
        if os.path.exists(filepath):
            data = read_xvg(filepath)
            ax.plot(data[:, 0], data[:, 1], color=color)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(title)
            ax.grid(True, alpha=0.3)
        else:
            ax.text(0.5, 0.5, f"{filename} not found",
                    ha='center', va='center', transform=ax.transAxes)

    plt.tight_layout()
    output_png = os.path.join(output_dir, "md_analysis.png")
    plt.savefig(output_png, dpi=150)
    print(f"Plot saved to {output_png}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plot_md_analysis.py <output_dir>")
        sys.exit(1)
    plot_md(sys.argv[1])
