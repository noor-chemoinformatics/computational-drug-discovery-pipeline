"""
Plot deep docking iteration metrics.

Usage:
    python plot_deep_docking.py <iteration_log.csv> <output.png>
"""
import sys
import pandas as pd
import matplotlib.pyplot as plt


def plot_metrics(csv_file, output_png):
    df = pd.read_csv(csv_file)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # R² per iteration
    axes[0].plot(df["iteration"], df["r2"], 'o-', color='blue')
    axes[0].set_xlabel("Iteration")
    axes[0].set_ylabel("R²")
    axes[0].set_title("Model R² per Iteration")
    axes[0].grid(True, alpha=0.3)

    # RMSE per iteration
    axes[1].plot(df["iteration"], df["rmse"], 's-', color='red')
    axes[1].set_xlabel("Iteration")
    axes[1].set_ylabel("RMSE (kcal/mol)")
    axes[1].set_title("Model RMSE per Iteration")
    axes[1].grid(True, alpha=0.3)

    # Library size per iteration
    axes[2].plot(df["iteration"], df["library_size"], '^-', color='green')
    axes[2].set_xlabel("Iteration")
    axes[2].set_ylabel("Library Size")
    axes[2].set_title("Library Size per Iteration")
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_png, dpi=150)
    print(f"Plot saved to {output_png}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python plot_deep_docking.py <iteration_log.csv> <output.png>")
        sys.exit(1)
    plot_metrics(sys.argv[1], sys.argv[2])
