"""
Plot cross-validation fold R2 and RMSE distributions.

Usage:
    python plot_cv_results.py <cv_results.csv> <output.png>
"""
import sys
import pandas as pd
import matplotlib.pyplot as plt


def plot_cv(csv_file, output_png):
    df = pd.read_csv(csv_file)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].bar(range(1, len(df) + 1), df["r2"], color='steelblue', alpha=0.7)
    axes[0].axhline(df["r2"].mean(), color='red', linestyle='--',
                    label=f'Mean: {df["r2"].mean():.3f}')
    axes[0].set_xlabel("Fold")
    axes[0].set_ylabel("R2")
    axes[0].set_title("Cross-Validation R2 per Fold")
    axes[0].legend()

    axes[1].bar(range(1, len(df) + 1), df["rmse"], color='coral', alpha=0.7)
    axes[1].axhline(df["rmse"].mean(), color='red', linestyle='--',
                    label=f'Mean: {df["rmse"].mean():.3f}')
    axes[1].set_xlabel("Fold")
    axes[1].set_ylabel("RMSE (kcal/mol)")
    axes[1].set_title("Cross-Validation RMSE per Fold")
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(output_png, dpi=150)
    print(f"Plot saved to {output_png}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python plot_cv_results.py <cv_results.csv> <output.png>")
        sys.exit(1)
    plot_cv(sys.argv[1], sys.argv[2])
