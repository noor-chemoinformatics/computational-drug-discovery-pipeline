"""
Calculate speed-up of deep docking vs brute force.

Usage:
    python calculate_speedup.py <iteration_log.csv> <library_size>
"""
import sys
import pandas as pd


def calculate_speedup(iteration_log_csv, library_size):
    df = pd.read_csv(iteration_log_csv)
    total_docked = df["sampled"].sum()
    speedup = library_size / total_docked if total_docked > 0 else 0

    print(f"Library size: {library_size}")
    print(f"Total docked: {total_docked}")
    print(f"Speed-up: {speedup:.2f}x")
    print(f"Compounds never docked: {library_size - total_docked}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python calculate_speedup.py <iteration_log.csv> <library_size>")
        sys.exit(1)
    calculate_speedup(sys.argv[1], int(sys.argv[2]))
