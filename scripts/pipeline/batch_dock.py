"""
Batch docking with AutoDock Vina.

Usage:
    python batch_dock.py <ligand_dir> <receptor.pdbqt> <output_csv>
                        --center_x X --center_y Y --center_z Z
"""
import sys
import os
import subprocess
import argparse
import pandas as pd
from tqdm import tqdm


def dock_single_ligand(ligand_pdbqt, receptor_pdbqt, output_pdbqt,
                        center, size=(30, 30, 30), exhaustiveness=8):
    """Run Vina for a single ligand. Returns best affinity or None."""
    cmd = [
        "vina",
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--out", output_pdbqt,
        "--exhaustiveness", str(exhaustiveness),
        "--seed", "42",
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            return None

        # Parse affinity from output PDBQT
        with open(output_pdbqt, 'r') as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT"):
                    parts = line.split()
                    return float(parts[3])

    except subprocess.TimeoutExpired:
        return None
    except Exception as e:
        print(f"Error: {e}")
        return None

    return None


def batch_dock(ligand_dir, receptor_pdbqt, output_csv,
                center, size=(30, 30, 30), exhaustiveness=8):
    """Dock all ligands in a directory."""
    ligand_files = sorted([
        f for f in os.listdir(ligand_dir)
        if f.endswith('.pdbqt')
    ])

    print(f"Found {len(ligand_files)} ligands to dock")

    output_dir = "data/processed/docked"
    os.makedirs(output_dir, exist_ok=True)

    results = []

    for lig_file in tqdm(ligand_files, desc="Docking"):
        ligand_path = os.path.join(ligand_dir, lig_file)
        output_path = os.path.join(output_dir, lig_file.replace('.pdbqt', '_out.pdbqt'))

        affinity = dock_single_ligand(
            ligand_path, receptor_pdbqt, output_path,
            center, size, exhaustiveness
        )

        results.append({
            "ligand": lig_file,
            "affinity_kcal_mol": affinity,
        })

    df = pd.DataFrame(results)
    df = df.sort_values('affinity_kcal_mol', na_position='last')
    df.to_csv(output_csv, index=False)

    print(f"\nResults saved to {output_csv}")
    print(f"Successfully docked: {df['affinity_kcal_mol'].notna().sum()}/{len(df)}")
    print(f"\nTop 10:")
    print(df.head(10).to_string(index=False))

    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("ligand_dir")
    parser.add_argument("receptor_pdbqt")
    parser.add_argument("output_csv")
    parser.add_argument("--center_x", type=float, required=True)
    parser.add_argument("--center_y", type=float, required=True)
    parser.add_argument("--center_z", type=float, required=True)
    parser.add_argument("--exhaustiveness", type=int, default=8)

    args = parser.parse_args()

    batch_dock(
        args.ligand_dir,
        args.receptor_pdbqt,
        args.output_csv,
        center=(args.center_x, args.center_y, args.center_z),
        exhaustiveness=args.exhaustiveness,
    )
