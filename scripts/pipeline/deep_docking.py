"""
Deep Docking: ML-accelerated virtual screening.

Usage:
    python deep_docking.py <library.smi> <receptor.pdbqt> <output_dir>
                        --center_x X --center_y Y --center_z Z
                        [--sample_size 500] [--keep_fraction 0.1]
                        [--n_iterations 3]
"""
import sys
import os
import subprocess
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error


def smiles_to_fingerprint(smi, radius=2, n_bits=2048):
    """Convert SMILES to Morgan fingerprint."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    return np.array(fp)


def smiles_to_3d_pdbqt(smi, output_pdbqt):
    """Generate 3D conformer and convert to PDBQT."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return False

    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, randomSeed=42)
    if result != 0:
        return False

    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass

    tmp_sdf = output_pdbqt.replace('.pdbqt', '.sdf')
    writer = Chem.SDWriter(tmp_sdf)
    writer.write(mol)
    writer.close()

    cmd = ["mk_prepare_ligand.py", "-i", tmp_sdf, "-o", output_pdbqt]
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        os.remove(tmp_sdf)
        return True
    except subprocess.CalledProcessError:
        if os.path.exists(tmp_sdf):
            os.remove(tmp_sdf)
        return False


def dock_single(ligand_pdbqt, receptor_pdbqt, output_pdbqt,
                center, size=(30, 30, 30), exhaustiveness=4):
    """Dock a single ligand, return best affinity."""
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
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            return None
        with open(output_pdbqt, 'r') as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT"):
                    return float(line.split()[3])
    except Exception:
        return None
    return None


def dock_batch(smiles_list, receptor_pdbqt, work_dir, center,
                size=(30, 30, 30), exhaustiveness=4):
    """Dock a batch of SMILES. Returns list of (smiles, score)."""
    os.makedirs(work_dir, exist_ok=True)
    results = []

    for i, smi in enumerate(tqdm(smiles_list, desc="Docking batch")):
        lig_pdbqt = os.path.join(work_dir, f"lig_{i:05d}.pdbqt")
        out_pdbqt = os.path.join(work_dir, f"lig_{i:05d}_out.pdbqt")

        if not smiles_to_3d_pdbqt(smi, lig_pdbqt):
            results.append((smi, None))
            continue

        score = dock_single(lig_pdbqt, receptor_pdbqt, out_pdbqt,
                            center, size, exhaustiveness)
        results.append((smi, score))

    return results


def train_model(smiles_list, scores):
    """Train Random Forest on fingerprints."""
    X, y = [], []
    for smi, score in zip(smiles_list, scores):
        if score is None:
            continue
        fp = smiles_to_fingerprint(smi)
        if fp is None:
            continue
        X.append(fp)
        y.append(score)

    X = np.array(X)
    y = np.array(y)

    if len(X) < 10:
        print("Not enough data to train")
        return None, None, None

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )

    model = RandomForestRegressor(n_estimators=100, random_state=42, n_jobs=-1)
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)
    r2 = r2_score(y_test, y_pred)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))

    return model, r2, rmse


def predict_scores(model, smiles_list):
    """Predict docking scores for all SMILES."""
    X = []
    valid_smiles = []

    for smi in smiles_list:
        fp = smiles_to_fingerprint(smi)
        if fp is not None:
            X.append(fp)
            valid_smiles.append(smi)

    X = np.array(X)
    predictions = model.predict(X)

    return valid_smiles, predictions


def deep_docking(library_file, receptor_pdbqt, output_dir, center,
                  sample_size=500, keep_fraction=0.1, n_iterations=3,
                  exhaustiveness=4):
    """Run full deep docking pipeline."""
    os.makedirs(output_dir, exist_ok=True)

    with open(library_file, 'r') as f:
        all_smiles = [line.strip().split()[0] for line in f if line.strip()]

    print(f"Total library size: {len(all_smiles)}")

    current_smiles = all_smiles.copy()
    iteration_results = []

    for iteration in range(n_iterations):
        print(f"\n{'='*60}")
        print(f"ITERATION {iteration + 1}")
        print(f"{'='*60}")
        print(f"Current library size: {len(current_smiles)}")

        n_sample = min(sample_size, len(current_smiles))
        sampled = np.random.choice(current_smiles, n_sample, replace=False).tolist()

        work_dir = os.path.join(output_dir, f"iter_{iteration}")
        print(f"Docking {n_sample} compounds...")
        results = dock_batch(sampled, receptor_pdbqt, work_dir, center,
                              exhaustiveness=exhaustiveness)

        sampled_smiles = [r[0] for r in results if r[1] is not None]
        sampled_scores = [r[1] for r in results if r[1] is not None]

        print(f"Successfully docked: {len(sampled_scores)}/{n_sample}")

        if len(sampled_scores) < 10:
            print("Not enough data, stopping")
            break

        print("Training ML model...")
        model, r2, rmse = train_model(sampled_smiles, sampled_scores)
        print(f"Model R²: {r2:.3f}, RMSE: {rmse:.3f}")

        print("Predicting scores for remaining compounds...")
        remaining = [s for s in current_smiles if s not in set(sampled_smiles)]
        valid_remaining, predictions = predict_scores(model, remaining)

        all_valid = list(valid_remaining) + sampled_smiles
        all_preds = list(predictions) + sampled_scores

        sorted_idx = np.argsort(all_preds)
        n_keep = max(int(len(sorted_idx) * keep_fraction), 100)

        top_idx = sorted_idx[:n_keep]
        current_smiles = [all_valid[i] for i in top_idx]

        print(f"Kept top {n_keep} compounds (keep fraction: {keep_fraction})")

        iteration_results.append({
            "iteration": iteration + 1,
            "library_size": len(current_smiles),
            "sampled": n_sample,
            "r2": r2,
            "rmse": rmse,
        })

    final_file = os.path.join(output_dir, "final_smiles.smi")
    with open(final_file, 'w') as f:
        for smi in current_smiles:
            f.write(smi + "\n")

    print(f"\nFinal library size: {len(current_smiles)}")
    print(f"Final SMILES saved to {final_file}")

    log_df = pd.DataFrame(iteration_results)
    log_df.to_csv(os.path.join(output_dir, "iteration_log.csv"), index=False)

    return current_smiles


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("library_file")
    parser.add_argument("receptor_pdbqt")
    parser.add_argument("output_dir")
    parser.add_argument("--center_x", type=float, required=True)
    parser.add_argument("--center_y", type=float, required=True)
    parser.add_argument("--center_z", type=float, required=True)
    parser.add_argument("--sample_size", type=int, default=500)
    parser.add_argument("--keep_fraction", type=float, default=0.1)
    parser.add_argument("--n_iterations", type=int, default=3)
    parser.add_argument("--exhaustiveness", type=int, default=4)

    args = parser.parse_args()

    deep_docking(
        args.library_file,
        args.receptor_pdbqt,
        args.output_dir,
        center=(args.center_x, args.center_y, args.center_z),
        sample_size=args.sample_size,
        keep_fraction=args.keep_fraction,
        n_iterations=args.n_iterations,
        exhaustiveness=args.exhaustiveness,
    )
