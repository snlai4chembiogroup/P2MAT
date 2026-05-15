######################################################
# Developer : Methun Kamruzzaman
# Date      : May 14, 2026
# Purpose   : Computes structural similarity and information-leakage metrics
#             between two SMILES datasets.
######################################################
"""
dataset_similarity.py
---------------------
Compute two complementary metrics between two SMILES datasets to quantify
information leakage / structural overlap:

  MMTS  - Mean Maximum Tanimoto Similarity (continuous, threshold-free)
  LF    - Leakage Fraction: % of compounds with nearest-neighbour
           Tanimoto >= threshold (default 0.8; Butina 1999)

Both input CSVs must contain a column named 'smiles'.

Usage
-----
    python dataset_similarity.py datasetA.csv datasetB.csv
    python dataset_similarity.py datasetA.csv datasetB.csv --threshold 0.8 --output results.csv
"""

import argparse
import os
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from tqdm import tqdm

RDLogger.DisableLog('rdApp.*')

SIMILARITY_THRESHOLD = 0.8


# ---------------------------------------------------------------------------
# Fingerprints
# ---------------------------------------------------------------------------

def compute_fps(smiles_list, label=""):
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps, valid_idx = [], []
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(str(smi))
        if mol:
            fps.append(gen.GetFingerprint(mol))
            valid_idx.append(i)
    n_invalid = len(smiles_list) - len(fps)
    if n_invalid:
        print(f"  [{label}] Skipped {n_invalid} invalid SMILES")
    return fps, valid_idx


# ---------------------------------------------------------------------------
# MMTS and leakage fraction
# ---------------------------------------------------------------------------

def mmts_one_direction(fps_query, fps_ref, label=""):
    """For each query fp, find max Tanimoto to any ref fp."""
    max_sims = []
    for fp in tqdm(fps_query, desc=label):
        sims = DataStructs.BulkTanimotoSimilarity(fp, fps_ref)
        max_sims.append(max(sims) if sims else 0.0)
    return np.array(max_sims)


def compute_mmts(fps_a, fps_b):
    """Return (a_to_b array, b_to_a array) of per-compound max similarities."""
    a_to_b = mmts_one_direction(fps_a, fps_b, label="A->B (each A vs all B)")
    b_to_a = mmts_one_direction(fps_b, fps_a, label="B->A (each B vs all A)")
    return a_to_b, b_to_a


def compute_leakage_fraction(max_sims, threshold):
    """Fraction (0-100) of compounds with nearest-neighbour >= threshold."""
    return (max_sims >= threshold).mean() * 100


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def print_report(a_to_b, b_to_a, name_a, name_b, threshold):
    sym_mmts = (a_to_b.mean() + b_to_a.mean()) / 2.0
    leak_ab  = compute_leakage_fraction(a_to_b, threshold)
    leak_ba  = compute_leakage_fraction(b_to_a, threshold)
    sym_lf   = (leak_ab + leak_ba) / 2.0

    W = 65
    print()
    print("=" * W)
    print("  DATASET SIMILARITY REPORT  (Morgan ECFP4, Tanimoto)")
    print("=" * W)
    print(f"  Dataset A : {name_a}  ({len(a_to_b)} valid compounds)")
    print(f"  Dataset B : {name_b}  ({len(b_to_a)} valid compounds)")
    print(f"  Threshold : {threshold}")
    print("-" * W)
    print(f"  {'Metric':<40} {'A->B':>8}  {'B->A':>8}")
    print(f"  {'-'*40} {'-'*8}  {'-'*8}")
    print(f"  {'Mean max similarity (MMTS)':<40} {a_to_b.mean():>8.4f}  {b_to_a.mean():>8.4f}")
    print(f"  {'Median max similarity':<40} {np.median(a_to_b):>8.4f}  {np.median(b_to_a):>8.4f}")
    print(f"  {'Min max similarity':<40} {a_to_b.min():>8.4f}  {b_to_a.min():>8.4f}")
    print(f"  {'Max max similarity':<40} {a_to_b.max():>8.4f}  {b_to_a.max():>8.4f}")
    print(f"  {'Leakage fraction (LF, % >= threshold)':<40} {leak_ab:>7.2f}%  {leak_ba:>7.2f}%")
    print("-" * W)
    print(f"  Symmetric MMTS                  : {sym_mmts:.4f}")
    print(f"  Symmetric leakage fraction (LF) : {sym_lf:.2f}%")
    print("=" * W)

    # Auto-generated paragraph
    overlap = "substantial" if sym_lf > 10 else "limited"
    risk    = "notable"     if sym_lf > 10 else "negligible"
    print("\n  SUMMARY PARAGRAPH")
    print("-" * W)
    print(
        f"Structural similarity between {os.path.basename(name_a)} and "
        f"{os.path.basename(name_b)} was assessed using two complementary metrics "
        f"based on Morgan fingerprints (ECFP4, radius = 2, 2048 bits). "
        f"The Mean Maximum Tanimoto Similarity (MMTS) was computed by identifying, "
        f"for each compound, its nearest neighbour in the opposite dataset. "
        f"In the A->B direction the MMTS was {a_to_b.mean():.4f} "
        f"(median {np.median(a_to_b):.4f}) and in the B->A direction "
        f"{b_to_a.mean():.4f} (median {np.median(b_to_a):.4f}), "
        f"giving a symmetric MMTS of {sym_mmts:.4f}. "
        f"The leakage fraction (LF) -- the percentage of compounds whose "
        f"nearest-neighbour Tanimoto similarity meets or exceeds the near-duplicate "
        f"threshold of {threshold} (Butina, 1999) -- was {leak_ab:.2f}% (A->B) "
        f"and {leak_ba:.2f}% (B->A), with a symmetric LF of {sym_lf:.2f}%. "
        f"Together, these metrics indicate that the two datasets share "
        f"{overlap} structural overlap, with {risk} risk of information "
        f"leakage between them."
    )
    print("=" * W + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compute MMTS and leakage fraction between two SMILES datasets."
    )
    parser.add_argument("dataset_a", help="Path to first CSV (must have 'smiles' column)")
    parser.add_argument("dataset_b", help="Path to second CSV (must have 'smiles' column)")
    parser.add_argument(
        "--threshold", type=float, default=SIMILARITY_THRESHOLD,
        help=f"Tanimoto threshold for leakage fraction (default: {SIMILARITY_THRESHOLD})"
    )
    parser.add_argument(
        "--output", default=None,
        help="Optional path to save per-compound max-similarity scores as CSV"
    )
    args = parser.parse_args()

    # Load
    df_a = pd.read_csv(args.dataset_a)
    df_b = pd.read_csv(args.dataset_b)
    assert "smiles" in df_a.columns, "Dataset A has no 'smiles' column."
    assert "smiles" in df_b.columns, "Dataset B has no 'smiles' column."

    print(f"\nDataset A : {args.dataset_a}  ({len(df_a)} rows)")
    print(f"Dataset B : {args.dataset_b}  ({len(df_b)} rows)")

    # Fingerprints
    print("\nComputing fingerprints...")
    fps_a, valid_a = compute_fps(df_a["smiles"].tolist(), label="A")
    fps_b, valid_b = compute_fps(df_b["smiles"].tolist(), label="B")

    # MMTS
    print("\nComputing MMTS and leakage fraction...")
    a_to_b, b_to_a = compute_mmts(fps_a, fps_b)

    # Report
    print_report(a_to_b, b_to_a, args.dataset_a, args.dataset_b, args.threshold)

    # Optional output
    if args.output:
        out_a = df_a.iloc[valid_a].copy().reset_index(drop=True)
        out_b = df_b.iloc[valid_b].copy().reset_index(drop=True)
        out_a["max_sim_to_B"] = np.round(a_to_b, 6)
        out_b["max_sim_to_A"] = np.round(b_to_a, 6)
        base, ext = os.path.splitext(args.output)
        path_a = f"{base}_A{ext}"
        path_b = f"{base}_B{ext}"
        out_a.to_csv(path_a, index=False)
        out_b.to_csv(path_b, index=False)
        print(f"  Per-compound scores saved:")
        print(f"    {path_a}")
        print(f"    {path_b}\n")


if __name__ == "__main__":
    main()
