######################################################
# Developer : Methun Kamruzzaman
# Date      : May 6, 2026
# Purpose   : Williams-plot-based outlier detection for the best melting point
#             (2D) stacking ensemble models.
######################################################
"""
outlier_detection.py
--------------------
Williams-plot-based outlier detection for the top-N stacking models ranked by
test MAE (loaded dynamically from results/stacking_results.csv).

Workflow per dataset
--------------------
1. Load OOF (train) and test predictions via the stacking pipeline cache.
2. Compute leverage h and standardized residuals for both splits.
   - Leverage uses the OOF meta-feature matrix as the training-set centroid.
   - Standardized residuals for both splits are scaled by the TRAINING residual
     std so that the scale is consistent (standard AD practice).
3. Classify training-set compounds:
     - Response outlier:    |std_res| > 3
     - Leverage outlier:    h > h*  (outside applicability domain)
     - Bad leverage point:  both    (structurally influential + large residual)
4. Flag test-set compounds using the same thresholds (do NOT remove them).
5. Report test performance:
     - All test compounds
     - Excluding flagged test outliers
     NOTE: High-leverage test compounds are extrapolations beyond the AD.
6. Save annotated Williams plots → results/outlier_williams_{dataset}.png
7. Save training-set outlier SMILES → results/outliers_smiles_{dataset}.csv

Usage
-----
    python outlier_detection.py           # top 3 by default
    python outlier_detection.py --top 5   # any N
"""

import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.config import DATASETS, MODELS_DIR, RESULTS_DIR, DATA_DIR
from src.data_loader import load_dataset
from src.stacking import BASE_MODELS, _load_oof_cache
from src.evaluator import compute_metrics

# Re-use leverage helper from williams_plot (avoids code duplication)
from williams_plot import _get_predictions, _compute_leverage

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
STD_RES_THRESH = 3.0   # |standardized residual| cutoff


# ---------------------------------------------------------------------------
# Top-combination loader  (sourced from retrain_cleaned.py)
# ---------------------------------------------------------------------------

def load_top_combinations(n: int) -> list:
    """Return the top-N (dataset, fs_method, transform) tuples ranked by test MAE."""
    csv_path = os.path.join(RESULTS_DIR, "stacking_results.csv")
    if not os.path.exists(csv_path):
        raise FileNotFoundError(
            f"Stacking results not found: {csv_path}\n"
            "Run the original training pipeline first."
        )
    df = pd.read_csv(csv_path)
    stacks = df[df[np.logical_and(df["model"] == "stack", df["fs_method"] == 'whole')].sort_values("test_mae")  # type: ignore[call-overload]
    top = stacks.head(n)[["dataset", "fs_method", "transform"]]
    return [tuple(row) for row in top.itertuples(index=False, name=None)]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_smiles(dataset: str):
    """
    Load raw train/test CSVs and return SMILES as numpy arrays.
    Row order matches load_dataset (no row-level reordering in data_loader).
    """
    ds_cfg = DATASETS[dataset]
    train_raw = pd.read_csv(os.path.join(DATA_DIR, ds_cfg["train"]))
    test_raw  = pd.read_csv(os.path.join(DATA_DIR, ds_cfg["test"]))

    def _get(df):
        if "smiles" in df.columns:
            return df["smiles"].values
        return np.array(["N/A"] * len(df))

    return _get(train_raw), _get(test_raw)


def _outlier_categories(std_res, h, h_star):
    """
    Return boolean arrays for outlier classification.

    Returns
    -------
    high_res : |std_res| > STD_RES_THRESH
    high_lev : h > h_star
    bad_lev  : both high_res and high_lev
    any_flag : high_res OR high_lev
    """
    high_res = np.abs(std_res) > STD_RES_THRESH
    high_lev = h > h_star
    bad_lev  = high_res & high_lev
    any_flag = high_res | high_lev
    return high_res, high_lev, bad_lev, any_flag


def _save_williams_plot(
    dataset, transform,
    h_tr, std_res_tr, h_ts, std_res_ts, h_star,
    n_train, n_test,
    tr_high_res, tr_high_lev, ts_any_flag,
):
    """Williams plot with outlier categories colour-coded."""
    fig, (ax_tr, ax_ts) = plt.subplots(1, 2, figsize=(14, 6))

    # ---- Training panel ------------------------------------------------
    tr_bad_lev  = tr_high_res & tr_high_lev
    tr_res_only = tr_high_res & ~tr_high_lev
    tr_lev_only = ~tr_high_res & tr_high_lev
    tr_normal   = ~tr_high_res & ~tr_high_lev

    ax_tr.scatter(
        h_tr[tr_normal],   std_res_tr[tr_normal],
        alpha=0.5, color="steelblue", s=18, label="Normal",
    )
    ax_tr.scatter(
        h_tr[tr_res_only], std_res_tr[tr_res_only],
        alpha=0.85, color="orange", s=45, marker="s", label=f"|e|>{STD_RES_THRESH:.0f} only",
    )
    ax_tr.scatter(
        h_tr[tr_lev_only], std_res_tr[tr_lev_only],
        alpha=0.85, color="purple", s=45, marker="D", label="h>h* only",
    )
    ax_tr.scatter(
        h_tr[tr_bad_lev],  std_res_tr[tr_bad_lev],
        alpha=0.95, color="red", s=70, marker="*", zorder=5, label="Both (bad leverage pt.)",
    )

    ax_tr.axhline( STD_RES_THRESH, color="gray",  linestyle="--", linewidth=1)
    ax_tr.axhline(-STD_RES_THRESH, color="gray",  linestyle="--", linewidth=1,
                  label=f"±{STD_RES_THRESH:.0f} std")
    ax_tr.axvline(h_star,          color="black", linestyle="--", linewidth=1,
                  label=f"h* = {h_star:.4f}")
    ax_tr.set_xlabel("Leverage (h)", fontsize=12)
    ax_tr.set_ylabel("Standardized Residuals", fontsize=12)
    ax_tr.set_title(f"Training Set  (n={n_train})", fontsize=12)
    ax_tr.legend(fontsize=8)
    ax_tr.grid(True, alpha=0.3)

    # ---- Test panel ----------------------------------------------------
    ts_normal = ~ts_any_flag

    ax_ts.scatter(
        h_ts[ts_normal],  std_res_ts[ts_normal],
        alpha=0.5, color="steelblue", s=18, label="In AD",
    )
    ax_ts.scatter(
        h_ts[ts_any_flag], std_res_ts[ts_any_flag],
        alpha=0.85, color="red", s=55, marker="^", zorder=5, label="Flagged",
    )

    ax_ts.axhline( STD_RES_THRESH, color="gray",  linestyle="--", linewidth=1)
    ax_ts.axhline(-STD_RES_THRESH, color="gray",  linestyle="--", linewidth=1,
                  label=f"±{STD_RES_THRESH:.0f} std")
    ax_ts.axvline(h_star,          color="black", linestyle="--", linewidth=1,
                  label=f"h* = {h_star:.4f}")
    ax_ts.set_xlabel("Leverage (h)", fontsize=12)
    ax_ts.set_ylabel("Standardized Residuals", fontsize=12)
    ax_ts.set_title(f"Test Set  (n={n_test})", fontsize=12)
    ax_ts.legend(fontsize=8)
    ax_ts.grid(True, alpha=0.3)

    fig.suptitle(
        f"Williams Plot — Outlier Detection  |  {dataset.upper()} / {transform}",
        fontsize=13, y=1.01,
    )

    out_path = os.path.join(RESULTS_DIR, f"outlier_williams_{dataset}.png")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Williams plot (annotated) → {out_path}")


def _save_outlier_smiles(
    dataset, smiles_arr, any_flag, high_res, high_lev, y_true, y_pred, h, std_res
):
    """Save flagged training-set SMILES and diagnostics to CSV."""
    indices = np.where(any_flag)[0]
    if len(indices) == 0:
        print(f"  No training outliers to save for {dataset.upper()}.")
        return

    rows = []
    for i in indices:
        kind = []
        if high_res[i]:
            kind.append("response_outlier")
        if high_lev[i]:
            kind.append("leverage_outlier")
        rows.append({
            "index":        int(i),
            "smiles":       smiles_arr[i] if i < len(smiles_arr) else "N/A",
            "y_true_K":     round(float(y_true[i]), 3),
            "y_pred_K":     round(float(y_pred[i]), 3),
            "residual_K":   round(float(y_true[i] - y_pred[i]), 3),
            "std_residual": round(float(std_res[i]), 3),
            "leverage_h":   round(float(h[i]), 6),
            "outlier_type": "+".join(kind),
        })

    out_df   = pd.DataFrame(rows)
    out_path = os.path.join(RESULTS_DIR, f"outliers_smiles_{dataset}.csv")
    out_df.to_csv(out_path, index=False)
    print(f"  Training outlier SMILES  → {out_path}  ({len(rows)} compounds)")


# ---------------------------------------------------------------------------
# Main detection function
# ---------------------------------------------------------------------------

def detect_outliers(dataset: str, fs_method: str, transform: str) -> dict:
    """
    Full outlier detection workflow for one (dataset, fs_method, transform).

    Returns a summary dict for the final report table.
    """
    ds_cfg = DATASETS[dataset]

    print(f"\n{'='*70}")
    print(f"  Outlier Detection — {dataset.upper()} / {fs_method} / {transform}")
    print(f"{'='*70}")

    # ---- Load true labels --------------------------------------------------
    _, _, y_train, y_test, _ = load_dataset(ds_cfg["train"], ds_cfg["test"])
    y_train_arr = y_train.values if hasattr(y_train, "values") else np.asarray(y_train)
    y_test_arr  = y_test.values  if hasattr(y_test,  "values") else np.asarray(y_test)

    # ---- Load SMILES from raw CSV ------------------------------------------
    train_smiles, _ = _load_smiles(dataset)

    # ---- Load stacking predictions -----------------------------------------
    print(f"  Loading cached OOF / test predictions ...")
    y_tr, y_ts, oof_matrix, test_matrix = _get_predictions(dataset, fs_method, transform)

    # ---- Residuals ---------------------------------------------------------
    res_tr = y_train_arr - y_tr
    res_ts = y_test_arr  - y_ts

    # Scale BOTH splits by training residual std (consistent AD scale)
    std_tr     = np.std(res_tr, ddof=1)
    std_res_tr = res_tr / std_tr
    std_res_ts = res_ts / std_tr   # intentional: training std as reference

    # ---- Leverage ----------------------------------------------------------
    h_tr = _compute_leverage(oof_matrix, oof_matrix)
    h_ts = _compute_leverage(oof_matrix, test_matrix)

    n_train = oof_matrix.shape[0]
    k       = oof_matrix.shape[1]   # number of base models
    h_star  = 3 * (k + 1) / n_train

    # ---- Classify outliers -------------------------------------------------
    tr_high_res, tr_high_lev, tr_bad_lev, tr_any = _outlier_categories(
        std_res_tr, h_tr, h_star
    )
    ts_high_res, ts_high_lev, ts_bad_lev, ts_any = _outlier_categories(
        std_res_ts, h_ts, h_star
    )

    # ---- Print training summary --------------------------------------------
    print(f"\n  Training set  (n={n_train}, h* = {h_star:.4f}):")
    print(f"    Response outliers   |std_res| > {STD_RES_THRESH:.0f} : {tr_high_res.sum():>4}")
    print(f"    Leverage outliers   h > h*                  : {tr_high_lev.sum():>4}")
    print(f"    Bad leverage pts    both                    : {tr_bad_lev.sum():>4}")
    print(f"    Total flagged                               : {tr_any.sum():>4}")

    # ---- Test set metrics --------------------------------------------------
    m_all   = compute_metrics(y_test_arr, y_ts)
    ts_keep = ~ts_any
    m_clean = (
        compute_metrics(y_test_arr[ts_keep], y_ts[ts_keep])
        if ts_keep.sum() > 0 else None
    )

    print(f"\n  Test set  (n={len(y_test_arr)}):")
    print(f"    Flagged (|std_res|>{STD_RES_THRESH:.0f} or h>h*)  : {ts_any.sum():>4}")
    print(f"    Bad leverage pts (both)             : {ts_bad_lev.sum():>4}")

    print(f"\n  Test performance — ALL compounds (n={len(y_test_arr)}):")
    print(f"    MAE  = {m_all['mae']:.3f}  K")
    print(f"    R²   = {m_all['r2']:.4f}")
    print(f"    RMSE = {m_all['rmse']:.3f}  K")

    if m_clean is not None and ts_any.sum() > 0:
        print(
            f"\n  Test performance — EXCLUDING {ts_any.sum()} flagged outliers"
            f"  (n={ts_keep.sum()}):"
        )
        print(f"    MAE  = {m_clean['mae']:.3f}  K")
        print(f"    R²   = {m_clean['r2']:.4f}")
        print(f"    RMSE = {m_clean['rmse']:.3f}  K")
        if ts_high_lev.sum() > 0:
            print(
                f"\n  ** NOTE: {ts_high_lev.sum()} test compound(s) have h > h* = {h_star:.4f}."
                f"\n     These lie OUTSIDE the training applicability domain."
                f"\n     Their predictions are extrapolations and should be treated"
                f"\n     with caution — they are flagged but NOT removed from the"
                f"\n     test set."
            )

    # ---- Save annotated Williams plot --------------------------------------
    _save_williams_plot(
        dataset, transform,
        h_tr, std_res_tr, h_ts, std_res_ts, h_star,
        n_train, len(y_test_arr),
        tr_high_res, tr_high_lev, ts_any,
    )

    # ---- Save training outlier SMILES --------------------------------------
    _save_outlier_smiles(
        dataset, train_smiles,
        tr_any, tr_high_res, tr_high_lev,
        y_train_arr, y_tr, h_tr, std_res_tr,
    )

    return {
        "dataset":               dataset,
        "n_train":               n_train,
        "n_test":                len(y_test_arr),
        "h_star":                round(h_star, 5),
        "tr_response_outliers":  int(tr_high_res.sum()),
        "tr_leverage_outliers":  int(tr_high_lev.sum()),
        "tr_bad_leverage_pts":   int(tr_bad_lev.sum()),
        "tr_total_flagged":      int(tr_any.sum()),
        "ts_flagged":            int(ts_any.sum()),
        "ts_bad_leverage_pts":   int(ts_bad_lev.sum()),
        "test_mae_all":          m_all["mae"],
        "test_r2_all":           m_all["r2"],
        "test_rmse_all":         m_all["rmse"],
        "test_mae_clean":        m_clean["mae"]  if m_clean else m_all["mae"],
        "test_r2_clean":         m_clean["r2"]   if m_clean else m_all["r2"],
        "test_rmse_clean":       m_clean["rmse"] if m_clean else m_all["rmse"],
    }


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--top", type=int, default=3,
                        help="Number of top combinations to analyse (default: 3)")
    args = parser.parse_args()

    best_models = load_top_combinations(args.top)

    print(f"\n{'='*70}")
    print(f"  OUTLIER DETECTION  (top {args.top} by test MAE)")
    print(f"{'='*70}")
    for i, (ds, fs, tr) in enumerate(best_models, 1):
        print(f"    {i}. {ds.upper()} / {fs} / {tr}")

    summary = []
    for dataset, fs_method, transform in best_models:
        res = detect_outliers(dataset, fs_method, transform)
        summary.append(res)

    # ---- Summary table -------------------------------------------------------
    print(f"\n\n{'='*100}")
    print("  OUTLIER DETECTION SUMMARY")
    print(f"{'='*100}")
    hdr = (
        f"{'Dataset':<8} {'n_tr':>6} {'n_ts':>6} {'h*':>8} "
        f"{'Tr|e|>3':>8} {'Tr_h>h*':>8} {'Tr_both':>8} {'Ts_flag':>8} "
        f"{'MAE_all':>9} {'MAE_clean':>10} {'R2_all':>8} {'R2_clean':>9}"
    )
    print(hdr)
    print("-" * 100)
    for r in summary:
        print(
            f"{r['dataset'].upper():<8} {r['n_train']:>6} {r['n_test']:>6} "
            f"{r['h_star']:>8.5f} "
            f"{r['tr_response_outliers']:>8} {r['tr_leverage_outliers']:>8} "
            f"{r['tr_bad_leverage_pts']:>8} {r['ts_flagged']:>8} "
            f"{r['test_mae_all']:>9.3f} {r['test_mae_clean']:>10.3f} "
            f"{r['test_r2_all']:>8.4f} {r['test_r2_clean']:>9.4f}"
        )
    print(f"{'='*100}")


if __name__ == "__main__":
    main()
