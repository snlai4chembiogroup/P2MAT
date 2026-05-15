######################################################
# Developer : Methun Kamruzzaman
# Date      : April 13, 2026
# Purpose   : Williams plots for applicability domain assessment of melting
#             point (2D) stacking ensemble models.
######################################################
"""
williams_plot.py
----------------
Williams plots (leverage vs. standardized residuals) for the three best
stacking ensemble models:

  cho  / whole / standard   (Standard Scaler)
  chon / whole / quantile   (Quantile Transformer)
  full / whole / standard   (Standard Scaler)

Workflow per model
------------------
1. Load the saved Ridge meta-learner from the stack artifact.
2. Load the cached OOF predictions (oof_{model}.npy) for each base model
   → stack into oof_matrix  [n_train × n_base]
   These are the same predictions used to compute train MAE in stacking_results.csv.
3. Load the cached test predictions (test_{model}.npy) for each base model
   → stack into test_matrix [n_test  × n_base]
4. y_tr = meta.predict(oof_matrix)    (OOF — matches CSV train MAE)
   y_ts = meta.predict(test_matrix)   (matches CSV test MAE)
5. Compute leverage in meta-feature space and standardized residuals.
6. Save a figure with two side-by-side subplots (train | test).
   Output: results/williams_plot_{dataset}.png
"""

import os
import sys
import warnings

import joblib
import numpy as np
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.config import DATASETS, MODELS_DIR, RESULTS_DIR
from src.data_loader import load_dataset
from src.stacking import BASE_MODELS, _load_oof_cache

# ---------------------------------------------------------------------------
# Best models (dataset, fs_method, transform)
# ---------------------------------------------------------------------------
BEST_MODELS = [
    ("cho",  "whole", "standard"),
    ("chon", "whole", "quantile"),
    ("full", "whole", "standard"),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _get_predictions(dataset, fs_method, transform):
    """
    Load cached OOF and test .npy files for every base model, stack into
    meta-feature matrices, then apply the saved Ridge meta-learner.

    Returns
    -------
    y_tr         : np.ndarray  OOF-based train predictions (consistent with CSV)
    y_ts         : np.ndarray  test predictions
    oof_matrix   : np.ndarray  [n_train × n_base]  — used for leverage
    test_matrix  : np.ndarray  [n_test  × n_base]  — used for leverage
    """
    save_dir = os.path.join(MODELS_DIR, dataset, fs_method, transform, "stack")

    oof_matrix  = None
    test_matrix = None

    for col, mn in enumerate(BASE_MODELS):
        cached = _load_oof_cache(save_dir, mn)
        if cached is None:
            raise FileNotFoundError(
                f"OOF cache missing for {mn} in {save_dir}. "
                "Re-run stack.py to regenerate."
            )
        oof, test_preds = cached

        if oof_matrix is None:
            oof_matrix  = np.zeros((len(oof),        len(BASE_MODELS)))
            test_matrix = np.zeros((len(test_preds), len(BASE_MODELS)))

        oof_matrix[:, col]  = oof
        test_matrix[:, col] = test_preds

    stack_path = os.path.join(save_dir, "pipeline.joblib")
    meta = joblib.load(stack_path)["meta_learner"]

    y_tr = meta.predict(oof_matrix)   # OOF — matches stacking_results.csv train MAE
    y_ts = meta.predict(test_matrix)  # matches stacking_results.csv test MAE

    return y_tr, y_ts, oof_matrix, test_matrix


def _compute_leverage(X_meta_tr, X_meta_query):
    """
    Leverage h_i = x_i^T (X_tr^T X_tr)^{-1} x_i for each row in X_meta_query,
    computed in the meta-feature space (base-model prediction space).
    """
    XtX_inv = np.linalg.pinv(X_meta_tr.T @ X_meta_tr)
    return np.einsum("ij,jk,ik->i", X_meta_query, XtX_inv, X_meta_query)


def _williams_panel(ax, h, std_res, h_star, color, marker, title, n):
    """Draw one Williams plot panel on ax."""
    ax.scatter(h, std_res, alpha=0.6, color=color, s=28, marker=marker,
               zorder=3, label=f"n = {n}")
    ax.axhline( 3, color="gray",  linestyle="--", linewidth=1, label="±3 std")
    ax.axhline(-3, color="gray",  linestyle="--", linewidth=1)
    ax.axvline(h_star, color="black", linestyle="--", linewidth=1,
               label=f"h* = {h_star:.3f}")
    ax.set_xlabel("Leverage (h)", fontsize=12)
    ax.set_ylabel("Standardized Residuals", fontsize=12)
    ax.set_title(title, fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)


# ---------------------------------------------------------------------------
# Main plot function
# ---------------------------------------------------------------------------

def make_williams_plot(dataset, fs_method, transform):
    ds_cfg = DATASETS[dataset]

    print(f"\nLoading {dataset.upper()} data ...")
    _, _, y_train, y_test, _ = load_dataset(ds_cfg["train"], ds_cfg["test"])

    y_train_arr = y_train.values if hasattr(y_train, "values") else np.asarray(y_train)
    y_test_arr  = y_test.values  if hasattr(y_test,  "values") else np.asarray(y_test)

    print(f"  Loading cached OOF / test predictions ...")
    y_tr, y_ts, oof_matrix, test_matrix = _get_predictions(
        dataset, fs_method, transform
    )

    # -----------------------------------------------------------------------
    # Residuals — standardize each split independently
    # -----------------------------------------------------------------------
    res_tr = y_train_arr - y_tr
    res_ts = y_test_arr  - y_ts

    std_tr = np.std(res_tr, ddof=1)
    std_ts = np.std(res_ts, ddof=1)

    std_res_tr = res_tr / std_tr
    std_res_ts = res_ts / std_ts

    # -----------------------------------------------------------------------
    # Leverage in meta-feature space
    # oof_matrix defines the training applicability domain
    # -----------------------------------------------------------------------
    h_tr = _compute_leverage(oof_matrix, oof_matrix)
    h_ts = _compute_leverage(oof_matrix, test_matrix)

    n_train = oof_matrix.shape[0]
    k       = oof_matrix.shape[1]   # number of base models
    h_star  = 3 * (k + 1) / n_train

    # -----------------------------------------------------------------------
    # Two side-by-side subplots: train (left) | test (right)
    # -----------------------------------------------------------------------
    fig, (ax_tr, ax_ts) = plt.subplots(1, 2, figsize=(13, 6))

    _williams_panel(
        ax_tr, h_tr, std_res_tr, h_star,
        color="steelblue", marker="o",
        title=f"Training Set  (n={len(y_train_arr)})",
        n=len(y_train_arr),
    )
    _williams_panel(
        ax_ts, h_ts, std_res_ts, h_star,
        color="tomato", marker="^",
        title=f"Test Set  (n={len(y_test_arr)})",
        n=len(y_test_arr),
    )

    fig.suptitle(
        f"Williams Plot — {dataset.upper()} / {fs_method} / {transform}  [stack]",
        fontsize=14, y=1.02,
    )

    out_path = os.path.join(RESULTS_DIR, f"williams_plot_{dataset}.png")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved → {out_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    for dataset, fs_method, transform in BEST_MODELS:
        make_williams_plot(dataset, fs_method, transform)


if __name__ == "__main__":
    main()
