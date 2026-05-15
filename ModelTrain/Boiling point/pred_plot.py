######################################################
# Developer : Methun Kamruzzaman
# Date      : April 13, 2026
# Purpose   : Real vs. predicted scatter plots for the best boiling point
#             stacking ensemble models.
######################################################
"""
pred_plot.py
------------
Real vs. predicted scatter plots for the three best stacking ensemble models:

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
5. Save a figure with two side-by-side subplots (train | test), each showing
   actual vs. predicted with a dashed diagonal perfect-prediction line.
   Output: results/pred_plot_{dataset}.png
"""

import os
import sys
import warnings

import joblib
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_absolute_error, r2_score

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
# Helper: load cached OOF / test predictions and feed through meta-learner
# ---------------------------------------------------------------------------

def _get_predictions(dataset, fs_method, transform):
    """
    Load the cached OOF and test .npy files for every base model, stack them
    into meta-feature matrices, then apply the saved Ridge meta-learner.

    Returns
    -------
    y_tr : np.ndarray   OOF-based train predictions (consistent with CSV MAE)
    y_ts : np.ndarray   test predictions
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

    return y_tr, y_ts


# ---------------------------------------------------------------------------
# Plot helper
# ---------------------------------------------------------------------------

def _scatter_panel(ax, y_pred, y_true, color, marker, title):
    """Draw one actual-vs-predicted panel on ax."""
    mae = mean_absolute_error(y_true, y_pred)
    r2  = r2_score(y_true, y_pred)

    all_vals = np.concatenate([y_true, y_pred])
    lo = all_vals.min()
    hi = all_vals.max()
    pad  = (hi - lo) * 0.04
    diag = np.array([lo - pad, hi + pad])

    ax.scatter(y_pred, y_true, alpha=0.6, color=color, s=28,
               marker=marker, zorder=3)
    ax.plot(diag, diag, color="black", linestyle="--", linewidth=1.2,
            zorder=2, label="Perfect prediction")

    ax.set_xlim(diag)
    ax.set_ylim(diag)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Predicted (K)", fontsize=12)
    ax.set_ylabel("Actual (K)", fontsize=12)
    ax.set_title(title, fontsize=12)
    ax.grid(True, alpha=0.3)

    # Annotation box
    ax.text(
        0.04, 0.96,
        f"MAE = {mae:.2f} K\n$R^2$ = {r2:.4f}",
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.8),
    )


# ---------------------------------------------------------------------------
# Main plot function
# ---------------------------------------------------------------------------

def make_pred_plot(dataset, fs_method, transform):
    ds_cfg = DATASETS[dataset]

    print(f"\nLoading {dataset.upper()} data ...")
    _, _, y_train, y_test, _ = load_dataset(ds_cfg["train"], ds_cfg["test"])

    y_train_arr = y_train.values if hasattr(y_train, "values") else np.asarray(y_train)
    y_test_arr  = y_test.values  if hasattr(y_test,  "values") else np.asarray(y_test)

    print(f"  Loading cached OOF / test predictions ...")
    y_tr, y_ts = _get_predictions(dataset, fs_method, transform)

    # Verify consistency with stored metrics (informational)
    mae_tr = mean_absolute_error(y_train_arr, y_tr)
    mae_ts = mean_absolute_error(y_test_arr,  y_ts)
    print(f"  MAE — train (OOF): {mae_tr:.3f}   test: {mae_ts:.3f}")

    # -----------------------------------------------------------------------
    # Two side-by-side subplots: train (left) | test (right)
    # -----------------------------------------------------------------------
    fig, (ax_tr, ax_ts) = plt.subplots(1, 2, figsize=(13, 6))

    _scatter_panel(
        ax_tr, y_tr, y_train_arr,
        color="steelblue", marker="o",
        title=f"Training Set  (n={len(y_train_arr)})",
    )
    _scatter_panel(
        ax_ts, y_ts, y_test_arr,
        color="tomato", marker="o",
        title=f"Test Set  (n={len(y_test_arr)})",
    )

    fig.suptitle(
        f"Real vs. Predicted — {dataset.upper()} / {fs_method} / {transform}  [stack]",
        fontsize=14, y=1.02,
    )

    out_path = os.path.join(RESULTS_DIR, f"pred_plot_{dataset}.png")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved → {out_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    for dataset, fs_method, transform in BEST_MODELS:
        make_pred_plot(dataset, fs_method, transform)


if __name__ == "__main__":
    main()
