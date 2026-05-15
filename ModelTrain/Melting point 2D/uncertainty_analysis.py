######################################################
# Developer : Methun Kamruzzaman
# Date      : April 17, 2026
# Purpose   : Prediction uncertainty estimation on the test set using melting
#             point (2D) stacking ensembles.
######################################################
"""
uncertainty_analysis.py
------------------------
Prediction uncertainty analysis on the test set using the original
stacking ensemble models (no outlier removal).

Two complementary uncertainty sources are estimated:

  1. Base-model disagreement  — std of the 4 base-model test predictions.
     Captures how much the individual learners (LightGBM, XGBoost, ExtraTrees,
     DNN) disagree for each test sample.  High disagreement → high uncertainty.

  2. Conformal prediction intervals  — split-conformal approach using OOF
     residuals as a calibration set.
       • Non-conformity score : |y_train – stack_oof_pred|
       • For a target coverage level α, threshold q_α = quantile(scores, α)
       • Test interval : [ŷ_test − q_α , ŷ_test + q_α]
       • Empirical coverage = fraction of test samples inside the interval

  Plots saved to results/ :
    uncertainty_analysis_{dataset}.png  (4-panel diagnostic)

Models used
-----------
  cho  / whole / standard
  chon / whole / quantile
  full / whole / standard
"""

import os
import sys
import warnings

import joblib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, r2_score

warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# ---------------------------------------------------------------------------
# Paths — use original (non-cleaned) models and results
# ---------------------------------------------------------------------------
BASE_DIR    = os.path.dirname(os.path.abspath(__file__))
MODELS_DIR  = os.path.join(BASE_DIR, "saved_models")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

import src.config as _cfg
_cfg.MODELS_DIR  = MODELS_DIR
_cfg.RESULTS_DIR = RESULTS_DIR

from src.config import DATASETS
from src.data_loader import load_dataset
from src.stacking import BASE_MODELS, _load_oof_cache

BEST_MODELS = [
    ("cho",  "whole", "standard"),
    ("chon", "whole", "quantile"),
    ("full", "whole", "standard"),
]

COVERAGE_LEVELS = [0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99]


# ---------------------------------------------------------------------------
# Prediction loader
# ---------------------------------------------------------------------------

def _load_predictions(dataset: str, fs_method: str, transform: str):
    """
    Load cached OOF and test predictions from the original model directory.

    Returns
    -------
    oof_matrix   : (n_train, n_base)  – OOF predictions per base model
    test_matrix  : (n_test,  n_base)  – test predictions per base model
    stack_oof    : (n_train,)         – stack meta-learner OOF predictions
    stack_test   : (n_test,)          – stack meta-learner test predictions
    meta_learner : fitted Ridge       – meta-learner object
    """
    save_dir = os.path.join(MODELS_DIR, dataset, fs_method, transform, "stack")

    oof_matrix  = None
    test_matrix = None

    for col, mn in enumerate(BASE_MODELS):
        cached = _load_oof_cache(save_dir, mn)
        if cached is None:
            raise FileNotFoundError(
                f"OOF cache missing for {mn} in {save_dir}. "
                "Re-run main.py / stack.py to regenerate."
            )
        oof, test_preds = cached

        if oof_matrix is None:
            oof_matrix  = np.zeros((len(oof),        len(BASE_MODELS)))
            test_matrix = np.zeros((len(test_preds), len(BASE_MODELS)))

        oof_matrix[:, col]  = oof
        test_matrix[:, col] = test_preds

    stack_path = os.path.join(save_dir, "pipeline.joblib")
    meta = joblib.load(stack_path)["meta_learner"]

    stack_oof  = meta.predict(oof_matrix)
    stack_test = meta.predict(test_matrix)

    return oof_matrix, test_matrix, stack_oof, stack_test, meta


# ---------------------------------------------------------------------------
# Conformal prediction intervals
# ---------------------------------------------------------------------------

def conformal_intervals(
    y_calib: np.ndarray,
    oof_preds: np.ndarray,
    test_preds: np.ndarray,
    coverage: float = 0.90,
):
    """
    Split-conformal prediction intervals.

    Parameters
    ----------
    y_calib    : true labels for the calibration set (full OOF training labels)
    oof_preds  : stack OOF predictions (same length as y_calib)
    test_preds : stack test predictions
    coverage   : desired coverage level (e.g. 0.90)

    Returns
    -------
    lo, hi : lower and upper bounds for each test point
    q      : the conformal quantile used
    """
    scores  = np.abs(y_calib - oof_preds)
    n       = len(scores)
    q_level = np.ceil((n + 1) * coverage) / n
    q_level = min(q_level, 1.0)
    q  = np.quantile(scores, q_level)
    lo = test_preds - q
    hi = test_preds + q
    return lo, hi, q


def coverage_calibration(
    y_calib: np.ndarray,
    oof_preds: np.ndarray,
    y_test: np.ndarray,
    test_preds: np.ndarray,
    levels: list = COVERAGE_LEVELS,
):
    nominal   = []
    empirical = []
    widths    = []

    for cov in levels:
        lo, hi, _ = conformal_intervals(y_calib, oof_preds, test_preds, coverage=cov)
        inside = ((y_test >= lo) & (y_test <= hi)).mean()
        nominal.append(cov)
        empirical.append(inside)
        widths.append((hi - lo).mean())

    return nominal, empirical, widths


# ---------------------------------------------------------------------------
# Main analysis + plot
# ---------------------------------------------------------------------------

def analyse_dataset(dataset: str, fs_method: str, transform: str):
    ds_cfg = DATASETS[dataset]
    print(f"\n{'='*60}")
    print(f"  Dataset: {dataset.upper()} / {fs_method} / {transform}")
    print(f"{'='*60}")

    _, _, y_train, y_test, _ = load_dataset(ds_cfg["train"], ds_cfg["test"])

    y_train_arr = y_train.values if hasattr(y_train, "values") else np.asarray(y_train)
    y_test_arr  = y_test.values  if hasattr(y_test,  "values") else np.asarray(y_test)

    print(f"  Train set: {len(y_train_arr)} samples (no outlier removal)")

    oof_matrix, test_matrix, stack_oof, stack_test, meta = \
        _load_predictions(dataset, fs_method, transform)

    n_test = len(y_test_arr)

    # -----------------------------------------------------------------------
    # 1. Base-model disagreement (uncertainty proxy)
    # -----------------------------------------------------------------------
    base_std = test_matrix.std(axis=1)

    # -----------------------------------------------------------------------
    # 2. Conformal prediction at 90 % coverage
    # -----------------------------------------------------------------------
    lo_90, hi_90, q_90 = conformal_intervals(
        y_train_arr, stack_oof, stack_test, coverage=0.90
    )
    inside_90 = ((y_test_arr >= lo_90) & (y_test_arr <= hi_90)).mean()

    # -----------------------------------------------------------------------
    # 3. Coverage calibration across multiple levels
    # -----------------------------------------------------------------------
    nominal, empirical, widths = coverage_calibration(
        y_train_arr, stack_oof, y_test_arr, stack_test
    )

    # -----------------------------------------------------------------------
    # 4. Metrics
    # -----------------------------------------------------------------------
    abs_errors = np.abs(y_test_arr - stack_test)
    rho, p_val = spearmanr(base_std, abs_errors)

    mae = mean_absolute_error(y_test_arr, stack_test)
    r2  = r2_score(y_test_arr, stack_test)

    print(f"  Stack MAE = {mae:.3f} K   R² = {r2:.4f}")
    print(f"  Conformal q(90%) = ±{q_90:.2f} K  |  empirical coverage = {inside_90:.3f}")
    print(f"  Spearman ρ (base_std vs |error|) = {rho:.4f}  (p={p_val:.4f})")

    print(f"\n  {'Nominal':>8}  {'Empirical':>10}  {'Interval width':>16}")
    print(f"  {'-'*38}")
    for nom, emp, w in zip(nominal, empirical, widths):
        flag = "  ✓" if abs(emp - nom) < 0.05 else ("  ↑" if emp > nom else "  ↓")
        print(f"  {nom:>7.0%}  {emp:>10.3f}  {w:>14.2f} K{flag}")

    # -----------------------------------------------------------------------
    # 5. Figure  (2 × 2 layout)
    # -----------------------------------------------------------------------
    fig = plt.figure(figsize=(15, 12))
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.32)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    title_kw = dict(fontsize=12, fontweight="bold")
    label_kw = dict(fontsize=10)

    # ── Panel 1: Base-model disagreement vs absolute error ──────────────────
    sc = ax1.scatter(
        base_std, abs_errors,
        c=stack_test, cmap="plasma", alpha=0.55, s=22, zorder=3,
    )
    plt.colorbar(sc, ax=ax1, label="Stack prediction (K)")
    ax1.set_xlabel("Base-model std (K)", **label_kw)
    ax1.set_ylabel("|True − Predicted| (K)", **label_kw)
    ax1.set_title("Base-Model Disagreement vs Absolute Error", **title_kw)
    ax1.text(
        0.97, 0.03,
        f"Spearman ρ = {rho:.3f}\n(p = {p_val:.3g})",
        transform=ax1.transAxes, fontsize=9,
        ha="right", va="bottom",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
    )
    ax1.grid(True, alpha=0.25)

    # ── Panel 2: Prediction intervals (sorted by prediction) ────────────────
    sort_idx = np.argsort(stack_test)
    x_plot   = np.arange(n_test)

    lo_80, hi_80, _ = conformal_intervals(
        y_train_arr, stack_oof, stack_test, coverage=0.80
    )

    ax2.fill_between(
        x_plot,
        lo_90[sort_idx], hi_90[sort_idx],
        alpha=0.25, color="steelblue", label="90% CI",
    )
    ax2.fill_between(
        x_plot,
        lo_80[sort_idx], hi_80[sort_idx],
        alpha=0.35, color="steelblue", label="80% CI",
    )
    ax2.plot(x_plot, stack_test[sort_idx], color="steelblue",
             lw=1.2, zorder=4, label="Stack prediction")
    ax2.scatter(
        x_plot, y_test_arr[sort_idx],
        s=6, color="tomato", alpha=0.55, zorder=5, label="True value",
    )
    ax2.set_xlabel("Test samples (sorted by prediction)", **label_kw)
    ax2.set_ylabel("Melting Point (K)", **label_kw)
    ax2.set_title("Conformal Prediction Intervals on Test Set", **title_kw)
    ax2.legend(fontsize=9, loc="upper left")
    ax2.grid(True, alpha=0.25)
    ax2.text(
        0.97, 0.03,
        f"Empirical 90% coverage = {inside_90:.3f}",
        transform=ax2.transAxes, fontsize=9,
        ha="right", va="bottom",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
    )

    # ── Panel 3: Coverage calibration curve ─────────────────────────────────
    ax3.plot([0, 1], [0, 1], "--", color="grey", lw=1.2, label="Perfect calibration")
    ax3.plot(nominal, empirical, "o-", color="steelblue", lw=2, ms=7,
             label="Conformal (stack)")
    ax3.fill_between(nominal, nominal, empirical, alpha=0.12, color="steelblue")
    ax3.set_xlim(0.45, 1.02)
    ax3.set_ylim(0.45, 1.02)
    ax3.set_xlabel("Nominal coverage level", **label_kw)
    ax3.set_ylabel("Empirical coverage on test set", **label_kw)
    ax3.set_title("Coverage Calibration Curve", **title_kw)
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.25)

    # ── Panel 4: Uncertainty distribution + per-base-model MAE inset ─────────
    model_labels = ["LightGBM", "XGBoost", "ExtraTrees", "DNN"]
    colors       = ["#4C72B0", "#DD8452", "#55A868", "#C44E52"]

    ax4.hist(base_std, bins=40, color="steelblue", alpha=0.65, density=True,
             label="Base-model std (test)")
    ax4.axvline(base_std.mean(), color="navy", lw=1.8, ls="--",
                label=f"Mean = {base_std.mean():.2f} K")
    ax4.axvline(np.median(base_std), color="darkorange", lw=1.8, ls=":",
                label=f"Median = {np.median(base_std):.2f} K")
    ax4.set_xlabel("Base-model std (K)", **label_kw)
    ax4.set_ylabel("Density", **label_kw)
    ax4.set_title("Distribution of Base-Model Disagreement", **title_kw)
    ax4.legend(fontsize=9, loc="lower right", framealpha=0.4)
    ax4.grid(True, alpha=0.25)

    ax_ins = ax4.inset_axes([0.58, 0.45, 0.40, 0.50])
    bm_maes = [mean_absolute_error(y_test_arr, test_matrix[:, i])
               for i in range(len(BASE_MODELS))]
    y_pos = np.arange(len(BASE_MODELS))
    bars  = ax_ins.barh(y_pos, bm_maes, color=colors, alpha=0.80)
    ax_ins.set_yticks(y_pos)
    ax_ins.set_yticklabels(model_labels)
    ax_ins.set_xlabel("MAE (K)", fontsize=8)
    ax_ins.set_title("Base MAE", fontsize=8)
    ax_ins.tick_params(labelsize=7)
    for bar, val in zip(bars, bm_maes):
        ax_ins.text(
            val + 0.5, bar.get_y() + bar.get_height() / 2,
            f"{val:.1f}", va="center", fontsize=7,
        )

    # ── Figure title ─────────────────────────────────────────────────────────
    fig.suptitle(
        f"Prediction Uncertainty Analysis — {dataset.upper()} / {fs_method} / {transform}  "
        f"[stack, no outlier removal]\n"
        f"Test MAE = {mae:.2f} K   R² = {r2:.4f}   n_test = {n_test}",
        fontsize=13, y=1.01,
    )

    out_path = os.path.join(RESULTS_DIR, f"uncertainty_analysis_{dataset}.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"\n  Saved → {out_path}")

    return {
        "dataset":         dataset,
        "fs_method":       fs_method,
        "transform":       transform,
        "mae":             mae,
        "r2":              r2,
        "conformal_q90":   q_90,
        "coverage_90":     inside_90,
        "spearman_rho":    rho,
        "spearman_p":      p_val,
        "mean_base_std":   base_std.mean(),
        "median_base_std": np.median(base_std),
    }


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    summary_rows = []
    for dataset, fs_method, transform in BEST_MODELS:
        row = analyse_dataset(dataset, fs_method, transform)
        summary_rows.append(row)

    print(f"\n{'='*70}")
    print("  SUMMARY ACROSS DATASETS")
    print(f"{'='*70}")
    hdr = (
        f"  {'Dataset':>6}  {'MAE':>7}  {'R²':>7}  "
        f"{'q(90%)':>8}  {'Cov@90%':>9}  "
        f"{'Spearman ρ':>11}  {'Mean Std':>9}"
    )
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for r in summary_rows:
        print(
            f"  {r['dataset']:>6}  {r['mae']:>7.2f}  {r['r2']:>7.4f}  "
            f"{r['conformal_q90']:>8.2f}  {r['coverage_90']:>9.3f}  "
            f"{r['spearman_rho']:>11.4f}  {r['mean_base_std']:>9.3f}"
        )


if __name__ == "__main__":
    main()
