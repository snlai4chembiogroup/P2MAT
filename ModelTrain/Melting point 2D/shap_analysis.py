######################################################
# Developer : Methun Kamruzzaman
# Date      : May 2, 2026
# Purpose   : SHAP feature-importance analysis for the best melting point (2D)
#             stacking ensemble models.
######################################################
"""
shap_analysis.py
----------------
SHAP feature-importance analysis for the three best stacking ensembles.

Best configurations
-------------------
  cho  / whole / standard
  chon / whole / quantile
  full / whole / standard

Method
------
The Ridge meta-learner is linear:
    y_stack = Σ_i  w_i · y_base_i  + intercept

Because the meta-learner is linear, ensemble SHAP values decompose as:
    SHAP_stack(x) = Σ_i  w_i · SHAP_base_i(x)

  lgbm, xgb             → shap.TreeExplainer   (exact)
  extratrees            → shap.TreeExplainer   (approximate=True, much faster)
  dnn                   → shap.DeepExplainer    (expected gradients,
                                                  background = 200 train samples)

  The test set is subsampled to MAX_EXPLAIN_SAMPLES (default 500) before SHAP
  computation.  This is standard practice for large datasets and reduces
  runtime from O(hours) to O(minutes) with negligible effect on the feature-
  importance ranking.

Outputs  (saved to results/shap/)
------
  {dataset}_shap_beeswarm.png  – top-20 beeswarm (test subsample)
  {dataset}_shap_bar.png       – top-20 mean |SHAP| bar chart
  {dataset}_shap_values.npz    – raw SHAP array + feature names + sample indices
"""

import os
import sys
import warnings

import joblib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shap
import torch

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.config import DATASETS, MODELS_DIR, RESULTS_DIR, RANDOM_STATE
from src.data_loader import load_dataset
from src.stacking import BASE_MODELS

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BEST_MODELS = [
    ("cho",  "whole", "standard"),
    ("chon", "whole", "quantile"),
    ("full", "whole", "standard"),
]

SHAP_DIR           = os.path.join(RESULTS_DIR, "shap")
TOP_N              = 20
DNN_BG_SIZE        = 200
MAX_EXPLAIN_SAMPLES = 500   # subsample test set; keeps runtime to ~minutes


# ---------------------------------------------------------------------------
# Artifact helpers
# ---------------------------------------------------------------------------

def _artifact_path(dataset, fs_method, transform, model_name):
    return os.path.join(MODELS_DIR, dataset, fs_method, transform, model_name, "pipeline.joblib")


def _preprocess(artifact, X: pd.DataFrame) -> np.ndarray:
    """Scale + select features using a base-model artifact. Returns float32 array."""
    scaler   = artifact["scaler"]
    selector = artifact["selector"]

    if scaler is not None:
        X_s = pd.DataFrame(scaler.transform(X), columns=X.columns, index=X.index)
    else:
        X_s = X.copy()

    X_f = selector.transform(X_s)
    if hasattr(X_f, "values"):
        X_f = X_f.values
    return X_f.astype(np.float32)


def _load_dnn(artifact, dataset, fs_method, transform, X_tr_f, y_train):
    """
    Reconstruct and reload the DNN from disk.
    The stored dnn_weights_path may point to a stale location; fall back to
    the canonical path under MODELS_DIR.
    """
    from src.dnn_model import DNNRegressor

    best_params = artifact.get("best_params", {})
    params = {k: v for k, v in best_params.items()
              if k not in ("max_epochs", "patience")}

    model = DNNRegressor(max_epochs=1, patience=1, **params)

    y_arr = y_train.values if hasattr(y_train, "values") else np.asarray(y_train)
    model.fit(X_tr_f[:2], y_arr[:2])  # build architecture with dummy fit

    # Prefer canonical path; fall back to the path stored in the artifact
    canonical = os.path.join(MODELS_DIR, dataset, fs_method, transform,
                             "dnn", "dnn_weights.pt")
    weights_path = canonical if os.path.exists(canonical) else artifact.get("dnn_weights_path")

    if weights_path is None or not os.path.exists(weights_path):
        raise FileNotFoundError(f"DNN weights not found for {dataset}/{fs_method}/{transform}")

    model.model_.load_state_dict(torch.load(weights_path, map_location=model.device_))
    return model


# ---------------------------------------------------------------------------
# SHAP computation per base model
# ---------------------------------------------------------------------------

def _shap_tree(model, X_arr: np.ndarray, approximate: bool = False) -> np.ndarray:
    """TreeSHAP — exact for lgbm/xgb; approximate mode for extratrees (faster)."""
    explainer = shap.TreeExplainer(model)
    sv = explainer.shap_values(X_arr, check_additivity=False,
                               approximate=approximate)
    if isinstance(sv, list):
        sv = sv[0]
    return sv  # [n_samples, n_features]


def _shap_dnn(dnn_model, X_train_arr: np.ndarray, X_explain_arr: np.ndarray) -> np.ndarray:
    """DeepExplainer SHAP for PyTorch DNN.

    The model is moved to CPU and set to eval mode before explanation so that
    BatchNorm uses running statistics (not batch statistics).

    DeepExplainer (SHAP 0.50) requires 2D output [batch, n_outputs], but
    _QSARNet.forward() calls .squeeze(-1) giving 1D [batch].  The thin
    _OutputWrapper below restores the missing dimension.
    """
    dnn_model.model_.cpu()
    dnn_model.model_.eval()

    class _OutputWrapper(torch.nn.Module):
        """Re-adds the output dimension removed by .squeeze(-1)."""
        def __init__(self, m: torch.nn.Module):
            super().__init__()
            self.m = m
        def forward(self, x: torch.Tensor) -> torch.Tensor:
            return self.m(x).unsqueeze(-1)   # [batch] → [batch, 1]

    rng    = np.random.RandomState(RANDOM_STATE)
    bg_idx = rng.choice(len(X_train_arr), size=min(DNN_BG_SIZE, len(X_train_arr)), replace=False)
    background = torch.tensor(X_train_arr[bg_idx], dtype=torch.float32)

    wrapper   = _OutputWrapper(dnn_model.model_)
    explainer = shap.DeepExplainer(wrapper, background)
    X_tensor  = torch.tensor(X_explain_arr, dtype=torch.float32)

    sv = explainer.shap_values(X_tensor)
    # DeepExplainer returns a list (one array per output); we have one output
    if isinstance(sv, list):
        sv = sv[0]
    sv = np.array(sv)
    if sv.ndim == 3:          # [n_samples, n_features, 1] → [n_samples, n_features]
        sv = sv[:, :, 0]
    return sv


# ---------------------------------------------------------------------------
# Ensemble SHAP
# ---------------------------------------------------------------------------

def compute_ensemble_shap(
    dataset: str,
    fs_method: str,
    transform: str,
    X_train: pd.DataFrame,
    X_test: pd.DataFrame,
    y_train,
) -> tuple[np.ndarray, list[str], np.ndarray]:
    """
    Compute weighted-sum SHAP values for the stacking ensemble.

    The test set is subsampled to MAX_EXPLAIN_SAMPLES before computation.

    Returns
    -------
    ensemble_shap : np.ndarray  shape [n_explain, n_features]
    feature_names : list[str]
    explain_idx   : np.ndarray  indices of the selected test rows
    """
    stack_path   = os.path.join(MODELS_DIR, dataset, fs_method, transform,
                                "stack", "pipeline.joblib")
    stack_art    = joblib.load(stack_path)
    meta_weights = stack_art["meta_weights"]

    fmt = {k: round(v, 4) for k, v in meta_weights.items()}
    print(f"  Meta-learner weights: {fmt}")

    # Subsample test set
    rng = np.random.RandomState(RANDOM_STATE)
    n_explain   = min(MAX_EXPLAIN_SAMPLES, len(X_test))
    explain_idx = rng.choice(len(X_test), size=n_explain, replace=False)
    explain_idx = np.sort(explain_idx)
    X_explain   = X_test.iloc[explain_idx]
    print(f"  Explaining {n_explain} / {len(X_test)} test samples")

    feature_names  = X_train.columns.tolist()
    ensemble_shap  = np.zeros((n_explain, len(feature_names)), dtype=np.float64)

    for mn in BASE_MODELS:
        weight = meta_weights.get(mn, 0.0)
        print(f"  [{mn}]  weight={weight:.4f}  computing SHAP ...", flush=True)

        art         = joblib.load(_artifact_path(dataset, fs_method, transform, mn))
        X_tr_arr    = _preprocess(art, X_train)
        X_expl_arr  = _preprocess(art, X_explain)

        if mn == "dnn":
            dnn = _load_dnn(art, dataset, fs_method, transform, X_tr_arr, y_train)
            sv  = _shap_dnn(dnn, X_tr_arr, X_expl_arr)
        else:
            # approximate=True for ExtraTrees: mean-path estimator, ~10× faster
            sv = _shap_tree(art["model"], X_expl_arr,
                            approximate=(mn == "extratrees"))

        ensemble_shap += weight * sv
        print(f"  [{mn}]  done — SHAP shape: {sv.shape}")

    return ensemble_shap, feature_names, explain_idx


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _save_beeswarm(shap_values, X_display, feature_names, dataset, out_dir):
    plt.figure(figsize=(10, 8))
    shap.summary_plot(
        shap_values,
        pd.DataFrame(X_display, columns=feature_names),
        max_display=TOP_N,
        show=False,
        plot_size=None,
    )
    plt.title(
        f"SHAP Beeswarm — {dataset.upper()} stacking ensemble  (test set, top {TOP_N})",
        fontsize=12,
    )
    plt.tight_layout()
    path = os.path.join(out_dir, f"{dataset}_shap_beeswarm.png")
    plt.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved → {path}")


def _save_bar(shap_values, feature_names, dataset, out_dir):
    plt.figure(figsize=(9, 7))
    shap.summary_plot(
        shap_values,
        feature_names=feature_names,
        plot_type="bar",
        max_display=TOP_N,
        show=False,
        plot_size=None,
    )
    plt.title(
        f"SHAP Feature Importance — {dataset.upper()} stacking ensemble  (top {TOP_N})",
        fontsize=12,
    )
    plt.tight_layout()
    path = os.path.join(out_dir, f"{dataset}_shap_bar.png")
    plt.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved → {path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(SHAP_DIR, exist_ok=True)

    for dataset, fs_method, transform in BEST_MODELS:
        print(f"\n{'='*64}")
        print(f"  SHAP  {dataset.upper()} / {fs_method} / {transform}")
        print(f"{'='*64}")

        ds_cfg = DATASETS[dataset]
        X_train, X_test, y_train, y_test, _ = load_dataset(
            ds_cfg["train"], ds_cfg["test"]
        )
        print(f"  Train: {X_train.shape}   Test: {X_test.shape}")

        shap_values, feature_names, explain_idx = compute_ensemble_shap(
            dataset, fs_method, transform, X_train, X_test, y_train
        )

        # Save raw arrays
        npz_path = os.path.join(SHAP_DIR, f"{dataset}_shap_values.npz")
        np.savez(npz_path,
                 shap_values=shap_values,
                 feature_names=np.array(feature_names),
                 explain_idx=explain_idx)
        print(f"  Saved raw SHAP → {npz_path}")

        # For display: original (unscaled) values of the explained test subset
        X_te_orig = X_test.iloc[explain_idx].values

        _save_beeswarm(shap_values, X_te_orig, feature_names, dataset, SHAP_DIR)
        _save_bar(shap_values, feature_names, dataset, SHAP_DIR)


if __name__ == "__main__":
    main()
