######################################################
# Developer : Methun Kamruzzaman
# Date      : April 12, 2026
# Purpose   : Level-2 stacking ensemble (LightGBM, XGBoost, ExtraTrees, DNN)
#             for melting point (2D) prediction.
######################################################
"""
stacking.py
-----------
Level-2 stacking ensemble over LightGBM, XGBoost, Extra Trees, and DNN base models.

Workflow
--------
For each (dataset, transform, fs_method):
  1. Load saved base-model artifacts (scaler, selector, best_params, model).
  2. Generate k-fold OOF predictions on X_train using fresh model instances
     initialised with the saved best_params.  Preprocessing (scale + select)
     is applied once on the full training set using each model's own fitted
     scaler/selector from its artifact.
  3. Train a Ridge meta-learner on the [n_train × n_base] OOF feature matrix.
  4. Obtain test predictions from each fully-trained base model.
  5. Evaluate the stacked model: train MAE uses OOF preds, test MAE uses
     full-model preds fed through the meta-learner.
  6. Save the meta-learner artifact and return a result dict.
"""

import os
import time
import warnings
import joblib

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold

from src.config import (
    MODELS_DIR,
    RANDOM_STATE,
    DNN_HPO_MAX_EPOCHS,
    DNN_HPO_PATIENCE,
)
from src.evaluator import compute_metrics

warnings.filterwarnings("ignore")

BASE_MODELS   = ("lgbm", "xgb", "extratrees", "dnn")
META_ALPHA    = 1.0   # Ridge regularisation strength for the meta-learner
N_STACK_FOLDS = 5


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _artifact_path(dataset: str, fs_method: str, transform: str, model_name: str) -> str:
    return os.path.join(
        MODELS_DIR, dataset, fs_method, transform, model_name, "pipeline.joblib"
    )


def _artifacts_exist(dataset: str, fs_method: str, transform: str, base_models) -> bool:
    return all(
        os.path.exists(_artifact_path(dataset, fs_method, transform, mn))
        for mn in base_models
    )


def _preprocess(artifact: dict, X_train, X_test):
    """
    Apply the artifact's fitted scaler + selector to X_train and X_test.
    Returns (X_tr_arr, X_te_arr) as numpy arrays.
    """
    scaler   = artifact["scaler"]
    selector = artifact["selector"]

    if scaler is not None:
        X_tr_s = pd.DataFrame(
            scaler.transform(X_train),
            columns=X_train.columns,
            index=X_train.index,
        )
        X_te_s = pd.DataFrame(
            scaler.transform(X_test),
            columns=X_test.columns,
            index=X_test.index,
        )
    else:
        X_tr_s, X_te_s = X_train.copy(), X_test.copy()

    X_tr_f = selector.transform(X_tr_s)
    X_te_f = selector.transform(X_te_s)

    if hasattr(X_tr_f, "values"):
        X_tr_f = X_tr_f.values
    if hasattr(X_te_f, "values"):
        X_te_f = X_te_f.values

    return X_tr_f.astype(np.float32), X_te_f.astype(np.float32)


def _fresh_model(model_name: str, best_params: dict):
    """
    Return a new, unfitted model initialised with the saved best_params.
    DNN uses the reduced HPO epoch budget to keep OOF CV fast.
    """
    if model_name == "dnn":
        from src.dnn_model import DNNRegressor
        params = {k: v for k, v in best_params.items()
                  if k not in ("max_epochs", "patience")}
        return DNNRegressor(
            max_epochs=DNN_HPO_MAX_EPOCHS,
            patience=DNN_HPO_PATIENCE,
            **params,
        )
    from src.models import get_model
    m = get_model(model_name)
    m.set_params(**best_params)
    return m


def _reload_full_model(artifact: dict, model_name: str, X_tr_f: np.ndarray, y_train):
    """
    Return the fully-trained model stored in the artifact.
    For DNN the weights are reloaded from the separate .pt file.
    """
    model       = artifact["model"]
    best_params = artifact.get("best_params", {})

    if model_name == "dnn" and model is None:
        import torch
        from src.dnn_model import DNNRegressor

        params = {k: v for k, v in best_params.items()
                  if k not in ("max_epochs", "patience")}
        model = DNNRegressor(max_epochs=1, patience=1, **params)

        # Build the network architecture via a dummy 2-row fit
        y_arr = y_train.values if hasattr(y_train, "values") else np.asarray(y_train)
        model.fit(X_tr_f[:2], y_arr[:2])

        weights_path = artifact["dnn_weights_path"]
        model.model_.load_state_dict(
            torch.load(weights_path, map_location=model.device_)
        )

    return model


def _oof_predictions(
    model_name:  str,
    best_params: dict,
    X_tr_f:      np.ndarray,
    y_arr:       np.ndarray,
    n_folds:     int,
    verbose:     bool,
    tag:         str,
) -> np.ndarray:
    """K-fold OOF predictions using fresh model instances (no data leakage)."""
    kf  = KFold(n_splits=n_folds, shuffle=True, random_state=RANDOM_STATE)
    oof = np.zeros(len(y_arr))

    for fold, (tr_idx, val_idx) in enumerate(kf.split(X_tr_f), 1):
        if verbose:
            print(f"    {mn_label(model_name)} fold {fold}/{n_folds} ...")
        m = _fresh_model(model_name, best_params)
        m.fit(X_tr_f[tr_idx], y_arr[tr_idx])
        oof[val_idx] = m.predict(X_tr_f[val_idx])

    return oof


def mn_label(model_name: str) -> str:
    return {
        "lgbm":       "LightGBM",
        "xgb":        "XGBoost",
        "extratrees": "ExtraTrees",
        "dnn":        "DNN",
    }.get(model_name, model_name)


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------

def stacking_exists(dataset: str, transform: str, fs_method: str) -> bool:
    path = os.path.join(
        MODELS_DIR, dataset, fs_method, transform, "stack", "pipeline.joblib"
    )
    return os.path.exists(path)


def _oof_cache_path(save_dir: str, model_name: str) -> tuple[str, str]:
    """Return (oof_path, test_path) for cached predictions of a base model."""
    return (
        os.path.join(save_dir, f"oof_{model_name}.npy"),
        os.path.join(save_dir, f"test_{model_name}.npy"),
    )


def _load_oof_cache(save_dir: str, model_name: str) -> tuple[np.ndarray, np.ndarray] | None:
    """Return (oof, test_preds) from disk, or None if cache is missing."""
    oof_path, test_path = _oof_cache_path(save_dir, model_name)
    if os.path.exists(oof_path) and os.path.exists(test_path):
        return np.load(oof_path), np.load(test_path)
    return None


def _save_oof_cache(save_dir: str, model_name: str, oof: np.ndarray, test_preds: np.ndarray):
    oof_path, test_path = _oof_cache_path(save_dir, model_name)
    np.save(oof_path, oof)
    np.save(test_path, test_preds)


# ---------------------------------------------------------------------------
# Main experiment function
# ---------------------------------------------------------------------------

def run_stacking_experiment(
    dataset:     str,
    transform:   str,
    fs_method:   str,
    X_train,
    X_test,
    y_train,
    y_test,
    base_models: tuple = BASE_MODELS,
    n_folds:     int   = N_STACK_FOLDS,
    retrain_oof: bool  = False,
    verbose:     bool  = True,
) -> dict | None:
    """
    Run a stacking ensemble for one (dataset, transform, fs_method) combo.

    OOF predictions for each base model are cached to disk. On subsequent
    runs only models whose cache is missing are recomputed; the meta-learner
    is always retrained (Ridge, negligible cost).

    Parameters
    ----------
    retrain_oof : force recomputation of all OOF caches even if they exist.

    Returns a result dict (same schema as trainer.run_experiment) or None if
    any base-model artifact is missing.
    """
    t0  = time.time()
    tag = f"{dataset}/{fs_method}/{transform}/stack"

    # Guard: every base-model artifact must exist
    if not _artifacts_exist(dataset, fs_method, transform, base_models):
        missing = [
            mn for mn in base_models
            if not os.path.exists(_artifact_path(dataset, fs_method, transform, mn))
        ]
        if verbose:
            print(f"  [{tag}] Skipping — missing artifacts for: {missing}")
        return None

    artifacts = {
        mn: joblib.load(_artifact_path(dataset, fs_method, transform, mn))
        for mn in base_models
    }

    save_dir = os.path.join(MODELS_DIR, dataset, fs_method, transform, "stack")
    os.makedirs(save_dir, exist_ok=True)

    if verbose:
        print(
            f"\n  [{tag}] Stacking ensemble  "
            f"({', '.join(mn_label(m) for m in base_models)}, "
            f"{n_folds}-fold OOF)"
        )

    y_arr      = y_train.values if hasattr(y_train, "values") else np.asarray(y_train)
    y_test_arr = y_test.values  if hasattr(y_test,  "values") else np.asarray(y_test)

    oof_matrix  = np.zeros((len(y_arr),      len(base_models)))
    test_matrix = np.zeros((len(y_test_arr), len(base_models)))
    n_features  = None

    for col, mn in enumerate(base_models):
        art         = artifacts[mn]
        best_params = art.get("best_params", {})

        X_tr_f, X_te_f = _preprocess(art, X_train, X_test)
        if n_features is None:
            n_features = X_tr_f.shape[1]

        # Load OOF cache if available; otherwise compute and save
        cached = None if retrain_oof else _load_oof_cache(save_dir, mn)
        if cached is not None:
            oof, test_preds = cached
            if verbose:
                print(f"  [{tag}] OOF — {mn_label(mn)} ... CACHED")
        else:
            if verbose:
                print(f"  [{tag}] OOF — {mn_label(mn)} ...")
            kf  = KFold(n_splits=n_folds, shuffle=True, random_state=RANDOM_STATE)
            oof = np.zeros(len(y_arr))
            for fold, (tr_idx, val_idx) in enumerate(kf.split(X_tr_f), 1):
                if verbose:
                    print(f"    fold {fold}/{n_folds}")
                m = _fresh_model(mn, best_params)
                m.fit(X_tr_f[tr_idx], y_arr[tr_idx])
                oof[val_idx] = m.predict(X_tr_f[val_idx])

            full_model = _reload_full_model(art, mn, X_tr_f, y_train)
            test_preds = full_model.predict(X_te_f)
            _save_oof_cache(save_dir, mn, oof, test_preds)

        oof_matrix[:, col]  = oof
        test_matrix[:, col] = test_preds

    # -----------------------------------------------------------------------
    # Train Ridge meta-learner on OOF predictions
    # -----------------------------------------------------------------------
    meta = Ridge(alpha=META_ALPHA)
    meta.fit(oof_matrix, y_arr)

    weights = dict(zip(base_models, meta.coef_.tolist()))
    if verbose:
        fmt = {mn_label(k): round(v, 4) for k, v in weights.items()}
        print(f"  [{tag}] Meta-learner weights: {fmt}  intercept: {meta.intercept_:.4f}")

    # -----------------------------------------------------------------------
    # Evaluate
    # -----------------------------------------------------------------------
    y_train_pred  = meta.predict(oof_matrix)   # OOF — unbiased train estimate
    y_test_pred   = meta.predict(test_matrix)

    train_metrics = compute_metrics(y_arr,      y_train_pred)
    test_metrics  = compute_metrics(y_test_arr, y_test_pred)

    elapsed = time.time() - t0

    if verbose:
        print(
            f"  [{tag}] Done — "
            f"Test MAE={test_metrics['mae']:.3f}  "
            f"R²={test_metrics['r2']:.3f}  "
            f"RMSE={test_metrics['rmse']:.3f}  "
            f"({elapsed:.0f}s)"
        )

    # -----------------------------------------------------------------------
    # Save artifact
    # -----------------------------------------------------------------------
    joblib.dump(
        {
            "base_models":    list(base_models),
            "meta_learner":   meta,
            "n_folds":        n_folds,
            "meta_weights":   weights,
            "meta_intercept": float(meta.intercept_),
            "train_metrics":  train_metrics,
            "test_metrics":   test_metrics,
        },
        os.path.join(save_dir, "pipeline.joblib"),
        compress=3,
    )

    return {
        "dataset":       dataset,
        "transform":     transform,
        "fs_method":     fs_method,
        "model":         "stack",
        "n_features":    n_features,
        "meta_weights":  weights,
        "train_mae":     train_metrics["mae"],
        "train_r2":      train_metrics["r2"],
        "train_rmse":    train_metrics["rmse"],
        "test_mae":      test_metrics["mae"],
        "test_r2":       test_metrics["r2"],
        "test_rmse":     test_metrics["rmse"],
        "elapsed_sec":   round(elapsed, 1),
    }
