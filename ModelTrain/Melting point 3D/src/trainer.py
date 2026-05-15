######################################################
# Developer : Methun Kamruzzaman
# Date      : April 9, 2026
# Purpose   : Runs a single training experiment for one pipeline configuration.
######################################################
"""
trainer.py
----------
Runs a single experiment for one combination of
    (dataset, transformation, feature_selection, model)

and returns the best model, its hyperparameters, and train/test metrics.

Pipeline per experiment
-----------------------
1. Apply scaler  : fitted on X_train, applied to X_train and X_test
2. Apply FS      : fitted on X_train_scaled + y_train,
                   applied to X_train_scaled and X_test_scaled
3. BayesSearchCV : CV on X_train_selected with neg_mean_absolute_error.
                   DNN uses a reduced budget (5-fold, 15 iter, 50 epochs)
                   to keep HPO fast; all other models use 10-fold, 20 iter.
4. Refit best estimator on full X_train_selected.
                   DNN is always refitted from scratch with full budget
                   (200 epochs, patience=20) after HPO completes.
5. Save model artefacts (scaler + FS + model) with joblib / torch.save
6. Return metrics dict

Saved artefact layout
---------------------
    saved_models/{dataset}/{fs}/{transform}/{model}/
        pipeline.joblib    – dict {scaler, selector, model, params, metrics}
"""

import os
import time
import warnings
import joblib
import numpy as np

from sklearn.preprocessing import StandardScaler, QuantileTransformer
from sklearn.model_selection import cross_val_score, KFold
from skopt import BayesSearchCV

from src.config import (
    MODELS_DIR,
    N_FOLDS,
    N_BAYES_ITER,
    N_BAYES_INIT,
    SCORING,
    RANDOM_STATE,
    N_JOBS,
    DNN_BAYES_N_JOBS,
    DNN_MAX_EPOCHS,
    DNN_PATIENCE,
    DNN_HPO_MAX_EPOCHS,
    DNN_HPO_PATIENCE,
    DNN_HPO_N_FOLDS,
    DNN_HPO_N_ITER,
    SEARCH_SPACES,
)
from src.models import get_model
from src.feature_selection import get_feature_selector
from src.evaluator import compute_metrics

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Scaler factory
# ---------------------------------------------------------------------------
def get_scaler(transform: str):
    if transform == "standard":
        return StandardScaler()
    if transform == "quantile":
        return QuantileTransformer(
            n_quantiles=1000,
            output_distribution="normal",
            random_state=RANDOM_STATE,
        )
    return None   # "none" → no transformation


# ---------------------------------------------------------------------------
# Load a saved artefact and re-evaluate on train + test
# ---------------------------------------------------------------------------
def load_and_evaluate(
    dataset:    str,
    transform:  str,
    fs_method:  str,
    model_name: str,
    X_train,
    X_test,
    y_train,
    y_test,
    verbose:    bool = True,
) -> dict:
    """
    Load a previously saved pipeline.joblib and evaluate it on the supplied
    data.  No HPO or retraining is performed.

    Returns the same result dict as run_experiment().
    """
    save_dir  = os.path.join(MODELS_DIR, dataset, fs_method, transform, model_name)
    artefact  = joblib.load(os.path.join(save_dir, "pipeline.joblib"))
    tag       = f"{dataset}/{fs_method}/{transform}/{model_name}"

    scaler    = artefact["scaler"]
    selector  = artefact["selector"]
    model     = artefact["model"]
    best_params = artefact.get("best_params", {})

    # For DNN the weights were saved separately; reload them
    if model_name == "dnn" and model is None:
        import torch
        from src.dnn_model import DNNRegressor
        model = DNNRegressor(**{k: v for k, v in best_params.items()
                                if k not in ("max_epochs", "patience")},
                             max_epochs=1, patience=1)  # dummy — weights overwritten
        # Re-build the network architecture so we can load state_dict
        import pandas as pd
        X_tmp = X_train.iloc[:2] if hasattr(X_train, "iloc") else X_train[:2]
        if scaler:
            X_tmp = scaler.transform(X_tmp)
            X_tmp = pd.DataFrame(X_tmp, columns=X_train.columns)
        X_tmp = selector.transform(X_tmp)
        model.fit(X_tmp, y_train.iloc[:2] if hasattr(y_train, "iloc") else y_train[:2])
        weights_path = artefact.get("dnn_weights_path",
                                    os.path.join(save_dir, "dnn_weights.pt"))
        model.model_.load_state_dict(
            torch.load(weights_path, map_location=model.device_)
        )

    # Apply the stored scaler
    import pandas as pd
    if scaler is not None:
        X_tr_s = scaler.transform(X_train)
        X_te_s = scaler.transform(X_test)
        X_tr_s = pd.DataFrame(X_tr_s, columns=X_train.columns, index=X_train.index)
        X_te_s = pd.DataFrame(X_te_s, columns=X_test.columns,  index=X_test.index)
    else:
        X_tr_s, X_te_s = X_train.copy(), X_test.copy()

    # Apply the stored selector
    X_tr_f = selector.transform(X_tr_s)
    X_te_f = selector.transform(X_te_s)

    n_features = X_tr_f.shape[1]

    y_train_pred  = model.predict(X_tr_f)
    y_test_pred   = model.predict(X_te_f)
    train_metrics = compute_metrics(y_train.values, y_train_pred)
    test_metrics  = compute_metrics(y_test.values,  y_test_pred)

    if verbose:
        print(f"  [{tag}] LOADED — "
              f"Test MAE={test_metrics['mae']:.3f}  "
              f"R²={test_metrics['r2']:.3f}  "
              f"RMSE={test_metrics['rmse']:.3f}")

    return {
        "dataset":           dataset,
        "transform":         transform,
        "fs_method":         fs_method,
        "model":             model_name,
        "n_features":        n_features,
        "selected_features": artefact.get("selected_features"),
        "best_params":       best_params,
        "train_mae":         train_metrics["mae"],
        "train_r2":          train_metrics["r2"],
        "train_rmse":        train_metrics["rmse"],
        "test_mae":          test_metrics["mae"],
        "test_r2":           test_metrics["r2"],
        "test_rmse":         test_metrics["rmse"],
        "elapsed_sec":       0.0,
    }


def _pipeline_path(dataset, fs_method, transform, model_name) -> str:
    return os.path.join(
        MODELS_DIR, dataset, fs_method, transform, model_name, "pipeline.joblib"
    )


# ---------------------------------------------------------------------------
# Main experiment function
# ---------------------------------------------------------------------------
def run_experiment(
    dataset:    str,
    transform:  str,
    fs_method:  str,
    model_name: str,
    X_train,
    X_test,
    y_train,
    y_test,
    verbose:    bool = True,
) -> dict:
    """
    Execute one full experiment and return a result dict.

    Returns
    -------
    dict with keys:
        best_params, train_mae, train_r2, train_rmse,
        test_mae,  test_r2,  test_rmse,
        n_features, elapsed_sec
    """
    t0 = time.time()
    tag = f"{dataset}/{fs_method}/{transform}/{model_name}"

    # ------------------------------------------------------------------
    # 1.  Scaling
    # ------------------------------------------------------------------
    scaler = get_scaler(transform)
    if scaler is not None:
        X_tr_s = scaler.fit_transform(X_train)
        X_te_s = scaler.transform(X_test)
        # Restore as DataFrames so FS methods can use column names
        import pandas as pd
        X_tr_s = pd.DataFrame(X_tr_s, columns=X_train.columns, index=X_train.index)
        X_te_s = pd.DataFrame(X_te_s, columns=X_test.columns,  index=X_test.index)
    else:
        X_tr_s, X_te_s = X_train.copy(), X_test.copy()

    # ------------------------------------------------------------------
    # 2.  Feature selection
    # ------------------------------------------------------------------
    selector = get_feature_selector(fs_method)
    X_tr_f   = selector.fit_transform(X_tr_s, y_train)
    X_te_f   = selector.transform(X_te_s)

    n_features = X_tr_f.shape[1] if hasattr(X_tr_f, "shape") else len(X_tr_f[0])

    # ------------------------------------------------------------------
    # 3.  Bayesian hyper-parameter optimisation
    # ------------------------------------------------------------------
    is_dnn       = model_name == "dnn"
    search_space = SEARCH_SPACES[model_name]
    n_jobs       = DNN_BAYES_N_JOBS if is_dnn else N_JOBS
    n_iter       = DNN_HPO_N_ITER   if is_dnn else N_BAYES_ITER
    n_cv_folds   = DNN_HPO_N_FOLDS  if is_dnn else N_FOLDS

    # For DNN: use a reduced epoch budget during search so each candidate
    # is evaluated quickly.  The winner is fully retrained below.
    if is_dnn:
        model = get_model("dnn")
        model.set_params(max_epochs=DNN_HPO_MAX_EPOCHS, patience=DNN_HPO_PATIENCE)
    else:
        model = get_model(model_name)

    cv = KFold(n_splits=n_cv_folds, shuffle=True, random_state=RANDOM_STATE)

    opt = BayesSearchCV(
        estimator        = model,
        search_spaces    = search_space,
        n_iter           = n_iter,
        optimizer_kwargs = {"n_initial_points": N_BAYES_INIT},
        scoring          = SCORING,
        cv               = cv,
        n_jobs           = n_jobs,
        refit            = not is_dnn,   # DNN refit handled manually below
        random_state     = RANDOM_STATE,
        verbose          = 0,
    )

    if verbose:
        budget_note = (f"HPO: {n_cv_folds}-fold × {n_iter} iter × "
                       f"{DNN_HPO_MAX_EPOCHS} epochs" if is_dnn else
                       f"{n_cv_folds}-fold × {n_iter} iter")
        print(f"  [{tag}] BayesSearchCV started  "
              f"(n_features={n_features}, {budget_note})")

    opt.fit(X_tr_f, y_train)
    best_params = dict(opt.best_params_)

    # DNN: refit from scratch with the full epoch budget on all training data
    if is_dnn:
        from src.dnn_model import DNNRegressor
        best_model = DNNRegressor(
            max_epochs = DNN_MAX_EPOCHS,
            patience   = DNN_PATIENCE,
            **{k: v for k, v in best_params.items()
               if k not in ("max_epochs", "patience")},
        )
        if verbose:
            print(f"  [{tag}] Refitting DNN with full budget "
                  f"({DNN_MAX_EPOCHS} epochs) ...")
        best_model.fit(X_tr_f, y_train)
    else:
        best_model = opt.best_estimator_

    if verbose:
        print(f"  [{tag}] Best CV MAE: {-opt.best_score_:.3f}  params: {best_params}")

    # ------------------------------------------------------------------
    # 4.  Evaluate on train + test
    # ------------------------------------------------------------------
    y_train_pred = best_model.predict(X_tr_f)
    y_test_pred  = best_model.predict(X_te_f)

    train_metrics = compute_metrics(y_train.values, y_train_pred)
    test_metrics  = compute_metrics(y_test.values,  y_test_pred)

    # ------------------------------------------------------------------
    # 5.  Save artefacts
    # ------------------------------------------------------------------
    save_dir = os.path.join(MODELS_DIR, dataset, fs_method, transform, model_name)
    os.makedirs(save_dir, exist_ok=True)

    # Collect selected feature names explicitly
    selected_features = selector.get_selected_features()

    artefact = {
        "scaler":            scaler,
        "selector":          selector,
        "model":             best_model,
        "best_params":       best_params,
        "train_metrics":     train_metrics,
        "test_metrics":      test_metrics,
        "n_features":        n_features,
        "selected_features": selected_features,   # explicit list of kept feature names
        "tag":               tag,
    }

    if model_name == "dnn":
        import torch
        model_path = os.path.join(save_dir, "dnn_weights.pt")
        torch.save(best_model.model_.state_dict(), model_path)
        artefact["model"] = None   # weights saved separately
        artefact["dnn_weights_path"] = model_path

    joblib.dump(artefact, os.path.join(save_dir, "pipeline.joblib"), compress=3)

    # Save feature names as a plain text file for easy inspection
    if selected_features is not None:
        feat_path = os.path.join(save_dir, "selected_features.txt")
        with open(feat_path, "w") as fh:
            fh.write("\n".join(selected_features))

    # ------------------------------------------------------------------
    # 6.  Return result dict
    # ------------------------------------------------------------------
    elapsed = time.time() - t0
    result = {
        "dataset":           dataset,
        "transform":         transform,
        "fs_method":         fs_method,
        "model":             model_name,
        "n_features":        n_features,
        "selected_features": selected_features,   # list of feature names (or PC labels for PCA)
        "best_params":       best_params,
        "train_mae":  train_metrics["mae"],
        "train_r2":   train_metrics["r2"],
        "train_rmse": train_metrics["rmse"],
        "test_mae":   test_metrics["mae"],
        "test_r2":    test_metrics["r2"],
        "test_rmse":  test_metrics["rmse"],
        "elapsed_sec": round(elapsed, 1),
    }

    if verbose:
        print(f"  [{tag}] Done — "
              f"Test MAE={test_metrics['mae']:.3f}  "
              f"R²={test_metrics['r2']:.3f}  "
              f"RMSE={test_metrics['rmse']:.3f}  "
              f"({elapsed:.0f}s)")

    return result
