######################################################
# Developer : Methun Kamruzzaman
# Date      : May 9, 2026
# Purpose   : Applies the best full-dataset stacking model to detect and remove
#             outliers from the complete compound set.
######################################################
"""
outlier_detection_full.py
--------------------------
Apply the saved best full-dataset stacking model (full / whole / standard)
to all compounds in MP_all_removed.csv, detect outliers via the Williams-plot
approach (leverage + standardised residuals), annotate each compound, and save
the result as Data/MP_all_outlier_removed.csv.

Outlier criteria (same as outlier_detection.py)
------------------------------------------------
  Response outlier  : |std_residual| > 3
  Leverage outlier  : h > h*  where h* = 3*(k+1)/n_train
  Bad leverage point: both criteria

Usage
-----
    python outlier_detection_full.py
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import joblib

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_DIR  = os.path.dirname(os.path.abspath(__file__))
DATA_DIR     = os.path.join(PROJECT_DIR, "Data")
MODELS_DIR   = os.path.join(PROJECT_DIR, "saved_models")
MP_ALL_PATH  = os.path.join(os.path.dirname(PROJECT_DIR), "Data", "MP_all_removed.csv")
OUTPUT_PATH  = os.path.join(os.path.dirname(PROJECT_DIR), "Data", "MP_all_outlier_removed.csv")

sys.path.insert(0, PROJECT_DIR)

# ---------------------------------------------------------------------------
# Configuration — best model for the full dataset (whole FS, standard scaler)
# ---------------------------------------------------------------------------
DATASET     = "full"
FS_METHOD   = "whole"
TRANSFORM   = "standard"
BASE_MODELS = ["lgbm", "xgb", "extratrees", "dnn"]
DROP_COLS   = ["#PubChem", "key", "name", "smiles", "mpC", "scaffold"]
TARGET_COL  = "MP_K"
STD_RES_THRESH = 3.0

# ---------------------------------------------------------------------------
# Load full dataset
# ---------------------------------------------------------------------------
print(f"Loading {MP_ALL_PATH} ...")
df_all = pd.read_csv(MP_ALL_PATH)
y_all  = df_all[TARGET_COL].values
print(f"  {len(df_all)} compounds, target column = {TARGET_COL}")

# ---------------------------------------------------------------------------
# Generate base-model predictions for every compound
# ---------------------------------------------------------------------------
meta_features = np.zeros((len(df_all), len(BASE_MODELS)), dtype=np.float64)

for col, mn in enumerate(BASE_MODELS):
    art_path = os.path.join(MODELS_DIR, DATASET, FS_METHOD, TRANSFORM, mn, "pipeline.joblib")
    print(f"\n[{mn}] Loading artifact: {art_path}")
    art         = joblib.load(art_path)
    scaler      = art["scaler"]
    selector    = art["selector"]
    model       = art["model"]
    best_params = art.get("best_params", {})

    # --- Align feature columns to what the scaler expects ---
    if scaler is not None and hasattr(scaler, "feature_names_in_"):
        feat_cols = list(scaler.feature_names_in_)
    else:
        feat_cols = [c for c in df_all.columns if c not in DROP_COLS + [TARGET_COL]
                     and pd.api.types.is_numeric_dtype(df_all[c])]

    X = df_all[feat_cols].copy()
    X = X.replace([np.inf, -np.inf], np.nan)
    X = X.fillna(X.median())

    # --- Scale ---
    if scaler is not None:
        X_scaled = pd.DataFrame(scaler.transform(X), columns=feat_cols, index=df_all.index)
    else:
        X_scaled = X

    # --- Select features (whole = pass-through) ---
    X_sel = selector.transform(X_scaled)
    if hasattr(X_sel, "values"):
        X_sel = X_sel.values
    X_sel = X_sel.astype(np.float32)
    print(f"  Feature shape after selection: {X_sel.shape}")

    # --- Predict ---
    if mn == "dnn" and model is None:
        import torch
        from src.dnn_model import DNNRegressor

        weights_path = os.path.join(MODELS_DIR, DATASET, FS_METHOD, TRANSFORM, mn, "dnn_weights.pt")
        params = {k: v for k, v in best_params.items() if k not in ("max_epochs", "patience")}
        model = DNNRegressor(max_epochs=1, patience=1, **params)
        model.fit(X_sel[:2], y_all[:2])          # build network architecture
        model.model_.load_state_dict(
            torch.load(weights_path, map_location=model.device_, weights_only=True)
        )
        print(f"  DNN weights loaded from {weights_path}")

    meta_features[:, col] = model.predict(X_sel)
    print(f"  [{mn}] predictions done")

# ---------------------------------------------------------------------------
# Stack through meta-learner
# ---------------------------------------------------------------------------
stack_dir  = os.path.join(MODELS_DIR, DATASET, FS_METHOD, TRANSFORM, "stack")
stack_art  = joblib.load(os.path.join(stack_dir, "pipeline.joblib"))
meta       = stack_art["meta_learner"]
y_pred_all = meta.predict(meta_features)
print(f"\nStacking done. Predicted {len(y_pred_all)} compounds.")

# ---------------------------------------------------------------------------
# Load OOF matrix (training reference for leverage & residual std)
# ---------------------------------------------------------------------------
oof_matrix = np.stack(
    [np.load(os.path.join(stack_dir, f"oof_{mn}.npy")) for mn in BASE_MODELS],
    axis=1,
)   # shape: (n_train, 4)

train_df = pd.read_csv(os.path.join(DATA_DIR, "train_set.csv"))
y_train  = train_df[TARGET_COL].values

y_oof     = meta.predict(oof_matrix)     # OOF-based train predictions
res_train = y_train - y_oof
std_train = np.std(res_train, ddof=1)    # reference std for standardisation

n_train = oof_matrix.shape[0]
k       = oof_matrix.shape[1]            # 4 base models
h_star  = 3 * (k + 1) / n_train

print(f"\nTraining reference: n_train={n_train}, k={k}, h*={h_star:.6f}")
print(f"Training residual std (OOF): {std_train:.4f} K")

# ---------------------------------------------------------------------------
# Standardised residuals & leverage for all compounds
# ---------------------------------------------------------------------------
res_all     = y_all - y_pred_all
std_res_all = res_all / std_train

XtX_inv = np.linalg.pinv(oof_matrix.T @ oof_matrix)
h_all   = np.einsum("ij,jk,ik->i", meta_features, XtX_inv, meta_features)

# ---------------------------------------------------------------------------
# Classify outliers
# ---------------------------------------------------------------------------
high_res = np.abs(std_res_all) > STD_RES_THRESH
high_lev = h_all > h_star
bad_lev  = high_res & high_lev
any_flag = high_res | high_lev

print(f"\n{'='*60}")
print(f"  Outlier Detection — {DATASET.upper()} / {FS_METHOD} / {TRANSFORM}")
print(f"  Total compounds          : {len(df_all)}")
print(f"  Response outliers        : {high_res.sum()}  (|std_res| > {STD_RES_THRESH})")
print(f"  Leverage outliers        : {high_lev.sum()}  (h > h*={h_star:.6f})")
print(f"  Bad leverage points      : {bad_lev.sum()}  (both)")
print(f"  Total flagged            : {any_flag.sum()}")
print(f"{'='*60}")

# ---------------------------------------------------------------------------
# Annotate dataframe
# ---------------------------------------------------------------------------
df_all["y_pred_K"]              = np.round(y_pred_all, 3)
df_all["residual_K"]            = np.round(res_all, 3)
df_all["std_residual"]          = np.round(std_res_all, 4)
df_all["leverage_h"]            = np.round(h_all, 6)
df_all["h_star"]                = round(h_star, 6)
df_all["is_response_outlier"]   = high_res
df_all["is_leverage_outlier"]   = high_lev
df_all["is_bad_leverage_point"] = bad_lev
df_all["is_outlier"]            = any_flag


def _outlier_type(row):
    t = []
    if row["is_response_outlier"]:  t.append("response_outlier")
    if row["is_leverage_outlier"]:  t.append("leverage_outlier")
    return "+".join(t) if t else "none"


df_all["outlier_type"] = df_all.apply(_outlier_type, axis=1)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
df_all.to_csv(OUTPUT_PATH, index=False)
print(f"\nSaved → {OUTPUT_PATH}")
print(f"  Rows: {len(df_all)}  |  Columns: {df_all.shape[1]}")
print(f"  New columns: y_pred_K, residual_K, std_residual, leverage_h, h_star,")
print(f"               is_response_outlier, is_leverage_outlier,")
print(f"               is_bad_leverage_point, is_outlier, outlier_type")
