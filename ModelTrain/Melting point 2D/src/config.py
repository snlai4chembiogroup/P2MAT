######################################################
# Developer : Methun Kamruzzaman
# Date      : April 8, 2026
# Purpose   : Central configuration for the melting point (2D) QSAR pipeline:
#             paths, experiment axes, and hyperparameter search spaces.
######################################################
"""
config.py
---------
Central configuration: paths, experiment axes, and hyperparameter search spaces.
All search-space objects use scikit-optimize (skopt) types.
"""

import os
from skopt.space import Real, Integer, Categorical

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR    = os.path.join(BASE_DIR, "Data")
MODELS_DIR  = os.path.join(BASE_DIR, "saved_models")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

for _d in (MODELS_DIR, RESULTS_DIR):
    os.makedirs(_d, exist_ok=True)

# ---------------------------------------------------------------------------
# Column definitions
# ---------------------------------------------------------------------------
DROP_COLS  = ["#PubChem", "key", "name", "smiles", "mpC", "scaffold"]
TARGET_COL = "MP_K"

# ---------------------------------------------------------------------------
# Experiment axes
# ---------------------------------------------------------------------------
DATASETS = {
    "full":  {"train": "train_set.csv",  "test": "test_set.csv"},
    "cho":   {"train": "train_cho.csv",  "test": "test_cho.csv"},
    "chon":  {"train": "train_chon.csv", "test": "test_chon.csv"},
}

TRANSFORMATIONS     = ["none", "standard", "quantile"]
FEATURE_SELECTIONS  = ["whole", "boruta", "pca", "correlation"]
MODEL_NAMES         = ["linear", "elasticnet", "extratrees", "lgbm", "xgb", "dnn"]

# ---------------------------------------------------------------------------
# Cross-validation & Bayesian optimisation
# ---------------------------------------------------------------------------
N_FOLDS       = 10
N_BAYES_ITER  = 20   # Bayesian optimisation iterations per model
N_BAYES_INIT  = 5    # Random initial points
SCORING       = "neg_mean_absolute_error"
RANDOM_STATE  = 42

# ---------------------------------------------------------------------------
# M2 Max hardware
# ---------------------------------------------------------------------------
N_JOBS = -1          # All CPU cores for sklearn models

# ---------------------------------------------------------------------------
# DNN training budget
# ---------------------------------------------------------------------------
DNN_MAX_EPOCHS    = 200   # epochs for final refit after HPO
DNN_PATIENCE      = 20    # early-stopping patience for final refit
DNN_BAYES_N_JOBS  = 1     # MPS doesn't support multi-process GPU → serial CV

# DNN HPO search budget (applied during BayesSearchCV only)
# Strategy: evaluate candidates cheaply, refit the winner with full budget.
#   200 epochs × 10-fold × 20 iter  →  50 epochs × 5-fold × 15 iter
#   = 40 000 epoch-runs  →  3 750 epoch-runs  (~10× speedup)
DNN_HPO_MAX_EPOCHS = 50   # max epochs per candidate during search
DNN_HPO_PATIENCE   = 10   # early-stopping patience during search
DNN_HPO_N_FOLDS    = 5    # CV folds for DNN HPO (others use N_FOLDS=10)
DNN_HPO_N_ITER     = 15   # Bayesian iterations for DNN (others use N_BAYES_ITER=20)

# ---------------------------------------------------------------------------
# Correlation-based feature selection thresholds
# ---------------------------------------------------------------------------
CORR_TARGET_THRESHOLD = 0.05   # |corr(feature, target)| must exceed this
CORR_INTER_THRESHOLD  = 0.85   # |corr(feat_i, feat_j)| → remove redundant

# ---------------------------------------------------------------------------
# Hyperparameter search spaces
# ---------------------------------------------------------------------------
SEARCH_SPACES = {
    # Ridge (L2-regularised linear model)
    "linear": {
        "alpha": Real(1e-4, 1e3, prior="log-uniform"),
    },

    # ElasticNet replaces SVR: L1+L2, handles multicollinearity, very fast
    "elasticnet": {
        "alpha":   Real(1e-4, 1e1, prior="log-uniform"),
        "l1_ratio": Real(0.0, 1.0),
    },

    # ExtraTreesRegressor replaces RF: randomised splits, no best-split search
    "extratrees": {
        "n_estimators":      Integer(100, 600),
        "max_depth":         Integer(3,   25),
        "min_samples_split": Integer(2,   20),
        "min_samples_leaf":  Integer(1,   10),
        "max_features":      Categorical(["sqrt", "log2", 0.3, 0.5]),
    },

    # LightGBM replaces GBR: histogram-based, leaf-wise, native multi-threading
    "lgbm": {
        "n_estimators":   Integer(100, 1000),
        "learning_rate":  Real(0.01, 0.3, prior="log-uniform"),
        "max_depth":      Integer(3,   12),
        "num_leaves":     Integer(15, 255),
        "subsample":      Real(0.5, 1.0),
        "colsample_bytree": Real(0.4, 1.0),
        "reg_alpha":      Real(1e-6, 10.0, prior="log-uniform"),
        "reg_lambda":     Real(1e-6, 10.0, prior="log-uniform"),
        "min_child_samples": Integer(5, 50),
    },

    "xgb": {
        "n_estimators":     Integer(100, 600),
        "learning_rate":    Real(0.01, 0.3, prior="log-uniform"),
        "max_depth":        Integer(3,  10),
        "subsample":        Real(0.5, 1.0),
        "colsample_bytree": Real(0.4, 1.0),
        "reg_alpha":        Real(1e-6, 10.0, prior="log-uniform"),
        "reg_lambda":       Real(1e-6, 10.0, prior="log-uniform"),
    },

    "dnn": {
        "n_layers":    Integer(2, 5),
        "hidden_size": Integer(64, 512),
        "dropout_rate": Real(0.0, 0.5),
        "learning_rate": Real(1e-4, 1e-2, prior="log-uniform"),
        "batch_size":  Integer(32, 256),
        "weight_decay": Real(1e-6, 1e-2, prior="log-uniform"),
        "activation":  Categorical(["relu", "tanh", "leaky_relu"]),
    },
}

# Human-readable labels for the results table
MODEL_DISPLAY_NAMES = {
    "linear":     "Linear Regression",
    "elasticnet": "ElasticNet",
    "extratrees": "Extra Trees Regressor",
    "lgbm":       "LightGBM",
    "xgb":        "XGBoost",
    "dnn":        "Deep Neural Network",
}

TRANSFORM_DISPLAY_NAMES = {
    "none":     "No Transformation",
    "standard": "Standard Scaler",
    "quantile": "Quantile Transformer",
}
