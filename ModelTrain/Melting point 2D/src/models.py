######################################################
# Developer : Methun Kamruzzaman
# Date      : April 3, 2026
# Purpose   : Factory returning unfitted sklearn-compatible estimators for the pipeline.
######################################################
"""
models.py
---------
Factory that returns an *unfitted* sklearn-compatible estimator for each
model name supported by the pipeline.

Model mapping
-------------
"linear"      →  Ridge regression (L2-regularised linear model)
"elasticnet"  →  ElasticNet (L1+L2, replaces SVR — fast, handles multicollinearity)
"extratrees"  →  Extra Trees Regressor (replaces RF — randomised splits, 2-4× faster)
"lgbm"        →  LightGBM (replaces GBR — histogram-based, leaf-wise, 10-20× faster)
"xgb"         →  XGBoost Regressor
"dnn"         →  Deep Neural Network (PyTorch, MPS-aware, sklearn wrapper)

All tree-based and sklearn models are configured with n_jobs=-1 where
applicable to exploit all CPU cores on M2 Max.
"""

from sklearn.linear_model import Ridge, ElasticNet
from sklearn.ensemble import ExtraTreesRegressor
from lightgbm import LGBMRegressor
from xgboost import XGBRegressor

from src.config import RANDOM_STATE, N_JOBS, DNN_MAX_EPOCHS, DNN_PATIENCE


def get_model(name: str):
    """Return a fresh, *unfitted* estimator for `name`."""
    builders = {
        "linear":     _build_linear,
        "elasticnet": _build_elasticnet,
        "extratrees": _build_extratrees,
        "lgbm":       _build_lgbm,
        "xgb":        _build_xgb,
        "dnn":        _build_dnn,
    }
    if name not in builders:
        raise ValueError(f"Unknown model: {name!r}. Choose from {list(builders)}")
    return builders[name]()


# ---------------------------------------------------------------------------
# Individual builders
# ---------------------------------------------------------------------------

def _build_linear():
    return Ridge(random_state=RANDOM_STATE)


def _build_elasticnet():
    # Replaces SVR: L1+L2 regularisation, handles high-dimensional descriptor
    # spaces with multicollinearity, trains in seconds instead of minutes
    return ElasticNet(
        max_iter=5000,
        random_state=RANDOM_STATE,
    )


def _build_extratrees():
    # Replaces RF: randomised split thresholds remove the expensive best-split
    # search, giving 2-4× speedup with comparable or better generalisation
    return ExtraTreesRegressor(
        n_jobs=N_JOBS,
        random_state=RANDOM_STATE,
    )


def _build_lgbm():
    # Replaces GBR: histogram binning + leaf-wise growth + native OpenMP
    # threading — typically 10-20× faster than sklearn GradientBoostingRegressor
    return LGBMRegressor(
        n_jobs=N_JOBS,
        random_state=RANDOM_STATE,
        verbose=-1,          # suppress LightGBM training output
    )


def _build_xgb():
    return XGBRegressor(
        tree_method="hist",
        device="cpu",        # Use all CPU cores; MPS reserved for DNN
        n_jobs=N_JOBS,
        random_state=RANDOM_STATE,
        verbosity=0,
        eval_metric="mae",
    )


def _build_dnn():
    from src.dnn_model import DNNRegressor
    return DNNRegressor(
        max_epochs=DNN_MAX_EPOCHS,
        patience=DNN_PATIENCE,
        random_state=RANDOM_STATE,
    )
