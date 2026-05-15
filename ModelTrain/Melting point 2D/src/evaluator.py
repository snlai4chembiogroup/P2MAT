######################################################
# Developer : Methun Kamruzzaman
# Date      : April 3, 2026
# Purpose   : Computes regression metrics and generates results tables for the pipeline.
######################################################
"""
evaluator.py
------------
Metric computation and results-table generation.

Metrics
-------
  MAE  – Mean Absolute Error
  R²   – Coefficient of Determination
  RMSE – Root Mean Squared Error

Results table
-------------
Rows    : (dataset, data_transformation)
Columns : model names
Cell    : best test-MAE across all feature-selection methods for that combo
"""

import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, r2_score, mean_squared_error

from src.config import (
    DATASETS,
    TRANSFORMATIONS,
    MODEL_NAMES,
    MODEL_DISPLAY_NAMES,
    TRANSFORM_DISPLAY_NAMES,
)


# ---------------------------------------------------------------------------
# Metric computation
# ---------------------------------------------------------------------------
def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    """Return dict with mae, r2, rmse."""
    mae  = mean_absolute_error(y_true, y_pred)
    r2   = r2_score(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    return {"mae": round(mae, 4), "r2": round(r2, 4), "rmse": round(rmse, 4)}


# ---------------------------------------------------------------------------
# Results table
# ---------------------------------------------------------------------------
def build_results_table(results: list[dict], phase: str = "test") -> pd.DataFrame:
    """
    Build a summary table from the list of experiment result dicts.

    Parameters
    ----------
    results : list of dicts returned by trainer.run_experiment
    phase   : 'train' or 'test' – which MAE to display

    Returns
    -------
    pd.DataFrame with MultiIndex rows (Dataset, Data Transformation)
    and model columns showing best MAE for that combo.
    """
    mae_key = f"{phase}_mae"

    records = []
    for r in results:
        records.append({
            "Dataset":           r["dataset"].upper(),
            "Data Transformation": TRANSFORM_DISPLAY_NAMES.get(
                                       r["transform"], r["transform"]),
            "Model":             MODEL_DISPLAY_NAMES.get(r["model"], r["model"]),
            "MAE":               r[mae_key],
        })

    df = pd.DataFrame(records)

    # Pivot: rows = (Dataset, Transform), columns = Model, values = min MAE
    pivot = df.groupby(
        ["Dataset", "Data Transformation", "Model"]
    )["MAE"].min().reset_index()

    table = pivot.pivot_table(
        index   = ["Dataset", "Data Transformation"],
        columns = "Model",
        values  = "MAE",
        aggfunc = "min",
    )

    # Order columns by MODEL_DISPLAY_NAMES order
    ordered_cols = [MODEL_DISPLAY_NAMES[m] for m in MODEL_NAMES
                    if MODEL_DISPLAY_NAMES[m] in table.columns]
    table = table[ordered_cols]
    table.columns.name = None

    return table


def build_full_report(results: list[dict]) -> pd.DataFrame:
    """
    Full experiment report: one row per (dataset, transform, fs, model)
    with all metrics for train and test.
    """
    rows = []
    for r in results:
        rows.append({
            "Dataset":      r["dataset"].upper(),
            "Transform":    TRANSFORM_DISPLAY_NAMES.get(r["transform"], r["transform"]),
            "Feature Selection": r["fs_method"].capitalize(),
            "Model":        MODEL_DISPLAY_NAMES.get(r["model"], r["model"]),
            "N Features":   r["n_features"],
            "Train MAE":    r["train_mae"],
            "Train R²":     r["train_r2"],
            "Train RMSE":   r["train_rmse"],
            "Test MAE":     r["test_mae"],
            "Test R²":      r["test_r2"],
            "Test RMSE":    r["test_rmse"],
            "Time (s)":     r["elapsed_sec"],
        })
    return pd.DataFrame(rows)


def print_results_table(table: pd.DataFrame, title: str = "Best Test MAE"):
    """Pretty-print the pivot table."""
    divider = "=" * 120
    print(f"\n{divider}")
    print(f"  {title}")
    print(divider)
    print(table.to_string(float_format=lambda x: f"{x:.3f}"))
    print(divider + "\n")
