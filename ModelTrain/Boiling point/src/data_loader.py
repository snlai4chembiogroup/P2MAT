######################################################
# Developer : Methun Kamruzzaman
# Date      : April 3, 2026
# Purpose   : Loads and cleans train/test CSV pairs for the boiling point pipeline.
######################################################
"""
data_loader.py
--------------
Loads a train/test CSV pair, cleans the data, and returns (X_train, X_test, y_train, y_test)
as pandas DataFrames / Series of numeric features.

Cleaning steps
--------------
1. Drop metadata / target columns specified in config.DROP_COLS.
2. Drop any remaining non-numeric columns (e.g. 'scaffold').
3. Drop constant (zero-variance) columns – they add no information.
4. Replace ±inf with NaN then impute with column median (robustness guard).
5. Return feature names for downstream interpretability.
"""

import os
import numpy as np
import pandas as pd
from typing import Tuple, List

from src.config import DATA_DIR, DROP_COLS, TARGET_COL


def load_dataset(
    train_file: str,
    test_file: str,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series, pd.Series, List[str]]:
    """
    Load and clean a train/test dataset pair.

    Returns
    -------
    X_train, X_test : pd.DataFrame  – numeric feature matrices
    y_train, y_test : pd.Series     – target (MP_K)
    feature_names   : List[str]     – column names of X_train
    """
    train_path = os.path.join(DATA_DIR, train_file)
    test_path  = os.path.join(DATA_DIR, test_file)

    train_df = pd.read_csv(train_path)
    test_df  = pd.read_csv(test_path)

    y_train = train_df[TARGET_COL].copy()
    y_test  = test_df[TARGET_COL].copy()

    # Drop metadata + target from both splits
    cols_to_drop = [c for c in DROP_COLS + [TARGET_COL] if c in train_df.columns]
    train_df = train_df.drop(columns=cols_to_drop, errors="ignore")
    test_df  = test_df.drop(columns=cols_to_drop, errors="ignore")

    # Keep only columns that appear in both splits and are numeric
    common_cols = train_df.columns.intersection(test_df.columns).tolist()
    train_df = train_df[common_cols].select_dtypes(include="number")
    test_df  = test_df[train_df.columns]          # align to train column order

    # Replace ±inf with NaN, then median-impute
    train_df = _impute(train_df)
    test_df  = _impute(test_df)

    # Drop constant columns (fit on train only)
    non_const = train_df.columns[train_df.std() > 0].tolist()
    train_df  = train_df[non_const]
    test_df   = test_df[non_const]

    feature_names = train_df.columns.tolist()

    return train_df, test_df, y_train, y_test, feature_names


def _impute(df: pd.DataFrame) -> pd.DataFrame:
    """Replace ±inf → NaN, then fill NaN with column median."""
    df = df.replace([np.inf, -np.inf], np.nan)
    medians = df.median()
    return df.fillna(medians)
