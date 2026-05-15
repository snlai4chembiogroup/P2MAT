######################################################
# Developer : Methun Kamruzzaman
# Date      : April 4, 2026
# Purpose   : Feature selection strategies (Boruta, variance, correlation, mutual info)
#             with sklearn-style interfaces for the QSAR pipeline.
######################################################
"""
feature_selection.py
--------------------
Four feature selection strategies, each implemented as a class with a
sklearn-style fit / transform / fit_transform interface.

Strategies
----------
WholeSelector      – no selection (pass-through)
BorutaSelector     – Boruta with XGBoost estimator
PCASelector        – top k PCs covering 90 % of explained variance
CorrelationSelector – Pearson correlation: high corr with target,
                      low inter-feature correlation (VIF-light)

All selectors are **fitted on training data only** then applied to test data.
"""

import sys
import os
import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.decomposition import PCA

# Boruta lives at the project root
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from boruta import BorutaPy  # noqa: E402

from src.config import CORR_TARGET_THRESHOLD, CORR_INTER_THRESHOLD, RANDOM_STATE


# ---------------------------------------------------------------------------
# Whole (pass-through)
# ---------------------------------------------------------------------------
class WholeSelector(BaseEstimator, TransformerMixin):
    """Returns all features unchanged."""

    def fit(self, X, y=None):
        if isinstance(X, pd.DataFrame):
            self.feature_names_in_ = X.columns.tolist()
        self.n_features_in_ = X.shape[1]
        return self

    def transform(self, X):
        return X

    def get_selected_features(self):
        return getattr(self, "feature_names_in_", None)


# ---------------------------------------------------------------------------
# Boruta
# ---------------------------------------------------------------------------
class BorutaSelector(BaseEstimator, TransformerMixin):
    """
    Boruta feature selection using XGBoost as the base estimator.

    Parameters
    ----------
    n_estimators : int or 'auto'
    perc         : int   – percentile of shadow-feature importance as threshold
    max_iter     : int   – max Boruta iterations
    random_state : int
    """

    def __init__(self, n_estimators="auto", perc=100, max_iter=100,
                 random_state=RANDOM_STATE):
        self.n_estimators = n_estimators
        self.perc         = perc
        self.max_iter     = max_iter
        self.random_state = random_state

    def fit(self, X, y):
        from xgboost import XGBRegressor

        estimator = XGBRegressor(
            n_estimators=100,
            max_depth=5,
            tree_method="hist",
            device="cpu",        # Boruta runs many fits – keep on CPU
            n_jobs=-1,
            random_state=self.random_state,
            verbosity=0,
        )
        selector = BorutaPy(
            estimator=estimator,
            n_estimators=self.n_estimators,
            perc=self.perc,
            max_iter=self.max_iter,
            random_state=self.random_state,
            verbose=0,
            early_stopping=True,
            n_iter_no_change=20,
        )

        X_arr = X.values if isinstance(X, pd.DataFrame) else X
        y_arr = y.values if isinstance(y, pd.Series) else y

        selector.fit(X_arr, y_arr)
        self._support = selector.support_ | selector.support_weak_

        if isinstance(X, pd.DataFrame):
            self.feature_names_in_ = X.columns.tolist()
            self.selected_features_ = [
                f for f, s in zip(self.feature_names_in_, self._support) if s
            ]
        return self

    def transform(self, X):
        if isinstance(X, pd.DataFrame):
            return X.iloc[:, self._support]
        return X[:, self._support]

    def get_selected_features(self):
        return getattr(self, "selected_features_", None)


# ---------------------------------------------------------------------------
# PCA
# ---------------------------------------------------------------------------
class PCASelector(BaseEstimator, TransformerMixin):
    """
    Retains principal components covering at least `variance_threshold` of
    total variance, with a minimum of `min_components` PCs enforced.

    Logic
    -----
    1. Fit full PCA and compute cumulative explained variance.
    2. Find k = smallest number of PCs that reach `variance_threshold`.
    3. If k >= min_components  → use k as-is.
       If k <  min_components  → force k = min_components (relaxing the
       variance constraint) and report the variance those extra PCs explain.

    Attributes (set after fit)
    --------------------------
    n_components_     : int   – final number of PCs selected
    actual_variance_  : float – cumulative variance explained by those PCs
    variance_relaxed_ : bool  – True if min_components forced k above threshold

    Input features **must already be scaled** before PCA.
    Returns a dense numpy array (PCA components are anonymous).
    """

    def __init__(self, variance_threshold: float = 0.90,
                 min_components: int = 10,
                 random_state: int = RANDOM_STATE):
        self.variance_threshold = variance_threshold
        self.min_components     = min_components
        self.random_state       = random_state

    def fit(self, X, y=None):
        X_arr = X.values if isinstance(X, pd.DataFrame) else X

        # Full PCA to get explained variance for all components
        pca_full = PCA(random_state=self.random_state)
        pca_full.fit(X_arr)
        cumvar = np.cumsum(pca_full.explained_variance_ratio_)

        # k from variance threshold (1-based → searchsorted + 1)
        k_var = int(np.searchsorted(cumvar, self.variance_threshold)) + 1
        k_var = max(k_var, 1)

        # Enforce minimum number of components
        k = max(k_var, self.min_components)
        # Never exceed available components
        k = min(k, X_arr.shape[1])

        self.n_components_     = k
        self.actual_variance_  = float(cumvar[k - 1])
        self.variance_relaxed_ = k > k_var

        if self.variance_relaxed_:
            print(f"  [PCASelector] 90% variance → {k_var} PCs (< {self.min_components} min). "
                  f"Relaxed to {k} PCs — explains {self.actual_variance_:.1%} variance.")
        else:
            print(f"  [PCASelector] {k} PCs selected — explains {self.actual_variance_:.1%} variance.")

        self._pca = PCA(n_components=k, random_state=self.random_state)
        self._pca.fit(X_arr)
        return self

    def transform(self, X):
        X_arr = X.values if isinstance(X, pd.DataFrame) else X
        return self._pca.transform(X_arr)

    def get_selected_features(self):
        return [f"PC{i+1}" for i in range(self.n_components_)]


# ---------------------------------------------------------------------------
# Correlation-based
# ---------------------------------------------------------------------------
class CorrelationSelector(BaseEstimator, TransformerMixin):
    """
    Two-stage Pearson correlation filter.

    Stage 1 – relevance   : keep features with |corr(f, y)| > target_threshold
    Stage 2 – redundancy  : iterate over features sorted by |corr(f, y)|
                            (high → low); discard f_i if it is highly
                            correlated with any already-accepted feature
                            (|corr(f_i, f_j)| ≥ inter_threshold).

    This avoids multicollinearity while preserving informative features.
    """

    def __init__(self,
                 target_threshold: float = CORR_TARGET_THRESHOLD,
                 inter_threshold:  float = CORR_INTER_THRESHOLD):
        self.target_threshold = target_threshold
        self.inter_threshold  = inter_threshold

    def fit(self, X, y):
        if isinstance(X, pd.DataFrame):
            feature_names = X.columns.tolist()
            X_arr = X.values
        else:
            feature_names = [str(i) for i in range(X.shape[1])]
            X_arr = X

        y_arr = y.values if isinstance(y, pd.Series) else y

        n_feat = X_arr.shape[1]

        # Stage 1: target correlation
        target_corr = np.array([
            abs(np.corrcoef(X_arr[:, i], y_arr)[0, 1])
            for i in range(n_feat)
        ])
        # Guard against NaN (constant features already removed upstream)
        target_corr = np.nan_to_num(target_corr)

        relevant_idx = np.where(target_corr >= self.target_threshold)[0]

        if len(relevant_idx) == 0:
            # Fallback: keep top-10 most correlated
            relevant_idx = np.argsort(target_corr)[-10:]

        # Sort by |corr with target| descending
        relevant_idx = relevant_idx[np.argsort(target_corr[relevant_idx])[::-1]]

        # Stage 2: inter-feature redundancy removal
        selected_idx = []
        X_rel = X_arr[:, relevant_idx]

        for pos, idx in enumerate(relevant_idx):
            if not selected_idx:
                selected_idx.append(idx)
                continue
            # Correlation with already-selected features
            sel_cols = X_arr[:, selected_idx]
            corr_with_selected = np.array([
                abs(np.corrcoef(X_arr[:, idx], sel_cols[:, j])[0, 1])
                for j in range(sel_cols.shape[1])
            ])
            corr_with_selected = np.nan_to_num(corr_with_selected)
            if corr_with_selected.max() < self.inter_threshold:
                selected_idx.append(idx)

        self._selected_idx    = np.array(selected_idx, dtype=int)
        self.feature_names_in_ = feature_names
        self.selected_features_ = [feature_names[i] for i in self._selected_idx]
        return self

    def transform(self, X):
        if isinstance(X, pd.DataFrame):
            return X.iloc[:, self._selected_idx]
        return X[:, self._selected_idx]

    def get_selected_features(self):
        return getattr(self, "selected_features_", None)


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------
def get_feature_selector(method: str) -> BaseEstimator:
    """Return an *unfitted* feature selector for the given method name."""
    selectors = {
        "whole":       WholeSelector,
        "boruta":      BorutaSelector,
        "pca":         PCASelector,
        "correlation": CorrelationSelector,
    }
    if method not in selectors:
        raise ValueError(f"Unknown feature selection method: {method!r}. "
                         f"Choose from {list(selectors)}")
    return selectors[method]()
