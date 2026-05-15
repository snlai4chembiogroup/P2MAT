######################################################
# Developer : Methun Kamruzzaman
# Date      : May 2, 2026
# Purpose   : Self-contained stacking ensemble predictor that combines four
#             base models (LightGBM, XGBoost, ExtraTrees, DNN) via a Ridge
#             meta-learner for melting-point prediction.
######################################################
"""
utils/stacking_predictor.py
---------------------------
Serialisable stacking ensemble for deployment inside P2MAT.

:class:`StackingPredictor` wraps four base models (LightGBM, XGBoost,
ExtraTrees, DNN) plus a Ridge meta-learner into a single sklearn-compatible
estimator.  It is the object saved as ``bestmodel_stack.sav`` and loaded by
:meth:`~utils.mspp.MaterialStructuralPropertyPrediction._predict_MP`.

DNN weights are stored as CPU numpy arrays so the artifact is
device-independent and can be deserialised on any platform.  The ``_QSARNet``
architecture is reconstructed on the fly at predict time.
"""
from __future__ import annotations

import numpy as np


class ScaleSelectPipeline:
    """Lightweight preprocessor that applies an optional scaler then boolean
    column selection.

    Used when the training pipeline scaled ALL features first (e.g. with a
    QuantileTransformer) and then applied Boruta feature selection on the scaled
    output.  Storing both steps together avoids any dependency on the training
    codebase's ``src.feature_selection.BorutaSelector`` at inference time.

    Attributes:
        _scaler:      A fitted sklearn transformer (or ``None`` for no scaling).
        _mask:        Boolean numpy array selecting the Boruta-confirmed columns
                      from the scaled feature matrix.

    Example::

        pipeline = ScaleSelectPipeline(scaler=qt, support_mask=boruta._support)
        X_ready  = pipeline.transform(X_full)   # scale then select
    """

    def __init__(self, scaler, support_mask: np.ndarray) -> None:
        self._scaler = scaler
        self._mask   = np.asarray(support_mask, dtype=bool)

    def transform(self, X: np.ndarray) -> np.ndarray:
        if self._scaler is not None:
            X = self._scaler.transform(X).astype(np.float32)
        return X[:, self._mask].astype(np.float32)


class StackingPredictor:
    """Stacking ensemble: scaler → base-model predictions → Ridge meta-learner.

    Each base model component is stored as a dict with keys:

    - ``"type"``   : ``"sklearn"`` or ``"dnn"``
    - ``"scaler"`` : fitted ``StandardScaler`` (or ``None``)
    - ``"model"``  : fitted estimator (sklearn models only)
    - ``"params"`` : DNN hyperparameter dict (DNN only)
    - ``"weights"``: ``{layer: np.ndarray}`` state dict (DNN only)

    Example::

        predictor = StackingPredictor(base_components, meta_learner)
        preds = predictor.predict(X_df[features])
    """

    def __init__(
        self,
        base_components: list[dict],
        meta_learner,
    ) -> None:
        """
        Args:
            base_components: List of per-base-model config dicts (see class doc).
            meta_learner:    Fitted ``sklearn.linear_model.Ridge`` meta-learner.
        """
        self.base_components = base_components
        self.meta_learner = meta_learner

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _apply_scaler(self, comp: dict, X: np.ndarray) -> np.ndarray:
        """Scale X if a fitted scaler is present, otherwise return X unchanged."""
        scaler = comp.get("scaler")
        return scaler.transform(X).astype(np.float32) if scaler is not None else X.astype(np.float32)

    def _predict_sklearn(self, comp: dict, X: np.ndarray) -> np.ndarray:
        X_s = self._apply_scaler(comp, X)
        return np.asarray(comp["model"].predict(X_s), dtype=np.float32)

    def _predict_dnn(self, comp: dict, X: np.ndarray) -> np.ndarray:
        """Reconstruct the DNN from stored weights and run inference."""
        import torch
        from utils.dnn_model import _QSARNet, _get_device

        X_s = self._apply_scaler(comp, X)
        params = comp["params"]
        weights_np: dict[str, np.ndarray] = comp["weights"]

        device = _get_device()
        model = _QSARNet(
            input_dim=X_s.shape[1],
            n_layers=int(params["n_layers"]),
            hidden_size=int(params["hidden_size"]),
            dropout_rate=float(params["dropout_rate"]),
            activation=str(params["activation"]),
        ).to(device)
        state_dict = {k: torch.tensor(v, device=device) for k, v in weights_np.items()}
        model.load_state_dict(state_dict)
        model.eval()

        with torch.no_grad():
            X_t = torch.tensor(X_s, device=device)
            return model(X_t).cpu().numpy().astype(np.float32)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def predict(self, X) -> np.ndarray:
        """Predict using the full stacking ensemble.

        Args:
            X: Feature matrix (DataFrame or ndarray) with the 1 141 descriptor
               columns expected by the base models.

        Returns:
            1-D numpy array of predicted melting points (K).
        """
        if hasattr(X, "values"):
            X = X.values
        X = np.asarray(X, dtype=np.float32)

        col_preds: list[np.ndarray] = []
        for comp in self.base_components:
            if comp["type"] == "dnn":
                p = self._predict_dnn(comp, X)
            else:
                p = self._predict_sklearn(comp, X)
            col_preds.append(p)

        meta_X = np.column_stack(col_preds)
        return self.meta_learner.predict(meta_X)
