######################################################
# Developer : Methun Kamruzzaman
# Date      : April 9, 2026
# Purpose   : PyTorch DNN wrapped as a scikit-learn estimator for melting point
#             (2D) QSAR regression.
######################################################
"""
dnn_model.py
------------
PyTorch Deep Neural Network for QSAR regression, wrapped in a
scikit-learn–compatible estimator so it can be used with BayesSearchCV.

Architecture
------------
  Input → [Linear → BatchNorm1d → Activation → Dropout] × n_layers → Linear(1)

Hardware
--------
Uses Apple MPS (Metal Performance Shaders) when available on M2 Max,
falling back to CPU automatically.

Hyperparameters (tunable via BayesSearchCV)
-------------------------------------------
  n_layers      : number of hidden blocks          (2 – 5)
  hidden_size   : units in each hidden layer        (64 – 512)
  dropout_rate  : dropout probability               (0.0 – 0.5)
  learning_rate : Adam learning rate                (1e-4 – 1e-2)
  batch_size    : mini-batch size                   (32 – 256)
  weight_decay  : L2 regularisation coefficient     (1e-6 – 1e-2)
  activation    : 'relu' | 'tanh' | 'leaky_relu'
"""

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.utils.validation import check_is_fitted

from src.config import RANDOM_STATE


# ---------------------------------------------------------------------------
# Detect best available device (prefer MPS on Apple Silicon)
# ---------------------------------------------------------------------------
def _get_device() -> torch.device:
    if torch.backends.mps.is_available():
        return torch.device("mps")
    if torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


class _SafeBatchSampler(torch.utils.data.Sampler):
    """
    Batch sampler that reshuffles indices each epoch and merges a singleton
    last batch into the preceding batch to prevent BatchNorm1d crashes.

    Example: 9 samples, batch_size=2
        → shuffled indices split into [2, 2, 2, 2, 1]
        → last singleton merged   → [2, 2, 2, 3]
    """

    def __init__(self, n_samples: int, batch_size: int):
        self.n_samples  = n_samples
        self.batch_size = batch_size

    def __iter__(self):
        indices = torch.randperm(self.n_samples).tolist()
        batches = [
            indices[i : i + self.batch_size]
            for i in range(0, self.n_samples, self.batch_size)
        ]
        # Merge a singleton last batch into the preceding one
        if len(batches) > 1 and len(batches[-1]) == 1:
            batches[-2].extend(batches[-1])
            batches = batches[:-1]
        yield from batches

    def __len__(self) -> int:
        n = self.n_samples // self.batch_size
        remainder = self.n_samples % self.batch_size
        # If remainder==1 it is merged into the last full batch, not a new one
        return n if remainder == 1 else n + (1 if remainder > 0 else 0)


def _make_loader(dataset, batch_size: int, n_samples: int) -> DataLoader:
    """Return a DataLoader backed by _SafeBatchSampler."""
    sampler = _SafeBatchSampler(n_samples, batch_size)
    return DataLoader(dataset, batch_sampler=sampler)


# ---------------------------------------------------------------------------
# PyTorch Module
# ---------------------------------------------------------------------------
class _QSARNet(nn.Module):
    """Variable-depth MLP with BatchNorm and Dropout."""

    _ACT_MAP = {
        "relu":       nn.ReLU,
        "tanh":       nn.Tanh,
        "leaky_relu": nn.LeakyReLU,
    }

    def __init__(self, input_dim: int, n_layers: int, hidden_size: int,
                 dropout_rate: float, activation: str):
        super().__init__()
        act_cls = self._ACT_MAP.get(activation, nn.ReLU)

        layers = []
        in_dim = input_dim
        for _ in range(n_layers):
            layers += [
                nn.Linear(in_dim, hidden_size),
                nn.BatchNorm1d(hidden_size),
                act_cls(),
                nn.Dropout(p=dropout_rate),
            ]
            in_dim = hidden_size
        layers.append(nn.Linear(in_dim, 1))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x).squeeze(-1)


# ---------------------------------------------------------------------------
# sklearn wrapper
# ---------------------------------------------------------------------------
class DNNRegressor(BaseEstimator, RegressorMixin):
    """
    sklearn-compatible wrapper around _QSARNet.

    Implements fit / predict / score for use with BayesSearchCV.
    score() returns negative MAE (higher = better, consistent with
    'neg_mean_absolute_error' convention).
    """

    def __init__(
        self,
        n_layers:      int   = 3,
        hidden_size:   int   = 256,
        dropout_rate:  float = 0.2,
        learning_rate: float = 1e-3,
        batch_size:    int   = 64,
        weight_decay:  float = 1e-4,
        activation:    str   = "relu",
        max_epochs:    int   = 200,
        patience:      int   = 20,
        random_state:  int   = RANDOM_STATE,
    ):
        self.n_layers      = n_layers
        self.hidden_size   = hidden_size
        self.dropout_rate  = dropout_rate
        self.learning_rate = learning_rate
        self.batch_size    = batch_size
        self.weight_decay  = weight_decay
        self.activation    = activation
        self.max_epochs    = max_epochs
        self.patience      = patience
        self.random_state  = random_state

    # ---- internal helpers -----------------------------------------------

    def _to_tensor(self, X, y=None, device=None):
        if hasattr(X, "values"):
            X = X.values
        X_t = torch.tensor(X.astype(np.float32), device=device)
        if y is not None:
            if hasattr(y, "values"):
                y = y.values
            y_t = torch.tensor(y.astype(np.float32), device=device)
            return X_t, y_t
        return X_t

    # ---- sklearn API ----------------------------------------------------

    def fit(self, X, y):
        torch.manual_seed(self.random_state)
        np.random.seed(self.random_state)

        device = _get_device()

        X_t, y_t = self._to_tensor(X, y, device=device)
        n_samples, input_dim = X_t.shape

        self.model_ = _QSARNet(
            input_dim    = input_dim,
            n_layers     = int(self.n_layers),
            hidden_size  = int(self.hidden_size),
            dropout_rate = float(self.dropout_rate),
            activation   = str(self.activation),
        ).to(device)

        optimiser = torch.optim.Adam(
            self.model_.parameters(),
            lr           = float(self.learning_rate),
            weight_decay = float(self.weight_decay),
        )
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimiser, mode="min", factor=0.5, patience=5
        )
        loss_fn = nn.HuberLoss()

        batch_size = min(int(self.batch_size), n_samples)
        dataset    = TensorDataset(X_t, y_t)
        loader     = _make_loader(dataset, batch_size, n_samples)

        best_loss    = float("inf")
        best_weights = None
        no_improve   = 0

        self.model_.train()
        for epoch in range(int(self.max_epochs)):
            epoch_loss = 0.0
            for xb, yb in loader:
                optimiser.zero_grad()
                pred = self.model_(xb)
                loss = loss_fn(pred, yb)
                loss.backward()
                nn.utils.clip_grad_norm_(self.model_.parameters(), 1.0)
                optimiser.step()
                epoch_loss += loss.item() * len(xb)

            epoch_loss /= n_samples
            scheduler.step(epoch_loss)

            if epoch_loss < best_loss - 1e-4:
                best_loss    = epoch_loss
                best_weights = {k: v.cpu().clone()
                                for k, v in self.model_.state_dict().items()}
                no_improve   = 0
            else:
                no_improve += 1
                if no_improve >= int(self.patience):
                    break

        if best_weights is not None:
            self.model_.load_state_dict(
                {k: v.to(device) for k, v in best_weights.items()}
            )

        self.device_   = device
        self.input_dim_ = input_dim
        return self

    def predict(self, X):
        check_is_fitted(self, "model_")
        self.model_.eval()
        device = self.device_
        X_t    = self._to_tensor(X, device=device)

        with torch.no_grad():
            out = self.model_(X_t).cpu().numpy()
        return out

    def score(self, X, y):
        """Returns negative MAE (higher is better)."""
        from sklearn.metrics import mean_absolute_error
        y_pred = self.predict(X)
        if hasattr(y, "values"):
            y = y.values
        return -mean_absolute_error(y, y_pred)
