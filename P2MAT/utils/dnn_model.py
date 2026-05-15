######################################################
# Developer : Methun Kamruzzaman
# Date      : May 2, 2026
# Purpose   : PyTorch DNN regressor (sklearn-compatible) used as a base model
#             inside the stacking ensemble for melting-point prediction.
######################################################
"""
utils/dnn_model.py
------------------
Minimal, self-contained copy of the DNNRegressor used in the MP stacking
ensemble.  Dependency on the training ``src.config`` module has been removed;
the only external constant is ``RANDOM_STATE``.

Exported symbols: ``DNNRegressor``, ``_QSARNet``, ``_get_device``
"""
from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.utils.validation import check_is_fitted
from torch.utils.data import DataLoader, TensorDataset

RANDOM_STATE: int = 42


def _get_device() -> torch.device:
    """Return the best available compute device (MPS > CUDA > CPU)."""
    if torch.backends.mps.is_available():
        return torch.device("mps")
    if torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


class _SafeBatchSampler(torch.utils.data.Sampler):
    """Reshuffles each epoch and merges a singleton last batch to avoid
    BatchNorm1d crashes on batches of size 1."""

    def __init__(self, n_samples: int, batch_size: int) -> None:
        self.n_samples = n_samples
        self.batch_size = batch_size

    def __iter__(self):
        indices = torch.randperm(self.n_samples).tolist()
        batches = [
            indices[i : i + self.batch_size]
            for i in range(0, self.n_samples, self.batch_size)
        ]
        if len(batches) > 1 and len(batches[-1]) == 1:
            batches[-2].extend(batches[-1])
            batches = batches[:-1]
        yield from batches

    def __len__(self) -> int:
        n = self.n_samples // self.batch_size
        remainder = self.n_samples % self.batch_size
        return n if remainder == 1 else n + (1 if remainder > 0 else 0)


def _make_loader(dataset, batch_size: int, n_samples: int) -> DataLoader:
    sampler = _SafeBatchSampler(n_samples, batch_size)
    return DataLoader(dataset, batch_sampler=sampler)


class _QSARNet(nn.Module):
    """Variable-depth MLP with BatchNorm and Dropout."""

    _ACT_MAP = {
        "relu": nn.ReLU,
        "tanh": nn.Tanh,
        "leaky_relu": nn.LeakyReLU,
    }

    def __init__(
        self,
        input_dim: int,
        n_layers: int,
        hidden_size: int,
        dropout_rate: float,
        activation: str,
    ) -> None:
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


class DNNRegressor(BaseEstimator, RegressorMixin):
    """sklearn-compatible wrapper around ``_QSARNet``.

    Implements ``fit`` / ``predict`` for use in stacking pipelines.
    """

    def __init__(
        self,
        n_layers: int = 3,
        hidden_size: int = 256,
        dropout_rate: float = 0.2,
        learning_rate: float = 1e-3,
        batch_size: int = 64,
        weight_decay: float = 1e-4,
        activation: str = "relu",
        max_epochs: int = 200,
        patience: int = 20,
        random_state: int = RANDOM_STATE,
    ) -> None:
        self.n_layers = n_layers
        self.hidden_size = hidden_size
        self.dropout_rate = dropout_rate
        self.learning_rate = learning_rate
        self.batch_size = batch_size
        self.weight_decay = weight_decay
        self.activation = activation
        self.max_epochs = max_epochs
        self.patience = patience
        self.random_state = random_state

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

    def fit(self, X, y):
        torch.manual_seed(self.random_state)
        np.random.seed(self.random_state)

        device = _get_device()
        X_t, y_t = self._to_tensor(X, y, device=device)
        n_samples, input_dim = X_t.shape

        self.model_ = _QSARNet(
            input_dim=input_dim,
            n_layers=int(self.n_layers),
            hidden_size=int(self.hidden_size),
            dropout_rate=float(self.dropout_rate),
            activation=str(self.activation),
        ).to(device)

        optimiser = torch.optim.Adam(
            self.model_.parameters(),
            lr=float(self.learning_rate),
            weight_decay=float(self.weight_decay),
        )
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimiser, mode="min", factor=0.5, patience=5
        )
        loss_fn = nn.HuberLoss()

        batch_size = min(int(self.batch_size), n_samples)
        dataset = TensorDataset(X_t, y_t)
        loader = _make_loader(dataset, batch_size, n_samples)

        best_loss = float("inf")
        best_weights = None
        no_improve = 0

        self.model_.train()
        for _ in range(int(self.max_epochs)):
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
                best_loss = epoch_loss
                best_weights = {k: v.cpu().clone() for k, v in self.model_.state_dict().items()}
                no_improve = 0
            else:
                no_improve += 1
                if no_improve >= int(self.patience):
                    break

        if best_weights is not None:
            self.model_.load_state_dict({k: v.to(device) for k, v in best_weights.items()})

        self.device_ = device
        self.input_dim_ = input_dim
        return self

    def predict(self, X) -> np.ndarray:
        check_is_fitted(self, "model_")
        self.model_.eval()
        device = self.device_
        X_t = self._to_tensor(X, device=device)
        with torch.no_grad():
            out = self.model_(X_t).cpu().numpy()
        return out
