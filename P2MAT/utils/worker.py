######################################################
# Developer : Methun Kamruzzaman
# Date      : May 2, 2026
# Purpose   : Background QThread worker that runs property prediction off
#             the main UI thread, emitting Qt signals for progress and results.
######################################################
"""
utils/worker.py
---------------
Background QThread worker for property prediction.

Separates the blocking computation from the PyQt5 main (UI) thread.
The worker receives a list of SMILES strings, runs
:class:`~utils.mspp.MaterialStructuralPropertyPrediction` for each one, and
communicates progress and results back to the main thread via Qt signals.

This module intentionally contains no UI code.  All GUI updates are handled by
slot methods in ``p2mat.MatPropPred`` that connect to the signals below.

Signals emitted by :class:`PredictionWorker`:

- ``progress(int)``         — current count of processed molecules (1..N)
- ``result(object)``        — final ``pd.DataFrame`` of predictions (or ``None``)
- ``invalid_smiles(list)``  — SMILES strings that failed RDKit validation
- ``error(str)``            — human-readable error message for unexpected exceptions
"""
from __future__ import annotations

from PyQt5.QtCore import QThread, pyqtSignal

from utils.mspp import MaterialStructuralPropertyPrediction


class PredictionWorker(QThread):
    """Run :class:`~utils.mspp.MaterialStructuralPropertyPrediction` in a
    background thread so the PyQt5 event loop stays responsive.

    Each SMILES is processed independently.  Errors on individual molecules are
    caught and emitted via :attr:`error` rather than aborting the whole batch.

    Example::

        worker = PredictionWorker(smiles_list, props={"MP": True})
        worker.progress.connect(progress_bar.setValue)
        worker.result.connect(on_result)
        worker.start()
    """

    # ---- Qt signals -------------------------------------------------------
    progress = pyqtSignal(int)
    """Emitted after each molecule is processed.  Value = molecules done so far."""

    valid_count = pyqtSignal(int)
    """Emitted once after SMILES validation with the count of valid molecules.
    The main thread uses this to correct the progress bar maximum."""

    status = pyqtSignal(str)
    """Emitted at each pipeline stage with a human-readable status message."""

    result = pyqtSignal(object)
    """Emitted once when the run completes.  Carries the combined
    ``pd.DataFrame`` of all successful predictions, or ``None`` if every
    molecule failed."""

    invalid_smiles = pyqtSignal(list)
    """Emitted once at the end.  Carries a (possibly empty) list of SMILES
    strings that failed RDKit structural validation."""

    error = pyqtSignal(str)
    """Emitted whenever an unexpected exception occurs for a single molecule."""

    # -----------------------------------------------------------------------

    def __init__(
        self,
        smiles: list[str],
        props: dict,
        batch: int = 128,
        verbose: int = 0,
        parent=None,
    ) -> None:
        """
        Args:
            smiles:  List of SMILES strings to process (one per molecule).
            props:   Property-enable dict, e.g. ``{"MP": True, "BP": False}``.
            batch:   PaDEL batch size passed to
                     :class:`~utils.mspp.MaterialStructuralPropertyPrediction`.
            verbose: Verbosity level forwarded to the prediction class.
            parent:  Optional Qt parent object.
        """
        super().__init__(parent)
        self._smiles = smiles
        self._props = props
        self._batch = batch
        self._verbose = verbose

    def run(self) -> None:
        """Entry point executed in the background thread.

        Passes all SMILES to a single MSPP instance so PaDEL starts once and
        models are loaded once, regardless of input size.  Progress and status
        messages are reported back to the main thread via Qt signals.

        The main thread must never call this method directly — use
        :meth:`QThread.start` instead.
        """
        try:
            mspp = MaterialStructuralPropertyPrediction(
                props=self._props,
                batch=self._batch,
                verbose=self._verbose,
                progress_callback=lambda n: self.progress.emit(n),
                status_callback=lambda msg: self.status.emit(msg),
                valid_count_callback=lambda n: self.valid_count.emit(n),
            )
            wrong_smiles = mspp.prediction_from_smile(self._smiles)
        except Exception as exc:
            self.error.emit(str(exc))
            self.invalid_smiles.emit([])
            self.result.emit(None)
            return

        self.status.emit("Prediction complete.")
        self.invalid_smiles.emit(wrong_smiles)
        self.result.emit(mspp.get_data())
