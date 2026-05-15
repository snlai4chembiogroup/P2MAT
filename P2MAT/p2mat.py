#!/usr/bin/env python
######################################################
# Developer : Methun Kamruzzaman
# Date      : August, 2024
# Purpose   : Main GUI application entry point for P2MAT — accepts SMILES
#             input, triggers property prediction, and displays/exports results.
######################################################
"""
p2mat.py
--------
Entry point and main GUI for P2MAT — Material Property Prediction.

Presents a PyQt5 window where the user can:
1. Type one or more SMILES strings (one per line).
2. Tick any combination of property checkboxes populated from ``config.json``.
3. Click **Predict properties** to run predictions in a background thread
   (no UI freezing).
4. Inspect the results table (one column per selected property) and export
   to CSV.

The window delegates all computation to :class:`utils.worker.PredictionWorker`
so the Qt event loop is never blocked.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from PyQt5.Qt import Qt, QTimer
from PyQt5.QtGui import QIcon, QPixmap
from PyQt5.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QCheckBox,
    QDialog,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from utils import helper
from utils.worker import PredictionWorker

# Maps config.json property keys → human-readable combo-box labels.
_LABEL_MAP: dict[str, str] = {
    "MP": "Melting point (K)",
    "BP": "Boiling point (K)",
    "HC": "Heat capacity (J/K)",
    "dH": "Heat of hydrogenation",
    "h2_uptake": "H₂ uptake",
}

class Popup(QDialog):
    """A minimal modal dialog used for informational pop-up messages."""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.setLayout(QVBoxLayout())
        self.label = QLabel(self)
        self.layout().addWidget(self.label)


class MatPropPred(QWidget):
    """Main application window for P2MAT.

    Responsibilities:
    - Build and lay out all Qt widgets.
    - Load ``config.json`` to determine which properties are enabled.
    - Launch :class:`~utils.worker.PredictionWorker` on button click.
    - Receive worker signals and update the UI accordingly.
    - Export prediction results to CSV.
    """

    def __init__(
        self,
        width: int = 800,
        height: int = 600,
        batch_size: int = 128,
        verbose: int = 0,
    ) -> None:
        """
        Args:
            width:      Initial window width in pixels.
            height:     Initial window height in pixels.
            batch_size: PaDEL descriptor batch size forwarded to the worker.
            verbose:    Verbosity level forwarded to the prediction classes.
        """
        super().__init__()
        self.setWindowTitle("Material property prediction")
        self.resize(width, height)
        self.setWindowIcon(QIcon("icon.png"))

        self._width = width
        self._height = height
        self._batch = batch_size
        self._verbose = verbose
        self._config: dict = helper.load_config()
        self._worker: PredictionWorker | None = None
        self.prop_data: pd.DataFrame | None = None

        self._create_UI()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _create_UI(self) -> None:
        """Build and arrange all widgets.  Called once from ``__init__``."""
        # --- Logo -------------------------------------------------------
        logo_lbl = QLabel("")
        logo_path = helper.resource_path(os.path.join("design", "snl-logo-inline.svg"))
        print(f"Logo path: {logo_path}")
        pixmap = QPixmap(logo_path)
        logo_lbl.setPixmap(pixmap)
        logo_lbl.resize(pixmap.width(), pixmap.height())

        # --- SMILES input -----------------------------------------------
        smile_lbl = QLabel("Enter a SMILES string\n(Multiple SMILEs in new line)")
        self.smile_txt = QTextEdit()

        # --- Property checkboxes ----------------------------------------
        # Only properties marked true in config.json are shown, all pre-checked.
        # Each checkbox is wired to _guard_min_one_checked so the last checked
        # box cannot be unchecked.
        self._prop_checkboxes: dict[str, QCheckBox] = {}
        config_props = self._config.get("properties", {})
        for key, label in _LABEL_MAP.items():
            if config_props.get(key, False):
                cb = QCheckBox(label)
                cb.setChecked(True)
                cb.stateChanged.connect(self._guard_min_one_checked)
                self._prop_checkboxes[key] = cb

        # --- Status label + Progress bar --------------------------------
        self.status_lbl = QLabel("")
        self.status_lbl.setAlignment(Qt.AlignLeft)
        self.pbar = QProgressBar()

        # --- Buttons ----------------------------------------------------
        self.button1 = QPushButton("Predict properties")
        self.button1.clicked.connect(self._prop_pred_btn)

        self.button2 = QPushButton("Save data")
        self.button2.clicked.connect(self._save_file)

        # --- Result widgets ---------------------------------------------
        self.result_lbl = QLabel("Result")
        self.result_details_lbl = QPlainTextEdit("")
        self.result_details_lbl.setReadOnly(True)

        self.result_table = QTableWidget()
        self.result_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.result_table.setSelectionMode(QAbstractItemView.NoSelection)
        self.result_table.setFocusPolicy(Qt.NoFocus)

        # Initially hidden
        self.result_lbl.setHidden(True)
        self.result_details_lbl.setHidden(True)
        self.result_table.setHidden(True)
        self.button2.setHidden(True)
        self.status_lbl.setHidden(True)
        self.pbar.setHidden(True)

        # --- Layout assembly --------------------------------------------
        row_logo = QHBoxLayout()
        row_logo.addWidget(logo_lbl)

        row_smile = QHBoxLayout()
        row_smile.addWidget(smile_lbl)
        row_smile.addWidget(self.smile_txt)

        prop_group = QGroupBox("Properties to predict")
        row_prop = QHBoxLayout()
        for cb in self._prop_checkboxes.values():
            row_prop.addWidget(cb)
        row_prop.addStretch()
        prop_group.setLayout(row_prop)

        row_button_1 = QHBoxLayout()
        row_button_1.addWidget(self.button1, alignment=Qt.AlignRight)

        row_pb = QVBoxLayout()
        row_pb.addWidget(self.status_lbl)
        row_pb.addWidget(self.pbar)

        result_form = QFormLayout()
        result_form.addRow(self.result_lbl)
        result_form.addRow(self.result_details_lbl)
        result_form.addRow(self.result_table)

        row_button_2 = QHBoxLayout()
        row_button_2.addWidget(self.button2, alignment=Qt.AlignRight)

        master_layout = QVBoxLayout()
        master_layout.addLayout(row_logo)
        master_layout.addLayout(row_smile)
        master_layout.addWidget(prop_group)
        master_layout.addLayout(row_button_1)
        master_layout.addLayout(row_pb)
        master_layout.addLayout(result_form)
        master_layout.addLayout(row_button_2)
        self.setLayout(master_layout)

    # ------------------------------------------------------------------
    # Helper methods
    # ------------------------------------------------------------------

    def _get_smiles(self) -> list[str]:
        """Return non-blank SMILES strings from the text input area."""
        return [s for s in self.smile_txt.toPlainText().split("\n") if s.strip()]

    def _guard_min_one_checked(self, _) -> None:
        """Schedule an enforcement check after the current click event finishes.

        Connected to the ``stateChanged(int)`` signal of every checkbox.
        Using ``QTimer.singleShot(0, ...)`` defers the check until Qt has fully
        processed the click, preventing the re-check from being overridden by
        the in-flight mouse-release event.
        """
        QTimer.singleShot(0, self._enforce_min_one_checked)

    def _enforce_min_one_checked(self) -> None:
        """Re-check Melting point if all checkboxes are currently unchecked.

        Runs as a deferred callback so the Qt event loop has finished applying
        the user's last checkbox interaction before we inspect or change state.
        """
        if any(cb.isChecked() for cb in self._prop_checkboxes.values()):
            return

        # Fall back to MP, or the first available property if MP is absent.
        fallback = self._prop_checkboxes.get("MP") or next(
            iter(self._prop_checkboxes.values()), None
        )
        if fallback:
            fallback.blockSignals(True)
            fallback.setChecked(True)
            fallback.blockSignals(False)

    def _build_props_dict(self) -> dict:
        """Build the *props* dict to pass to the worker from the checkboxes.

        Returns a dict mapping each property key to ``True`` if its checkbox
        is currently checked, ``False`` otherwise.
        """
        return {key: cb.isChecked() for key, cb in self._prop_checkboxes.items()}

    # ------------------------------------------------------------------
    # Result table
    # ------------------------------------------------------------------

    def _create_prop_table(self) -> None:
        """Populate the result table from ``self.prop_data``."""
        if self.prop_data is None:
            return

        self.result_table.setHidden(False)
        self.result_table.setRowCount(self.prop_data.shape[0])
        self.result_table.setColumnCount(self.prop_data.shape[1])

        # Map raw column names to display headers.
        _col_labels: dict[str, str] = {
            "smiles": "SMILES",
            "MP_pred_K": "Melting point (K)",
            "BP_pred_K": "Boiling point (K)",
            "HC_pred_JK": "Heat capacity (J/K)",
            "dH_pred_JK": "Heat of hydrogenation",
            "Saturated_smiles": "Saturated SMILES",
            "H2_uptake": "H₂ uptake",
        }
        for idx, col in enumerate(self.prop_data.columns):
            header = _col_labels.get(col, col)
            self.result_table.setHorizontalHeaderItem(idx, QTableWidgetItem(header))

        # Numeric prediction columns (all except SMILES strings).
        str_cols = {"smiles", "Saturated_smiles"}
        for row in range(self.prop_data.shape[0]):
            for col_idx, col in enumerate(self.prop_data.columns):
                value = self.prop_data.iloc[row, col_idx]
                if col not in str_cols and isinstance(value, (int, float)):
                    cell_text = str(np.round(value, 3))
                else:
                    cell_text = str(value)
                self.result_table.setItem(row, col_idx, QTableWidgetItem(cell_text))

        self.result_table.horizontalHeader().setStretchLastSection(True)
        self.result_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

    # ------------------------------------------------------------------
    # Button handlers
    # ------------------------------------------------------------------

    def _prop_pred_btn(self) -> None:
        """Handle the *Predict properties* button click.

        Validates that there is input, then starts the
        :class:`~utils.worker.PredictionWorker` background thread.
        """
        # Reset result area
        self.result_details_lbl.setPlainText("")
        self.result_details_lbl.setHidden(True)
        self.result_table.setHidden(True)
        self.result_lbl.setHidden(True)
        self.button2.setHidden(True)
        self.pbar.setHidden(True)
        self.prop_data = None

        smiles = self._get_smiles()
        if not smiles:
            return

        self.pbar.setMaximum(len(smiles))
        self.pbar.setValue(0)
        self.status_lbl.setText("Starting…")
        self.status_lbl.setHidden(False)
        self.pbar.setHidden(False)
        self.resize(self._width, self._height + 100)
        self.button1.setEnabled(False)

        props = self._build_props_dict()
        self._worker = PredictionWorker(
            smiles, props, batch=self._batch, verbose=self._verbose
        )
        self._worker.progress.connect(self._on_progress)
        self._worker.valid_count.connect(self.pbar.setMaximum)
        self._worker.status.connect(self._on_status)
        self._worker.result.connect(self._on_result)
        self._worker.invalid_smiles.connect(self._on_invalid_smiles)
        self._worker.error.connect(self._on_error)
        self._worker.finished.connect(self._on_worker_finished)
        self._worker.start()

    def _save_file(self) -> None:
        """Handle the *Save data* button click.

        Opens a save-file dialog and writes ``self.prop_data`` to the chosen
        CSV path.  The ``.csv`` extension is appended automatically if missing.
        """
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getSaveFileName(
            self, "Save File", "", "CSV Files (*.csv)", options=options
        )

        if not file_name:
            return

        if not file_name.lower().endswith(".csv"):
            file_name += ".csv"

        if self.prop_data is not None and not self.prop_data.empty:
            self.prop_data.to_csv(file_name, index=False)
            print(f"File saved to {file_name}")
        else:
            print("No data to save.")

    # ------------------------------------------------------------------
    # Worker signal slots
    # ------------------------------------------------------------------

    def _on_progress(self, value: int) -> None:
        """Update the progress bar as each molecule finishes."""
        self.pbar.setValue(value)

    def _on_status(self, message: str) -> None:
        """Update the status label with the latest pipeline stage message."""
        self.status_lbl.setText(message)

    def _on_result(self, data) -> None:
        """Store the completed results DataFrame."""
        self.prop_data = data

    def _on_invalid_smiles(self, smiles: list) -> None:
        """Display any SMILES that failed RDKit validation."""
        if smiles:
            self.result_lbl.setHidden(False)
            self.result_details_lbl.setHidden(False)
            self.result_details_lbl.appendPlainText(
                "The following SMILES are in an invalid format:"
            )
            for i, s in enumerate(smiles):
                self.result_details_lbl.appendPlainText(f"  {i + 1}) {s}")

    def _on_error(self, message: str) -> None:
        """Print unexpected per-molecule errors to the console."""
        print(f"Prediction error: {message}")

    def _on_worker_finished(self) -> None:
        """Re-enable the button and display results when the worker is done."""
        self.button1.setEnabled(True)
        self.status_lbl.setHidden(True)

        if self.prop_data is not None and not self.prop_data.empty:
            self.button2.setHidden(False)
            self.result_lbl.setHidden(False)
            self.resize(self._width, self._height + 300)
            self._create_prop_table()


# ---------------------------------------------------------------------------
# Application entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon("icon.png"))

    css_file = helper.resource_path(os.path.join("design", "myapp.css"))
    print(f"CSS file path: {css_file}")
    app.setStyleSheet(helper.read_and_replace_css(css_file))

    window = MatPropPred()
    window.show()
    sys.exit(app.exec_())
