######################################################
# Developer : Methun Kamruzzaman
# Date      : August, 2024
# Purpose   : Compute 2D and 3D molecular descriptors from SMILES strings
#             using PaDEL-Descriptor via the padelpy wrapper.
######################################################
"""
utils/descriptor.py
-------------------
Molecular descriptor computation using PaDEL-Descriptor via ``padelpy``.

PaDEL-Descriptor is a Java-based tool that generates 2D and 3D QSPR
descriptors from SMILES strings.  This module wraps it with proper temp-file
management: a ``.smi`` input file and a ``.csv`` output file are created in
the OS temp directory, used for one PaDEL call, and then deleted.
"""
import logging
import os
import tempfile
import warnings
from collections import OrderedDict
from csv import DictReader
from time import sleep

from padelpy import padeldescriptor

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)

_TIMEOUT_SECONDS: int = 20   # per-molecule timeout passed to PaDEL
_MAX_ATTEMPTS: int = 1       # retry count on RuntimeError from PaDEL


class Descriptor:
    """Wraps PaDEL-Descriptor to compute molecular features from SMILES.

    Two temporary files are created at construction time and are guaranteed to
    be deleted after each :meth:`from_smiles` call, even on failure.

    Example::

        desc = Descriptor()
        rows = desc.from_smiles(["CCO", "c1ccccc1"])
        # rows is a list of OrderedDicts, one per molecule
    """

    def __init__(self) -> None:
        # delete=False so we control the file lifetime explicitly.
        with tempfile.NamedTemporaryFile(suffix=".smi", delete=False, mode="w") as f:
            self.smi_file: str = f.name

        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False, mode="w") as f:
            self.output_csv: str = f.name

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _write_smiles(self, smiles) -> None:
        """Write *smiles* to the temp ``.smi`` input file.

        Args:
            smiles: A single SMILES string or a list of SMILES strings.

        Raises:
            TypeError: If *smiles* is neither a ``str`` nor a ``list``.
        """
        with open(self.smi_file, "w") as fh:
            if isinstance(smiles, str):
                fh.write(smiles)
            elif isinstance(smiles, list):
                fh.write("\n".join(smiles))
            else:
                raise TypeError(
                    f"Expected str or list for smiles, got {type(smiles)}"
                )

    def _run_padel(self) -> None:
        """Invoke PaDEL-Descriptor and raise ``RuntimeError`` on failure.

        Reads from ``self.smi_file`` and writes to ``self.output_csv``.
        Both 2D and 3D descriptors are computed; fingerprints are disabled.

        Raises:
            RuntimeError: If PaDEL fails on all retry attempts.
            KeyboardInterrupt: Re-raised immediately so the user can cancel.
        """
        for attempt in range(_MAX_ATTEMPTS):
            try:
                padeldescriptor(
                    mol_dir=self.smi_file,
                    d_file=self.output_csv,
                    convert3d=True,
                    retain3d=True,
                    d_2d=True,
                    d_3d=True,
                    fingerprints=False,
                    sp_timeout=_TIMEOUT_SECONDS,
                    retainorder=True,
                    maxruntime=-1,
                    threads=-1,
                )
                return  # success
            except RuntimeError:
                if attempt == _MAX_ATTEMPTS - 1:
                    raise
                continue
            except KeyboardInterrupt:
                raise

    def _read_output(self) -> list:
        """Parse the PaDEL output CSV and return rows as a list of dicts.

        Returns:
            A list of ``dict`` objects, one per molecule, mapping descriptor
            name to its string value.
        """
        with open(self.output_csv, "r", encoding="utf-8") as fh:
            return [row for row in DictReader(fh)]

    def _cleanup(self) -> None:
        """Delete both temp files, silently ignoring missing-file errors."""
        for path in (self.smi_file, self.output_csv):
            try:
                os.remove(path)
            except FileNotFoundError:
                pass

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def from_smiles(self, smiles) -> list:
        """Compute PaDEL 2D+3D descriptors for one or more SMILES strings.

        Temp files are always cleaned up — even if PaDEL raises.

        Args:
            smiles: A single SMILES string or a list of SMILES strings.

        Returns:
            A list of ``OrderedDict`` objects (one per molecule) mapping each
            descriptor name to its computed value.  The ``"Name"`` field
            produced by PaDEL is removed from every row.

        Raises:
            RuntimeError: If PaDEL returns fewer rows than the number of input
                molecules, or if PaDEL itself raises after all retry attempts.
            TypeError: If *smiles* is not a ``str`` or ``list``.
        """
        try:
            self._write_smiles(smiles)
            self._run_padel()
            rows = self._read_output()
        except RuntimeError:
            self._cleanup()
            sleep(0.5)
            raise
        except Exception:
            self._cleanup()
            raise
        else:
            self._cleanup()

        expected = 1 if isinstance(smiles, str) else len(smiles)
        if len(rows) != expected:
            raise RuntimeError(
                f"PaDEL returned {len(rows)} row(s) for {expected} input "
                "SMILES. Ensure all input structures are valid."
            )

        for row in rows:
            row.pop("Name", None)

        return rows
