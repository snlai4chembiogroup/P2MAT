######################################################
# Developer : Methun Kamruzzaman
# Date      : August, 2024
# Purpose   : Core prediction pipeline — validates SMILES, computes PaDEL
#             descriptors, runs pre-trained ML models, and optionally
#             calculates H2 uptake via bond saturation.
######################################################
"""
utils/mspp.py
-------------
Core material property prediction pipeline.

:class:`MaterialStructuralPropertyPrediction` orchestrates the full workflow:

1. Validate SMILES with RDKit via :class:`~utils.smile.SMILE`.
2. Compute ~1 240 PaDEL molecular descriptors in batches via
   :class:`~utils.descriptor.Descriptor`.
3. Load pre-trained joblib models and run inference for each enabled property.
4. Optionally compute H₂ uptake by chemically saturating the molecule.
5. Optionally render molecule PNG images via :class:`~utils.depictsmile.DepictSmile`.

Supported properties (toggled via *props* or ``config.json``):

- **MP** — Melting point (K)
- **BP** — Boiling point (K)
- **HC** — Heat capacity (J/K)
- **dH** — Heat of hydrogenation (J/K)
- **h2_uptake** — Moles of H₂ required for full saturation (rule-based)
"""
import functools
import logging
import os
import warnings

import joblib
import numpy as np
import pandas as pd
from rdkit import Chem

from utils import helper
from utils.depictsmile import DepictSmile
from utils.descriptor import Descriptor
from utils.smile import SMILE

warnings.filterwarnings("ignore")
pd.options.mode.chained_assignment = None

# ---------------------------------------------------------------------------
# File logger — writes descriptor failures and skipped molecules to p2mat.log
# ---------------------------------------------------------------------------
_log_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "p2mat.log")
_handler  = logging.FileHandler(_log_path, encoding="utf-8")
_handler.setFormatter(logging.Formatter("%(asctime)s  %(levelname)-8s  %(message)s",
                                         datefmt="%Y-%m-%d %H:%M:%S"))
logger = logging.getLogger("p2mat")
logger.setLevel(logging.DEBUG)
if not logger.handlers:
    logger.addHandler(_handler)


@functools.lru_cache(maxsize=None)
def _load_model(model_path: str):
    """Load and cache a joblib model by path — deserialised once per process."""
    return joblib.load(model_path)


class MaterialStructuralPropertyPrediction:
    """Predict thermophysical properties of molecules from SMILES strings.

    The class is designed to be instantiated once per prediction call (as done
    in the GUI worker) so that ``df_mol`` state never leaks between calls.

    Example::

        mspp = MaterialStructuralPropertyPrediction()
        invalid = mspp.prediction_from_smile(["CCO", "c1ccccc1"])
        df = mspp.get_data()   # DataFrame with smiles + MP_pred_K columns
    """

    def __init__(
        self,
        props: dict = None,
        imgpath: str = None,
        batch: int = 128,
        verbose: int = 0,
        progress_callback=None,
        status_callback=None,
        valid_count_callback=None,
    ) -> None:
        """
        Args:
            props: Dict mapping property keys to booleans, e.g.
                ``{"MP": True, "BP": False, ...}``.  If ``None``, the value is
                loaded from ``config.json`` via :func:`~utils.helper.load_config`.
            imgpath: Directory path for molecule PNG files.  Only used when
                ``h2_uptake`` is enabled.  The directory is created if absent.
            batch: Number of SMILES sent to PaDEL in a single call.
            verbose: Logging verbosity (0 = silent).
            progress_callback: Optional callable(n: int) invoked each time a
                batch (or single molecule) finishes, where n is the cumulative
                count of molecules processed so far.
            status_callback: Optional callable(msg: str) invoked at each
                pipeline stage with a human-readable status string.
            valid_count_callback: Optional callable(n: int) invoked once after
                SMILES validation with the count of valid molecules, so the
                caller can adjust the progress bar maximum accordingly.
        """
        if props is None:
            props = helper.load_config().get("properties", {"MP": True})

        self.props = props
        self._smile = SMILE(verbose)
        self._model_path = helper.resource_path(os.path.join("include", "best_models"))
        self._feature_path = helper.resource_path(os.path.join("include", "feature"))
        self._batch = batch
        self._verbose = verbose
        self._on_progress = progress_callback
        self._on_status = status_callback
        self._on_valid_count = valid_count_callback
        self._depict: DepictSmile | None = None
        self.df_mol: pd.DataFrame | None = None
        self._padel_failed: list[str] = []
        self._processed_count: int = 0
        self._n_batches: int = 0

        if imgpath is not None:
            self._depict = DepictSmile(
                save_file_path=self._ensure_dir(imgpath)
            )

    # ------------------------------------------------------------------
    # Status helper
    # ------------------------------------------------------------------

    def _status(self, msg: str) -> None:
        """Emit *msg* via the status callback if one was provided."""
        if self._on_status is not None:
            self._on_status(msg)

    # ------------------------------------------------------------------
    # Directory / file helpers
    # ------------------------------------------------------------------

    def _ensure_dir(self, path: str) -> str:
        """Create *path* as a directory (including parents) and return it."""
        os.makedirs(path, exist_ok=True)
        return path

    def _read_features(self, filepath: str) -> list:
        """Read a newline-delimited feature-name file.

        Args:
            filepath: Absolute path to the feature list text file.

        Returns:
            A list of feature name strings (one per line, blank lines skipped).

        Raises:
            FileNotFoundError: If *filepath* does not exist.
        """
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Feature file not found: {filepath}")
        with open(filepath, "r") as fh:
            return [line.rstrip("\n") for line in fh if line.strip()]

    # ------------------------------------------------------------------
    # Descriptor computation
    # ------------------------------------------------------------------

    def _download_batch(self, smiles: list, start: int, end: int) -> None:
        """Compute PaDEL descriptors for ``smiles[start:end]`` and append to
        ``self.df_mol``.

        On batch failure the range is recursively halved (divide-and-conquer)
        until individual molecules are attempted.  Molecules that still fail
        individually are skipped with a warning.

        Args:
            smiles: The full list of validated SMILES strings.
            start:  Inclusive start index into *smiles*.
            end:    Exclusive end index into *smiles*.
        """
        if start >= end:
            return

        if self._verbose > 0:
            print(f"Fetching descriptors [{start}, {end})")

        # PaDEL call isolated in its own try block so that divide-and-conquer
        # only triggers on descriptor failures, not on downstream pandas errors.
        try:
            desc = Descriptor()
            rows = desc.from_smiles(smiles[start:end])
            del desc  # release temp file handles immediately
        except Exception as e:
            if end - start == 1:
                # Single molecule still fails — record and skip it.
                logger.warning("Skipped molecule '%s': %s", smiles[start], e)
                self._padel_failed.append(smiles[start])
                self._processed_count += 1
                if self._on_progress is not None:
                    self._on_progress(self._processed_count)
                return

            logger.warning(
                "Descriptor batch [%d:%d] (%d molecules) failed — splitting: %s",
                start, end, end - start, e,
            )
            mid = start + (end - start) // 2
            self._download_batch(smiles, start, mid)
            self._download_batch(smiles, mid, end)
            return

        # PaDEL succeeded — build DataFrame outside the PaDEL try so any
        # pandas error here surfaces as a real exception, not a split trigger.
        batch_df = pd.DataFrame(rows)
        batch_df["smiles"] = smiles[start:end]
        self.df_mol = (
            batch_df
            if self.df_mol is None
            else pd.concat([self.df_mol, batch_df], ignore_index=True)
        )
        self._processed_count += end - start
        if self._on_progress is not None:
            self._on_progress(self._processed_count)

    def _compute_features(self, smiles: list) -> None:
        """Compute descriptors for all SMILES in ``self._batch``-sized batches.

        Results are accumulated in ``self.df_mol``.  Descriptor columns that
        contain empty strings are coerced to ``0.0`` and cast to ``float32``.

        Args:
            smiles: List of validated SMILES strings.

        Raises:
            RuntimeError: If no descriptors could be computed for any molecule.
        """
        self.df_mol = None
        self._padel_failed = []
        self._processed_count = 0
        self._n_batches = (len(smiles) + self._batch - 1) // self._batch

        self._status(
            f"Splitting {len(smiles)} molecule(s) into "
            f"{self._n_batches} batch(es) of up to {self._batch}…"
        )

        for batch_num, start in enumerate(range(0, len(smiles), self._batch), 1):
            end = min(start + self._batch, len(smiles))
            self._status(
                f"Fetching descriptors — batch {batch_num} of {self._n_batches} "
                f"({end - start} molecule(s))…"
            )
            self._download_batch(smiles, start, end)
            if batch_num < self._n_batches:
                self._status(
                    f"Batch {batch_num} done — preparing batch {batch_num + 1} "
                    f"of {self._n_batches}…"
                )

        if self.df_mol is None or self.df_mol.empty:
            raise RuntimeError(
                "No descriptors could be computed for the provided SMILES."
            )

        for col in self.df_mol.columns:
            if col != "smiles":
                self.df_mol[col] = (
                    self.df_mol[col].replace("", "0.0").astype(np.float32)
                )

    # ------------------------------------------------------------------
    # Feature validation
    # ------------------------------------------------------------------

    def _validate_features(self, df: pd.DataFrame) -> pd.DataFrame:
        """Replace ``inf``/``-inf`` with ``NaN``, then fill ``NaN`` with ``0``.

        Args:
            df: A DataFrame slice containing only the model's input features.

        Returns:
            A cleaned copy of *df* safe to pass to a scikit-learn model.
        """
        df = df.replace([np.inf, -np.inf], np.nan)
        df = df.fillna(0)
        return df

    # ------------------------------------------------------------------
    # Per-property model inference
    # ------------------------------------------------------------------

    def _predict_property(
        self, model_rel_path: str, feature_filename: str,
        df: pd.DataFrame | None = None,
    ) -> np.ndarray:
        """Generic property predictor: load model → select features → predict.

        Args:
            model_rel_path: Path relative to ``include/best_models/``.
            feature_filename: Filename inside ``include/feature/``.
            df: DataFrame to predict on.  Defaults to ``self.df_mol``.

        Returns:
            A 1-D ``numpy`` array of predictions aligned with *df* rows.

        Raises:
            FileNotFoundError: If the model or feature file is missing.
            KeyError: If a required feature column is absent from *df*.
        """
        if df is None:
            df = self.df_mol
        assert df is not None

        model_path = os.path.join(self._model_path, model_rel_path)
        if not os.path.exists(model_path):
            raise FileNotFoundError(f"Model not found: {model_path}")

        model = _load_model(model_path)
        features = self._read_features(
            os.path.join(self._feature_path, feature_filename)
        )
        fv = self._validate_features(df[features].copy())
        return model.predict(fv)

    # Elements present in CHO and CHON compounds (atomic numbers).
    _CHO_ATOMS  = frozenset({1, 6, 8})     # H, C, O
    _CHON_ATOMS = frozenset({1, 6, 7, 8})  # H, C, O, N

    @staticmethod
    def _classify_elements(smiles: str) -> str:
        """Return the elemental domain of *smiles*: ``'cho'``, ``'chon'``, or
        ``'full'``.

        Classification is based on the set of heavy-atom atomic numbers:
        - Only C/H/O  → ``'cho'``
        - Only C/H/O/N → ``'chon'``
        - Anything else → ``'full'``
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "full"
        elements = frozenset(atom.GetAtomicNum() for atom in mol.GetAtoms())
        if elements <= MaterialStructuralPropertyPrediction._CHO_ATOMS:
            return "cho"
        if elements <= MaterialStructuralPropertyPrediction._CHON_ATOMS:
            return "chon"
        return "full"

    # Maps each domain to its (model file, feature file) pair.
    _MP_DOMAIN_MODELS = {
        "cho":  (os.path.join("MP", "bestmodel_stack_cho.sav"),  "MP_features_cho.txt"),
        "chon": (os.path.join("MP", "bestmodel_stack_chon.sav"), "MP_features_chon.txt"),
        "full": (os.path.join("MP", "bestmodel_stack_full.sav"), "MP_features_full.txt"),
    }

    def _predict_MP(self) -> np.ndarray:
        """Predict melting point (K) using a domain-specific stacking ensemble.

        Each SMILES is classified into 'cho', 'chon', or 'full' by its
        constituent elements, then routed to the model trained on that domain:
        - CHO-only  → bestmodel_stack_cho  (trained on C/H/O compounds)
        - CHON-only → bestmodel_stack_chon (trained on C/H/O/N compounds)
        - Other     → bestmodel_stack_full (trained on the full dataset)
        """
        assert self.df_mol is not None
        smiles_list = self.df_mol["smiles"].tolist()
        domains     = [self._classify_elements(s) for s in smiles_list]
        preds       = np.empty(len(smiles_list), dtype=np.float64)

        for domain, (model_rel, feat_file) in self._MP_DOMAIN_MODELS.items():
            idx = [i for i, d in enumerate(domains) if d == domain]
            if not idx:
                continue
            self._status(
                f"  Predicting MP — {domain.upper()} domain "
                f"({len(idx)} molecule(s))…"
            )
            sub_df      = pd.DataFrame(self.df_mol.iloc[idx]).reset_index(drop=True)
            group_preds = self._predict_property(model_rel, feat_file, df=sub_df)
            for j, i in enumerate(idx):
                preds[i] = group_preds[j]

        return preds

    # Maps each domain to its (model file, feature file) pair for BP.
    _BP_DOMAIN_MODELS = {
        "cho":  (os.path.join("BP", "bestmodel_stack_cho.sav"),  "BP_features_cho.txt"),
        "chon": (os.path.join("BP", "bestmodel_stack_chon.sav"), "BP_features_chon.txt"),
        "full": (os.path.join("BP", "bestmodel_stack_full.sav"), "BP_features_full.txt"),
    }

    def _predict_BP(self) -> np.ndarray:
        """Predict boiling point (K) using a domain-specific stacking ensemble.

        Each SMILES is classified into 'cho', 'chon', or 'full' by its
        constituent elements, then routed to the model trained on that domain:
        - CHO-only  → bestmodel_stack_cho  (trained on C/H/O compounds)
        - CHON-only → bestmodel_stack_chon (trained on C/H/O/N compounds)
        - Other     → bestmodel_stack_full (trained on the full dataset)
        """
        assert self.df_mol is not None
        smiles_list = self.df_mol["smiles"].tolist()
        domains     = [self._classify_elements(s) for s in smiles_list]
        preds       = np.empty(len(smiles_list), dtype=np.float64)

        for domain, (model_rel, feat_file) in self._BP_DOMAIN_MODELS.items():
            idx = [i for i, d in enumerate(domains) if d == domain]
            if not idx:
                continue
            self._status(
                f"  Predicting BP — {domain.upper()} domain "
                f"({len(idx)} molecule(s))…"
            )
            sub_df      = pd.DataFrame(self.df_mol.iloc[idx]).reset_index(drop=True)
            group_preds = self._predict_property(model_rel, feat_file, df=sub_df)
            for j, i in enumerate(idx):
                preds[i] = group_preds[j]

        return preds

    def _predict_HC(self) -> np.ndarray:
        """Predict heat capacity (J/K) for all rows in ``self.df_mol``."""
        return self._predict_property(
            os.path.join("HC", "bestmodel.sav"), "HC_features.txt"
        )

    def _predict_HH(self) -> np.ndarray:
        """Predict heat of hydrogenation for all rows in ``self.df_mol``."""
        return self._predict_property(
            os.path.join("HH", "bestmodel.sav"), "HH_features.txt"
        )

    # ------------------------------------------------------------------
    # H₂ uptake / saturation
    # ------------------------------------------------------------------

    @staticmethod
    def _add_hydrogens(emol, atom_idx: int, count: int) -> None:
        """Add *count* explicit hydrogen atoms bonded to *atom_idx* in *emol*.

        Args:
            emol:     An RDKit ``EditableMol`` object.
            atom_idx: Index of the heavy atom to which H atoms are bonded.
            count:    Number of hydrogen atoms to add.
        """
        for _ in range(count):
            h_idx = emol.AddAtom(Chem.Atom(1))
            emol.AddBond(atom_idx, h_idx, Chem.BondType.SINGLE)

    def _saturate_smiles(self, smiles: str) -> tuple:
        """Convert all double/triple bonds to single bonds and report H₂ needed.

        Aromatic bonds are left unchanged.  The moles-of-H₂ value is computed
        as ``number_of_unsaturated_bonds × 2``, matching the original formula.

        Args:
            smiles: A valid SMILES string.

        Returns:
            A 2-tuple ``(saturated_smiles, moles_of_H2)``.  Returns ``('', 0)``
            if *smiles* cannot be parsed.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "", 0

        # Count unsaturated bonds on the original (no explicit H) molecule.
        unsaturated_types = (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE)
        n_unsaturated = sum(
            1 for bond in mol.GetBonds() if bond.GetBondType() in unsaturated_types
        )
        moles_of_h2 = n_unsaturated * 2

        # Saturate on a copy with explicit hydrogens.
        mol_h = Chem.AddHs(mol)
        emol = Chem.EditableMol(mol_h)

        for bond in mol_h.GetBonds():
            btype = bond.GetBondType()
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()

            if btype == Chem.BondType.DOUBLE:
                emol.RemoveBond(a1, a2)
                emol.AddBond(a1, a2, Chem.BondType.SINGLE)
                self._add_hydrogens(emol, a1, 1)
                self._add_hydrogens(emol, a2, 1)
            elif btype == Chem.BondType.TRIPLE:
                emol.RemoveBond(a1, a2)
                emol.AddBond(a1, a2, Chem.BondType.SINGLE)
                self._add_hydrogens(emol, a1, 2)
                self._add_hydrogens(emol, a2, 2)

        saturated = emol.GetMol()
        Chem.SanitizeMol(saturated)
        return Chem.MolToSmiles(saturated), moles_of_h2

    # ------------------------------------------------------------------
    # Prediction orchestration
    # ------------------------------------------------------------------

    def _run_predictions(self) -> None:
        """Run all enabled property models and trim ``self.df_mol`` to outputs.

        Prediction column names are tracked explicitly so that the output
        DataFrame is always correct even if a property is disabled or its model
        file is missing.  On model error the exception propagates to the caller.
        """
        output_cols: list = ["smiles"]
        n = len(self.df_mol) if self.df_mol is not None else 0

        self._status(f"All descriptors ready — running models on {n} molecule(s)…")

        if self.props.get("MP"):
            self._status("Predicting melting point (MP)…")
            self.df_mol["MP_pred_K"] = self._predict_MP()
            output_cols.append("MP_pred_K")

        if self.props.get("BP"):
            self._status("Predicting boiling point (BP)…")
            self.df_mol["BP_pred_K"] = self._predict_BP()
            output_cols.append("BP_pred_K")

        if self.props.get("HC"):
            self._status("Predicting heat capacity (HC)…")
            self.df_mol["HC_pred_JK"] = self._predict_HC()
            output_cols.append("HC_pred_JK")

        if self.props.get("dH"):
            self._status("Predicting heat of hydrogenation (dH)…")
            self.df_mol["dH_pred_JK"] = self._predict_HH()
            output_cols.append("dH_pred_JK")

        # Drop all raw descriptor columns — keep only SMILES + predictions.
        self.df_mol = self.df_mol[output_cols].copy()

        # Discard rows with physically unreasonable melting point predictions.
        if "MP_pred_K" in self.df_mol.columns:
            mask = self.df_mol["MP_pred_K"] > 6000
            if mask.any():
                for smi in self.df_mol.loc[mask, "smiles"]:
                    mp_val = self.df_mol.loc[self.df_mol["smiles"] == smi, "MP_pred_K"].iloc[0]
                    logger.warning(
                        "Discarded outlier: '%s'  MP_pred=%.1f K (> 6000 K)", smi, mp_val
                    )
                n_discarded = mask.sum()
                self._status(
                    f"Discarded {n_discarded} compound(s) with MP > 6000 K — see p2mat.log."
                )
                self.df_mol = self.df_mol[~mask].reset_index(drop=True)

        if self.props.get("h2_uptake") and not self.df_mol.empty:
            pairs = [self._saturate_smiles(s) for s in self.df_mol["smiles"]]
            sat_smiles, h2_vals = zip(*pairs) if pairs else ([], [])
            self.df_mol["Saturated_smiles"] = list(sat_smiles)
            self.df_mol["H2_uptake"] = list(h2_vals)

            if self._depict is not None:
                for s in self.df_mol["smiles"]:
                    self._depict.save_svg_to_png(s)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def get_data(self) -> pd.DataFrame | None:
        """Return the prediction results DataFrame.

        Returns:
            A ``DataFrame`` with a ``smiles`` column plus one column per
            enabled property, or ``None`` if no predictions have been run.
        """
        return self.df_mol

    def prediction_from_smile(self, smiles: list) -> list:
        """Validate *smiles*, compute descriptors, and run enabled models.

        Args:
            smiles: A list of SMILES strings (may include invalid ones).

        Returns:
            A list of SMILES strings that failed RDKit validation (may be
            empty).  Results are stored in ``self.df_mol``.
        """
        self._status(f"Validating {len(smiles)} SMILES…")
        valid, invalid = self._smile.get_valid_smiles(smiles)
        self._status(
            f"{len(valid)} valid, {len(invalid)} invalid SMILES found."
        )
        if self._on_valid_count is not None:
            self._on_valid_count(len(valid))
        if valid:
            self._compute_features(valid)
            self._run_predictions()
            invalid = invalid + self._padel_failed
        return invalid

    def prediction_from_pid(self, pubchem_ids) -> list:
        """Predict properties from PubChem compound IDs.

        Args:
            pubchem_ids: A single ``int`` or a list of ``int`` PubChem CIDs.

        Returns:
            A list of SMILES strings that failed validation (should be empty
            for valid PubChem IDs).
        """
        if isinstance(pubchem_ids, int):
            pubchem_ids = [pubchem_ids]
        # get_valid_smiles_from_pid returns (valid, invalid) — unpack correctly.
        valid, invalid = self._smile.get_valid_smiles_from_pid(pubchem_ids)
        if invalid and self._verbose > 0:
            print(f"Warning: {len(invalid)} SMILES from PubChem failed validation.")
        if valid:
            self._compute_features(valid)
            self._run_predictions()
        return invalid
