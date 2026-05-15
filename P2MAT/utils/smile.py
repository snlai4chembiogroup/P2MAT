######################################################
# Developer : Methun Kamruzzaman
# Date      : August, 2024
# Purpose   : SMILES string validation using RDKit and canonical SMILES
#             retrieval from PubChem by compound ID.
######################################################
"""
utils/smile.py
--------------
SMILES string validation and retrieval utilities.

Provides:
- Structural validity checking via RDKit (``Chem.MolFromSmiles``).
- Aromaticity / 6-membered ring detection.
- PubChem compound ID → canonical SMILES lookup via PubChemPy.
"""
import warnings

import rdkit
from pubchempy import Compound
from rdkit import Chem

warnings.filterwarnings("ignore")


class SMILE:
    """Validates SMILES strings and optionally retrieves them from PubChem."""

    def __init__(self, verbose: int = 0) -> None:
        """
        Args:
            verbose: Logging verbosity level.  0 = silent, >0 = print errors.
        """
        self.verbose = verbose

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _check_validity(self, smiles_list: list) -> dict:
        """Check whether each SMILES in *smiles_list* is structurally valid.

        RDKit's ``Chem.MolFromSmiles`` returns ``None`` for an unparseable
        structure and may also raise in rare edge cases.  Both outcomes are
        treated as invalid.

        Args:
            smiles_list: A list of SMILES strings to validate.

        Returns:
            A dict mapping each input SMILES string to ``True`` (valid) or
            ``False`` (invalid).
        """
        results: dict = {}
        for s in smiles_list:
            rv = None  # initialised here so the finally block always has a binding
            try:
                rv = Chem.MolFromSmiles(s)
            except Exception as e:
                if self.verbose > 0:
                    print(f"RDKit error for '{s}': {e}")
            finally:
                results[s] = isinstance(rv, rdkit.Chem.rdchem.Mol)
        return results

    def _fetch_smile(self, pubchem_id: int) -> str:
        """Fetch the isomeric SMILES for a PubChem compound by its CID.

        Args:
            pubchem_id: A valid PubChem compound identifier (integer).

        Returns:
            The isomeric SMILES string for the compound.
        """
        return Compound.from_cid(pubchem_id).isomeric_smiles

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def is_aromatic(self, smile: str) -> bool:
        """Return ``True`` if the molecule contains at least one 6-membered ring.

        Note: this is a structural heuristic (ring size == 6) rather than a
        strict aromaticity check.  It uses RDKit ring-info rather than atom
        aromaticity flags so that non-aromatic 6-rings (e.g. cyclohexane) also
        return ``True``.

        Args:
            smile: A single SMILES string.

        Returns:
            ``True`` if a 6-membered ring is present, ``False`` otherwise or if
            the SMILES is invalid.
        """
        validity = self._check_validity([smile])  # wrap in list — fixes original bug
        if not validity.get(smile, False):
            return False
        mol = Chem.MolFromSmiles(smile)
        ring_info = mol.GetRingInfo()
        return any(
            ring_info.IsAtomInRingOfSize(atom.GetIdx(), 6)
            for atom in mol.GetAtoms()
        )

    def get_valid_smiles(self, smiles_list: list) -> tuple:
        """Split a list of SMILES into valid and invalid sublists.

        Args:
            smiles_list: Candidate SMILES strings.

        Returns:
            A 2-tuple ``(valid, invalid)`` where each element is a list of
            SMILES strings.
        """
        validity = self._check_validity(smiles_list)
        valid = [s for s, ok in validity.items() if ok]
        invalid = [s for s, ok in validity.items() if not ok]
        return valid, invalid

    def get_valid_smiles_from_pid(self, pubchem_id_list: list) -> tuple:
        """Fetch SMILES from PubChem by CID, then validate and split.

        Args:
            pubchem_id_list: A list of integer PubChem compound IDs.

        Returns:
            A 2-tuple ``(valid, invalid)`` — same contract as
            :meth:`get_valid_smiles`.
        """
        smiles = [self._fetch_smile(pid) for pid in pubchem_id_list]
        return self.get_valid_smiles(smiles)
