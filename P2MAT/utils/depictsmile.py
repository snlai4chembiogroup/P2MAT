######################################################
# Developer : Methun Kamruzzaman
# Date      : August, 2024
# Purpose   : Render 2D molecular structure images from SMILES strings
#             (SVG via RDKit, rasterised to PNG via cairosvg).
######################################################
"""
utils/depictsmile.py
--------------------
Molecule image rendering utilities.

Converts a SMILES string to a 2D structural image (SVG → PNG) using RDKit
for drawing and cairosvg for rasterisation.  Optionally overlays the SMILES
text at the bottom of the image in wrapped lines.
"""
import base64
import hashlib
import os

from cairosvg import svg2png
from PIL import Image, ImageDraw
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D


class DepictSmile:
    """Renders a SMILES string as a PNG file and optionally annotates it.

    Images are saved to *save_file_path* with filenames derived from the
    SHA-1 hash of the SMILES string, so the same molecule always maps to the
    same file.

    Example::

        renderer = DepictSmile(save_file_path="./images")
        renderer.save_svg_to_png("c1ccccc1")
    """

    def __init__(
        self,
        save_file_path: str = "./",
        size: tuple = (300, 300),
        dpi: int = 200,
        verbose: int = 0,
    ) -> None:
        """
        Args:
            save_file_path: Directory where PNG files will be written.
            size: ``(width, height)`` of the output image in pixels.
            dpi: Resolution used when converting SVG to PNG.
            verbose: Logging verbosity.  0 = silent.
        """
        self.mol_size = size
        self.file_path = save_file_path
        self.dpi = dpi
        self.verbose = verbose

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _mol_to_svg(self, mol, kekulize: bool = True) -> str:
        """Render an RDKit molecule object to an SVG string.

        Args:
            mol: An RDKit ``Mol`` object.
            kekulize: If ``True``, attempt Kekulé form before drawing.

        Returns:
            An SVG string with the ``svg:`` namespace prefix stripped.
        """
        mc = Chem.Mol(mol.ToBinary())
        if kekulize:
            try:
                Chem.Kekulize(mc)
            except Exception:
                mc = Chem.Mol(mol.ToBinary())

        if not mc.GetNumConformers():
            rdDepictor.Compute2DCoords(mc)

        drawer = rdMolDraw2D.MolDraw2DSVG(self.mol_size[0], self.mol_size[1])
        drawer.DrawMolecule(mc)
        drawer.FinishDrawing()
        return drawer.GetDrawingText().replace("svg:", "")

    def _sha1(self, text: str) -> str:
        """Return the SHA-1 hex digest of *text* (used as a stable filename)."""
        return hashlib.sha1(text.encode()).hexdigest()

    def _get_file_name(self, smiles: str) -> str:
        """Derive a deterministic filename (no extension) from the SMILES string."""
        return self._sha1(smiles)

    def _add_text_to_image(self, img_path: str, title: str) -> None:
        """Overlay *title* text at the bottom of the image at *img_path*.

        The title is split into lines of at most 40 characters and stacked
        upward from the bottom of the image with a small left margin.

        Args:
            img_path: Absolute path to an existing PNG file (modified in place).
            title: The text string to overlay (typically the SMILES string).
        """
        img = Image.open(img_path)
        draw = ImageDraw.Draw(img)
        img_w, img_h = img.size  # PIL: (width, height)

        line_height: int = 14   # pixels between consecutive text lines
        margin_left: int = 4    # horizontal offset from the left edge
        margin_bottom: int = 4  # gap between the last line and the image bottom

        # Wrap the title into 40-character chunks
        chunks = [title[i : i + 40] for i in range(0, len(title), 40)]

        for i, chunk in enumerate(chunks):
            # Stack lines from the bottom upward
            x = margin_left
            y = img_h - margin_bottom - (len(chunks) - i) * line_height
            draw.text((x, y), chunk, fill=(0, 0, 0))

        img.save(img_path)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def save_svg_to_png(self, smiles: str) -> None:
        """Render *smiles* to a PNG file and annotate it with the SMILES text.

        The output filename is ``<sha1(smiles)>.png`` inside ``self.file_path``.

        Args:
            smiles: A valid SMILES string.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            if self.verbose > 0:
                print(f"Invalid SMILES, cannot render: {smiles}")
            return

        svg = self._mol_to_svg(mol)
        filename = self._get_file_name(smiles)
        filepath = os.path.join(self.file_path, f"{filename}.png")

        svg2png(bytestring=svg, dpi=self.dpi, write_to=filepath)
        self._add_text_to_image(filepath, smiles)

        if self.verbose > 0:
            print(f"Image saved to {filepath}")

    def get_decoded_image(self, smiles: str) -> str:
        """Return a Base64-encoded PNG for *smiles*, rendering it if needed.

        Args:
            smiles: A valid SMILES string.

        Returns:
            A UTF-8 Base64 string of the PNG image bytes.
        """
        filename = self._get_file_name(smiles)
        img_path = os.path.join(self.file_path, f"{filename}.png")

        if not os.path.exists(img_path):
            self.save_svg_to_png(smiles)

        with open(img_path, "rb") as fh:
            return base64.b64encode(fh.read()).decode("utf-8")
