######################################################
# Developer : Methun Kamruzzaman
# Date      : August, 2024
# Purpose   : Utility functions for path resolution, CSS loading, and
#             configuration management. Supports both development and
#             PyInstaller-frozen execution environments.
######################################################
"""
utils/helper.py
---------------
Utility functions for path resolution, CSS loading, and configuration management.

Supports both development mode and PyInstaller-frozen executables by detecting
the presence of ``sys._MEIPASS`` at runtime.
"""
import json
import os
import sys


def get_pkg_dir() -> str:
    """Return the absolute base directory of the package.

    When the application is frozen by PyInstaller, ``sys._MEIPASS`` points to
    the temporary directory where bundled resources are extracted.  In normal
    development mode the current working directory is used instead.
    """
    if getattr(sys, "frozen", False):
        return sys._MEIPASS
    return os.path.abspath(".")


def resource_path(relative_path: str) -> str:
    """Resolve *relative_path* against the package base directory.

    Use this for all asset lookups (models, feature lists, images, CSS) so
    that paths work correctly in both development and frozen-binary contexts.
    """
    return os.path.join(get_pkg_dir(), relative_path)


def read_and_replace_css(css_path: str) -> str:
    """Read a CSS file and return its contents as a string.

    Args:
        css_path: Absolute path to the CSS file.

    Returns:
        The raw CSS text.
    """
    with open(css_path, "r", encoding="utf-8") as fh:
        return fh.read()


def load_config() -> dict:
    """Load ``config.json`` from the package root and return it as a dict.

    If the file is missing or cannot be parsed, a hard-coded fallback is
    returned so the application can always start without crashing.

    Returns:
        A dict with at least a ``"properties"`` key mapping property names to
        booleans, e.g. ``{"properties": {"MP": True, "BP": False, ...}}``.
    """
    _FALLBACK: dict = {
        "properties": {
            "MP": True,
            "BP": False,
            "HC": False,
            "dH": False,
            "h2_uptake": False,
        }
    }
    path = resource_path("config.json")
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:
        return _FALLBACK
