"""Configuration and runtime path handling for the TES tank model.

This module isolates filesystem/config side effects so that importing the model
module does not create directories or read JSON.

Legacy globals are intentionally preserved for backward compatibility.
"""

from __future__ import annotations

import json
import os
import logging
from pathlib import Path
from math import pi

logger = logging.getLogger(__name__)

# Base paths
BASE_DIR = Path(__file__).resolve().parent
SCRIPT_DIR = str(BASE_DIR)
CONFIG_PATH = BASE_DIR / "config.json"

# ---------------------------------------------------------------------------
# Legacy-compatible path globals (populated lazily)
# ---------------------------------------------------------------------------

# Legacy relative strings (kept; the original code mixes / and \\)
OUTPUTS_PFAD = "datei/"
FOLDER_TEMPERATURPROFILE = "datei/temperaturprofile/"
SZ_DAT = "sz/"

folder_temperaturProfile = FOLDER_TEMPERATURPROFILE
outputs_pfad = OUTPUTS_PFAD
sz_dat = SZ_DAT

abs_folder_temperaturProfile: str = os.path.join(SCRIPT_DIR, folder_temperaturProfile)
abs_outputs_pfad: str = os.path.join(SCRIPT_DIR, outputs_pfad)
sz_folder: str = os.path.join(abs_outputs_pfad, sz_dat)
config_path: str = str(CONFIG_PATH)

# Cache for storage parameters loaded from config.json
speicher_param: dict = {}


def ensure_output_dirs() -> None:
    """Create output folders if needed (lazy, called by `ensure_initialized`)."""
    global abs_folder_temperaturProfile, abs_outputs_pfad, sz_folder

    abs_folder_temperaturProfile = os.path.join(SCRIPT_DIR, FOLDER_TEMPERATURPROFILE)
    os.makedirs(abs_folder_temperaturProfile, exist_ok=True)

    abs_outputs_pfad = os.path.join(SCRIPT_DIR, OUTPUTS_PFAD)
    # Keep legacy behavior: code sometimes uses relative folder `datei/` too.
    os.makedirs(OUTPUTS_PFAD, exist_ok=True)

    sz_folder = os.path.join(abs_outputs_pfad, SZ_DAT)
    os.makedirs(sz_folder, exist_ok=True)


def load_speicher_param(*, force_reload: bool = False) -> dict:
    """Load and cache storage parameters from config.json."""
    global speicher_param
    if speicher_param and not force_reload:
        return speicher_param

    cfg_path = str(CONFIG_PATH)
    with open(cfg_path, encoding="utf-8") as f:
        data = json.load(f)

    speicher_param = data["SPEICHER_PARAMETER"]
    speicher_param["A_Quer"] = pi * speicher_param["R_innen"] ** 2
    return speicher_param


def ensure_initialized(*, force_config_reload: bool = False, **_ignored) -> dict:
    """Ensure folders exist and parameters are loaded; returns `speicher_param`.

    The legacy model sometimes passed extra arguments (e.g., ``T_DR``) into an
    initialization helper. These are ignored here to maintain compatibility.
    """
    ensure_output_dirs()
    return load_speicher_param(force_reload=force_config_reload)
