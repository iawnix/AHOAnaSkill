from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parents[1]
DEFAULT_DB = SKILL_DIR / "data" / "aho.sqlite"
DEFAULT_REPORTS_DIR = SKILL_DIR / "reports"
DEFAULT_FIGURES_DIR = DEFAULT_REPORTS_DIR / "figures"
DEFAULT_CUSTOM_SCRIPTS_DIR = DEFAULT_REPORTS_DIR / "scripts"
SDF_INDEX_PATH = SKILL_DIR / "data" / "sdf_index.json"

CSV_COLUMNS = [
    "DATA_ID",
    "CAT_NAME",
    "SOL_NAME",
    "PRO_R_NAME",
    "PRO_S_NAME",
    "REA_NAME",
    "TEMP",
    "PRESSURE",
    "EE",
]

ROLE_TO_COLUMN = {
    "REA": "REA_NAME",
    "CAT": "CAT_NAME",
    "SOL": "SOL_NAME",
    "PRO_R": "PRO_R_NAME",
    "PRO_S": "PRO_S_NAME",
}

ROLES = tuple(ROLE_TO_COLUMN.keys())

FUNCTIONAL_GROUP_SMARTS = {
    "alcohol": "[OX2H][CX4]",
    "phenol": "[OX2H]c",
    "ether": "[OD2]([#6])[#6]",
    "methoxy_aromatic": "c[OX2][CH3]",
    "aldehyde": "[CX3H1](=O)[#6]",
    "ketone": "[#6][CX3](=O)[#6]",
    "carboxylic_acid": "[CX3](=O)[OX2H1]",
    "ester": "[#6][CX3](=O)[OX2H0][#6]",
    "amide": "[NX3][CX3](=[OX1])",
    "amine_primary": "[NX3;H2][#6]",
    "amine_secondary": "[NX3;H1]([#6])[#6]",
    "amine_tertiary": "[NX3;H0]([#6])([#6])[#6]",
    "aniline": "c[NX3]",
    "nitrile": "[CX2]#N",
    "nitro": "[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
    "halide": "[F,Cl,Br,I]",
    "fluoro": "[F]",
    "chloro": "[Cl]",
    "bromo": "[Br]",
    "iodo": "[I]",
    "alkene": "C=C",
    "alkyne": "C#C",
    "aromatic_ring_atom": "a",
    "benzene_ring": "c1ccccc1",
    "heteroaromatic_atom": "[a;!c]",
    "pyridine_like_n": "n",
    "thiol": "[SX2H]",
    "thioether": "[SX2]([#6])[#6]",
    "sulfone": "S(=O)(=O)",
    "sulfoxide": "S=O",
    "phosphine": "[PX3]([#6])([#6])[#6]",
    "phosphine_oxide": "[PX4](=O)",
    "boronic_acid": "B(O)O",
    "silane": "[Si]",
    "isocyanate": "N=C=O",
    "urea": "NC(=O)N",
    "carbonate": "OC(=O)O",
    "carbamate": "NC(=O)O",
    "imine": "C=N",
    "enone": "C=CC=O",
}

SMARTS_ALIASES = {
    "c-OMe": "c[OX2][CH3]",
    "aryl-OMe": "c[OX2][CH3]",
    "phenol": FUNCTIONAL_GROUP_SMARTS["phenol"],
    "benzene": FUNCTIONAL_GROUP_SMARTS["benzene_ring"],
}

METAL_SYMBOLS = {
    "Li",
    "Na",
    "K",
    "Mg",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Ir",
    "Pt",
    "Au",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def resolve_smarts(pattern: str) -> str:
    return SMARTS_ALIASES.get(pattern, FUNCTIONAL_GROUP_SMARTS.get(pattern, pattern))
