#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set

import pandas as pd

LOGGER = logging.getLogger(__name__)

# Columns to return in the output (plus a leading "Header" column for context)
REQUESTED_OUTPUT_COLUMNS: List[str] = [
    "LOINC_NUM",
    "COMPONENT",
    "PROPERTY",
    "TIME_ASPCT",
    "SYSTEM",
    "SCALE_TYP",
    "METHOD_TYP",
    "CLASS",
    "CLASSTYPE",
    "LONG_COMMON_NAME",
    "SHORTNAME",
]

# LOINC classes strongly associated with OB ultrasound
OB_ULTRASOUND_CLASSES: Set[str] = {"OB.US", "PANEL.OB.US"}

# Header-specific PROPERTY hints to reduce noise
PROPERTY_HINTS_BY_HEADER: Dict[str, Set[str]] = {
    "bpd": {"Len"},
    "biparietal diameter": {"Len"},
    "hc": {"Len", "Prctl"},
    "head circumference": {"Len", "Prctl"},
    "ac": {"Len", "Prctl"},
    "abdominal circumference": {"Len", "Prctl"},
    "efw": {"Mass", "Prctl"},
    "estimated fetal weight": {"Mass", "Prctl"},
    "fhr": {"NRat"},
    "fetal heart rate": {"NRat"},
    "femur": {"Len", "LenRto"},
    "humerus": {"Len"},
    "radius": {"Len"},
    "ulna": {"Len"},
    "tibia": {"Len"},
    "fibula": {"Len"},
    "cerebellum": {"Len"},
    "tcd": {"Len"},
    "cisterna magna": {"Len"},
    "nuchal translucency": {"Len"},
    "nuchal fold": {"Len"},
    "ofd": {"Len"},
    "presentation": {"Find", "Type"},
    "placenta": {"Find"},
    "amniotic": {"Vol", "Find"},
    "gestational age": {"Time", "PrThr"},
}

# Helpful aliases for matching
HEADER_ALIASES: Dict[str, Sequence[str]] = {
    "bpd": ["bpd", "biparietal diameter"],
    "hc": ["hc", "head circumference"],
    "ac": ["ac", "abdominal circumference"],
    "efw": ["efw", "estimated fetal weight", "fetal weight"],
    "fhr": ["fhr", "fetal heart rate"],
    "fl": ["fl", "femur length", "femur"],
    "femur": ["femur", "femur length", "fl"],
    "humerus": ["humerus", "humeral length"],
    "radius": ["radius"],
    "ulna": ["ulna"],
    "tibia": ["tibia"],
    "fibula": ["fibula"],
    "cerebellum": ["cerebellum", "transcerebellar diameter", "tcd"],
    "cisterna magna": ["cisterna magna"],
    "nuchal translucency": ["nuchal translucency", "nt"],
    "nuchal fold": ["nuchal fold"],
    "ofd": ["ofd", "occipitofrontal diameter", "occipito-frontal diameter"],
    "presentation": ["presentation", "fetal presentation"],
    "placenta": ["placenta", "placental"],
    "amniotic": ["amniotic fluid", "amniotic"],
    "gestational age": ["gestational age", "ga"],
}

# Terms that indicate an OB/fetal system context when found in SYSTEM
OB_SYSTEM_KEYWORDS: Sequence[str] = [
    "fetus",
    "fetal",
    "placenta",
    "amniotic",
    "uterus",
    "cervix",
    "umbilical",
]


def parse_cli_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Search local LOINC for codes related to fetal biometry headers "
            "(OB ultrasound–focused)."
        )
    )
    parser.add_argument(
        "--loinc-csv",
        required=True,
        help=(
            "Path to LoincTableCore.csv "
            "(e.g., /Users/you/Downloads/Loinc_2.80/LoincTableCore.csv)"
        ),
    )
    parser.add_argument(
        "--excel",
        required=True,
        help="Path to the Excel workbook containing the fetal biometry sheet.",
    )
    parser.add_argument(
        "--sheet",
        default="fetal_biometry",
        help="Sheet name that contains fetal biometry headers (default: fetal_biometry).",
    )
    parser.add_argument(
        "--out",
        default="fetal_biometry_loinc_matches.filtered.csv",
        help="Output CSV path (default: fetal_biometry_loinc_matches.filtered.csv).",
    )
    parser.add_argument(
        "--max-per-header",
        type=int,
        default=50,
        help="Maximum rows to keep per header after scoring (default: 50).",
    )
    parser.add_argument(
        "--no-context-filter",
        dest="no_context_filter",
        action="store_true",
        help="Disable fetal/OB context filter (not recommended).",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Log per-header diagnostics (raw/filtered/kept counts) at DEBUG level.",
    )
    return parser.parse_args()



def configure_logging(debug: bool) -> None:
    """
    Configure the root logger.

    Parameters
    ----------
    debug : bool
        If True, sets level to DEBUG; otherwise INFO.
    """
    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")



def validate_paths(loinc_csv_path: Path, excel_path: Path) -> None:
    """
    Validate input file paths.

    Parameters
    ----------
    loinc_csv_path : Path
        Path to the `LoincTableCore.csv`.
    excel_path : Path
        Path to the Excel workbook.

    Raises
    ------
    FileNotFoundError
        If either path does not exist.
    """
    if not loinc_csv_path.exists():
        raise FileNotFoundError(f"LOINC CSV not found: {loinc_csv_path}")
    if not excel_path.exists():
        raise FileNotFoundError(f"Excel file not found: {excel_path}")




def load_loinc_dataframe(loinc_csv_path: Path) -> pd.DataFrame:
    """
    Load LOINC table with consistent string dtype and required columns present.

    Parameters
    ----------
    loinc_csv_path : Path
        Path to the `LoincTableCore.csv`.

    Returns
    -------
    pandas.DataFrame
        LOINC table with string columns and NaNs filled by empty strings.

    Raises
    ------
    SystemExit
        If expected columns are missing.
    """
    loinc_dataframe = pd.read_csv(loinc_csv_path, dtype=str, low_memory=False).fillna("")
    missing_columns = [c for c in REQUESTED_OUTPUT_COLUMNS if c not in loinc_dataframe.columns]
    if missing_columns:
        raise SystemExit(f"Missing expected columns in LOINC CSV: {missing_columns}")
    return loinc_dataframe


def read_excel_headers(excel_path: Path, sheet_name: str) -> List[str]:
    """
    Read the header row from the given Excel sheet.

    Parameters
    ----------
    excel_path : Path
        Path to the Excel workbook.
    sheet_name : str
        Sheet name containing fetal biometry headers.

    Returns
    -------
    list of str
        Column names used as search headers.

    Raises
    ------
    SystemExit
        If the sheet is not found.
    """
    excel_file = pd.ExcelFile(excel_path)
    if sheet_name not in excel_file.sheet_names:
        raise SystemExit(f"Sheet '{sheet_name}' not found in {excel_path}. "
                         f"Found: {excel_file.sheet_names}")
    # Read only the header row for efficiency
    df = pd.read_excel(excel_path, sheet_name=sheet_name, nrows=0)
    return list(df.columns)



def normalize_header_to_term(header_text: object) -> str:
    """
    Normalize a column header into a search-friendly base term.

    This:
    - Converts to lowercase.
    - Removes units or hints in parentheses (e.g., ``(cm)``).
    - Collapses multiple spaces.

    Parameters
    ----------
    header_text : object
        The header name (usually a string). Non-strings return an empty string.

    Returns
    -------
    str
        The normalized term (possibly empty).
    """
    if not isinstance(header_text, str):
        return ""
    text = re.sub(r"\(.*?\)", "", header_text)  # remove "(...)"
    text = re.sub(r"\s+", " ", text).strip().lower()
    return text


