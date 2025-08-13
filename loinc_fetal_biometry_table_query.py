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

