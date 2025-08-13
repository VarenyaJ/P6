#!/usr/bin/env python3

from __future__ import annotations

import os
import time
import math
import re
import argparse
import logging
from typing import List, Dict, Any, Optional, Tuple, Set

import requests
from requests.auth import HTTPBasicAuth

try:
    import pandas as pd
except ImportError:
    pd = None

# ----------------------------
# FHIR server configuration
# ----------------------------

FHIR_BASE = "https://fhir.loinc.org"
# LOINC's implicit ValueSet for all codes (NOT the older /vs/loinc)
VALUESET_ALL_LOINC = "http://loinc.org/vs"
DEFAULT_SLEEP_SEC = 0.2  # politeness throttle between $lookup calls

# Recall & breadth controls
DEFAULT_MAX_VARIANT_RESULTS = 200  # per-variant trim cap after lightweight scoring
GLOBAL_CANDIDATE_CAP = 1200  # guardrail across variants

# ----------------------------
# Synonyms and expectations
# ----------------------------

# Basic clinical abbreviations → plain words for search expansion
ABBREV_TO_WORDS = {
    "BPD": "biparietal diameter",
    "HC": "head circumference",
    "AC": "abdominal circumference",
    "FL": "femur length",
    "EFW": "estimated fetal weight",
    "FHR": "fetal heart rate",
}

# Soft-preference tokens (we prefer but never hard-restrict)
SOFT_CONTEXT_TOKENS = ["fetal", "fetus", "ultrasound", "us", "pregnancy", "obstetric"]

# Hand-curated hints for tricky displays that filter text alone often misses.
HARD_HINTS: Dict[str, List[str]] = {
    "hc/ac": [
        "Head circumference/Abdominal circumference derived by US",
        "Head circumference/Abdominal circumference derived by ultrasound",
    ],
    "fl/ac": [
        "Femur length/Abdominal circumference derived by US",
        "Femur length/Abdominal circumference derived by ultrasound",
    ],
    "fl/bpd": [
        "Femur length/Biparietal diameter derived by US",
        "Femur length/Biparietal diameter derived by ultrasound",
    ],
    "efw": [
        "Estimated fetal weight [Mass] US",
        "Fetal body weight [Mass] US",
        "Fetus body weight [Mass] US",
        "EFW [Mass] US",
    ],
    "placenta appearance": [
        "Placenta morphology [Text] US",
        "Placenta structure [Presence] fetus US",
        "Placenta abnormality [Presence] fetus US",
    ],
    "heart abnormal": [
        "Cardiac abnormality [Presence] fetus US",
        "Heart anomaly [Presence] fetus US",
    ],
    "head abnormal": [
        "Head abnormality [Presence] fetus US",
        "Cranial abnormality [Presence] fetus US",
    ],
    "face/neck abnormal": [
        "Face abnormality [Presence] fetus US",
        "Neck abnormality [Presence] fetus US",
        "Facial anomaly [Presence] fetus",
    ],
    "spine abnormal": [
        "Spinal abnormality [Presence] fetus US",
        "Spine anomaly [Presence] fetus",
    ],
    "genitalia normal": [
        "Genitalia normal [Presence] fetus",
        "Genitalia abnormality [Presence] fetus",
    ],
}

# Default term list requested earlier
DEFAULT_TERMS = [
    "Gestational Age",
    "BPD (cm)",
    "HC (cm)",
    "AC (cm)",
    "Femur (cm)",
    "Cerebellum (cm)",
    "Cisterna Magna (mm)",
    "Humerus (cm)",
    "Radius (cm)",
    "Ulna (cm)",
    "Tibia (cm)",
    "Fibula (cm)",
    "HC/AC",
    "FL/AC",
    "FL/BPD",
    "EFW (g)",
    "EFW Percentile",
    "FHR (bpm)",
    "Presentation",
    "Placenta Appearance",
    "Heart Abnormal",
    "Head Abnormal",
    "Face/Neck Abnormal",
    "Spine Abnormal",
    "Genitalia Normal",
]

