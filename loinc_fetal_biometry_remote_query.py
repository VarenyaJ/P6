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