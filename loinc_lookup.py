#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import time
import math
import re
import argparse
from typing import List, Dict, Any, Optional, Tuple

import requests
from requests.auth import HTTPBasicAuth

try:
    import pandas as pd
except ImportError:
    pd = None


FHIR_BASE = "https://fhir.loinc.org"
# LOINC's implicit ValueSet for all codes (NOT the older /vs/loinc)
VALUESET_ALL_LOINC = "http://loinc.org/vs"
DEFAULT_SLEEP_SEC = 0.2            # politeness throttle between $lookup calls

# >>> Recall & breadth controls (raised caps; still bounded for runtime sanity)
DEFAULT_MAX_VARIANT_RESULTS = 200  # per-variant trim cap after lightweight scoring (was 20)
GLOBAL_CANDIDATE_CAP = 1200        # guardrail across variants (was 400)