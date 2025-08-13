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

# ----------------------------
# Credentials & HTTP
# ----------------------------


def _tidy(s: str) -> str:
    """Strip whitespace and newlines.

    Parameters
    ----------
    s : str
        Input string (possibly with surrounding whitespace).

    Returns
    -------
    str
        Tidied string with trailing/leading whitespace and line breaks removed.
    """
    return (s or "").strip().strip("\r").strip("\n")


def read_creds_from_file(path: str) -> Optional[Tuple[str, str]]:
    """Read credentials from a two-line file.

    Parameters
    ----------
    path : str
        File path. Line 1 = username (email), line 2 = password.

    Returns
    -------
    tuple of (str, str) or None
        Username and password if file exists and is valid, else ``None``.

    Raises
    ------
    RuntimeError
        If the file is present but missing required lines.
    """
    if not path or not os.path.isfile(path):
        return None
    with open(path, "r", encoding="utf-8") as f:
        lines = [ln.strip() for ln in f.readlines()]
    lines = [ln for ln in lines if ln.strip()]
    if len(lines) < 2:
        raise RuntimeError(
            f"Credentials file '{path}' must have two non-empty lines:\n"
            "  line 1 = username (email)\n  line 2 = password"
        )
    return _tidy(lines[0]), _tidy(lines[1])


def resolve_credentials(creds_path: Optional[str]) -> Tuple[str, str, str]:
    """Resolve LOINC credentials from args, env, or local file.

    Resolution order
    ----------------
    1. ``--creds`` file (two lines)
    2. Environment variables ``LOINC_USER`` / ``LOINC_PASS``
    3. Local file ``./loinc_creds.txt`` (two lines)

    Parameters
    ----------
    creds_path : str or None
        Path to credential file, or ``None``.

    Returns
    -------
    tuple of (str, str, str)
        Username, password, and source descriptor.

    Raises
    ------
    SystemExit
        If no credentials can be found.
    """
    if creds_path:
        got = read_creds_from_file(creds_path)
        if got:
            u, p = got
            if not u or not p:
                raise SystemExit(f"Credentials in '{creds_path}' are missing/empty.")
            return u, p, f"file:{creds_path}"

    env_u = _tidy(os.environ.get("LOINC_USER", ""))
    env_p = _tidy(os.environ.get("LOINC_PASS", ""))
    if env_u and env_p:
        return env_u, env_p, "env:LOINC_USER/LOINC_PASS"

    got = read_creds_from_file("loinc_creds.txt")
    if got:
        u, p = got
        if not u or not p:
            raise SystemExit("Found 'loinc_creds.txt' but a line was empty.")
        return u, p, "file:loinc_creds.txt"

    raise SystemExit(
        "No LOINC credentials found.\n"
        "Provide --creds /path/to/loinc_creds.txt (two lines), or set env LOINC_USER/LOINC_PASS,\n"
        "or put 'loinc_creds.txt' in the working directory."
    )


def make_requests_session(user: str, pw: str) -> requests.Session:
    """Create an authenticated `requests.Session` for LOINC FHIR.

    Parameters
    ----------
    user : str
        LOINC username (email).
    pw : str
        LOINC password.

    Returns
    -------
    requests.Session
        Session with HTTP Basic auth and JSON/FHIR headers set.
    """
    s = requests.Session()
    s.auth = HTTPBasicAuth(user, pw)
    s.headers.update(
        {
            "Accept": "application/fhir+json",
            "User-Agent": "loinc-lookup/2.1 (+local usage)",
        }
    )
    return s


def http_get_json(
    session: requests.Session,
    url: str,
    params: Optional[Dict[str, Any]] = None,
    max_retries: int = 3,
    sleep_base: float = 0.5,
) -> Dict[str, Any]:
    """GET JSON with basic retry/backoff for common transient errors.

    Behavior
    --------
    - 401 → raise ``RuntimeError`` with a curl hint for debugging credentials.
    - 429/502/503/504 → exponential backoff using ``sleep_base``.
    - Other non-OK → ``raise_for_status``.
    - OK → return parsed JSON payload.

    Parameters
    ----------
    session : requests.Session
        Authenticated session.
    url : str
        Endpoint URL.
    params : dict, optional
        Query parameters.
    max_retries : int, default 3
        Maximum retry attempts for transient errors.
    sleep_base : float, default 0.5
        Base backoff (seconds).

    Returns
    -------
    dict
        Parsed JSON payload.

    Raises
    ------
    RuntimeError
        On 401 Unauthorized with guidance.
    requests.HTTPError
        For non-success status codes after retries.
    """
    attempt = 0
    while True:
        resp = session.get(url, params=params, timeout=60)
        if resp.status_code == 401:
            raise RuntimeError(
                "401 Unauthorized from LOINC FHIR.\n"
                "Check username/password; quick test:\n"
                "  curl -u '<user>' 'https://fhir.loinc.org/CodeSystem?url=http://loinc.org&_format=json'"
            )
        if resp.status_code in (429, 502, 503, 504):
            attempt += 1
            if attempt > max_retries:
                resp.raise_for_status()
            retry_after = resp.headers.get("Retry-After")
            if retry_after:
                try:
                    sleep_s = float(retry_after)
                except ValueError:
                    sleep_s = sleep_base * (2 ** (attempt - 1))
            else:
                sleep_s = sleep_base * (2 ** (attempt - 1))
            time.sleep(sleep_s)
            continue
        resp.raise_for_status()
        return resp.json()


def assert_auth_ok(session: requests.Session) -> None:
    """Quick probe to verify credentials are valid.

    This function performs a small authenticated request against the LOINC FHIR
    endpoint and **only** raises an exception on failure. It does not return a
    value and is intended purely for side-effect validation.

    Parameters
    ----------
    session : requests.Session
        Authenticated LOINC FHIR session.

    Raises
    ------
    requests.HTTPError
        If the probe endpoint does not return success.
    RuntimeError
        If the server reports ``401 Unauthorized``.
    """
    # Intentionally ignore the response object; we only care that no exception is raised.
    http_get_json(
        session,
        f"{FHIR_BASE}/CodeSystem",
        params={"url": "http://loinc.org", "_format": "json"},
    )
