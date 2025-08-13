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

# ----------------------------
# Query normalization
# ----------------------------


def strip_units(user_term: str) -> str:
    """Remove trailing unit hints like ``'(cm)'`` from a term, if present.

    Parameters
    ----------
    user_term : str
        Original input term.

    Returns
    -------
    str
        Term without any trailing ``'( ... )'`` unit hint.
    """
    t = user_term.strip()
    if "(" in t and ")" in t and t.index("(") < t.index(")"):
        return t[: t.index("(")].strip()
    return t


def normalize_user_term(user_term: str) -> str:
    """Normalize a user term (expand common abbreviations, ratios).

    Examples
    --------
    >>> normalize_user_term("HC (cm)")
    'head circumference'
    >>> normalize_user_term("FL/AC")
    'femur length ratio abdominal circumference'

    Parameters
    ----------
    user_term : str
        Input term, e.g., ``"HC (cm)"``, ``"FL/AC"``, ``"EFW (g)"``.

    Returns
    -------
    str
        Normalized canonical phrase used to seed search variants.
    """
    base = strip_units(user_term)
    upper_no_space = base.upper().replace(" ", "")
    for abbr, words in ABBREV_TO_WORDS.items():
        if upper_no_space == abbr or base.upper() == abbr:
            return words
    if "/" in base:
        a, b = [p.strip() for p in base.split("/", 1)]
        a = ABBREV_TO_WORDS.get(a.upper(), a)
        b = ABBREV_TO_WORDS.get(b.upper(), b)
        return f"{a} ratio {b}"
    return base


def user_term_implies_numeric(user_term: str, normalized: str) -> bool:
    """Heuristic: does the term imply a numeric (Qn) quantity?

    Signals
    -------
    - Explicit units like ``(cm)``, ``(mm)``, ``(g)``, ``(bpm)``.
    - Words such as *length*, *diameter*, *circumference*, *weight*, *mass*, *rate*.

    Parameters
    ----------
    user_term : str
        Original user-provided term.
    normalized : str
        Normalized representation.

    Returns
    -------
    bool
        ``True`` if numeric is likely; otherwise ``False``.
    """
    t = f"{user_term} {normalized}".lower()
    if any(u in user_term.lower() for u in ["(cm)", "(mm)", "(g)", "(bpm)"]):
        return True
    for k in ["length", "diameter", "circumference", "weight", "mass", "rate"]:
        if k in t:
            return True
    return False


def is_ratio_intent(user_term: str, normalized: str) -> bool:
    """Return True if the term likely represents a ratio (e.g., HC/AC)."""
    return (" ratio " in f"{user_term} {normalized}".lower()) or (
        "/" in strip_units(user_term)
    )


# ----------------------------
# Variant generation
# ----------------------------


def _loincish_variants_for_family(user_term: str, normalized: str) -> List[str]:
    """Family-specific seed phrases that align with LOINC display conventions.

    Parameters
    ----------
    user_term : str
        Original term as typed by the user.
    normalized : str
        Normalized term.

    Returns
    -------
    list of str
        Candidate phrase variants specific to certain measurement families.
    """
    terms: List[str] = []
    ln = normalized.lower()
    lu = user_term.lower()

    if "biparietal diameter" in ln or "bpd" in lu:
        terms += [
            "Head Diameter.biparietal [Length] fetus US",
            "Head biparietal diameter [Length] fetus US",
            "Fetal head biparietal diameter [Length] US",
        ]
    if "head circumference" in ln or "hc" in lu:
        terms += ["Head [Circumference] fetus US"]
    if "abdominal circumference" in ln or "ac" in lu:
        terms += [
            "Abdomen [Circumference] fetus US",
            "Abdominal circumference fetus US",
        ]

    for b in ["femur", "humerus", "radius", "ulna", "tibia", "fibula"]:
        if b in ln or b in lu:
            terms += [
                f"{b.title()} diaphysis fetus [Length] US",
                f"Fetal {b} length [Length] US",
            ]
            break

    if "cerebellum" in ln:
        terms += [
            "Cerebellum fetus [Diameter] US",
            "Cerebellum Diameter transverse fetus US",
        ]
    if "cisterna magna" in ln:
        terms += [
            "Cisterna magna fetus [Diameter] US",
            "Fetal cisterna magna sagittal diameter US",
        ]

    if "estimated fetal weight" in ln or "efw" in lu:
        terms += [
            "Estimated fetal weight [Mass] US",
            "Fetal body weight [Mass] US",
            "Fetus body weight [Mass] US",
            "EFW [Mass] US",
        ]

    if "fetal heart rate" in ln or "fhr" in lu:
        terms += [
            "Fetal heart rate [Rate] US",
            "Fetal heart rate US",
            "Fetal heart rate by auscultation",
            "Fetal heart rate mean 10 minutes",
            "Fetal heart rate reactivity",
        ]
    return terms


def _qualitative_variants(user_term: str) -> List[str]:
    """Qualitative presence/absence variants for common OB structures.

    Parameters
    ----------
    user_term : str
        Original term.

    Returns
    -------
    list of str
        Candidate qualitative variants if applicable, otherwise empty.
    """
    lu = user_term.lower()
    qual_map = {
        "placenta appearance": [
            "Placenta morphology [Text] US",
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
    out: List[str] = []
    for k, variants in qual_map.items():
        if k in lu:
            out.extend(variants)
    return out


def _ratio_variants(user_term: str, normalized: str) -> List[str]:
    """Ratio-specific variants (both literal and normalized forms).

    Parameters
    ----------
    user_term : str
        Raw term (e.g., ``"HC/AC"``).
    normalized : str
        Normalized form (e.g., ``"head circumference ratio abdominal circumference"``).

    Returns
    -------
    list of str
        Phrase variants for ratio expressions; empty list if not ratio intent.
    """
    if not is_ratio_intent(user_term, normalized):
        return []
    ln = normalized.lower()
    a, b = None, None
    if " ratio " in ln:
        parts = ln.split(" ratio ")
        if len(parts) == 2:
            a, b = parts[0], parts[1]
    else:
        raw = [p.strip() for p in strip_units(user_term).split("/")]
        if len(raw) == 2:
            a = ABBREV_TO_WORDS.get(raw[0].upper(), raw[0]).lower()
            b = ABBREV_TO_WORDS.get(raw[1].upper(), raw[1]).lower()
    if not (a and b):
        return []

    raw_lit = strip_units(user_term)
    return [
        f"{a}/{b}",
        f"{a}/{b} ratio",
        f"{a} to {b} ratio",
        f"ratio of {a} to {b}",
        f"{a} divided by {b}",
        f"{a.title()} / {b.title()} derived by US",
        f"{a.title()} / {b.title()} derived by ultrasound",
        raw_lit,
        raw_lit.upper(),
    ]


def _inject_hard_hints(user_term: str) -> List[str]:
    """Inject known-good hard hints for ambiguous terms.

    Parameters
    ----------
    user_term : str
        Original term.

    Returns
    -------
    list of str
        Extra curated phrases to increase recall for difficult cases.
    """
    key = strip_units(user_term).lower()
    key_compact = key.replace(" ", "")
    return HARD_HINTS.get(key, []) + HARD_HINTS.get(key_compact, [])


def build_text_variants(user_term: str, normalized: str) -> List[str]:
    """Build a prioritized list of text variants for server-side $expand.

    Strategy
    --------
    - Base normalized text + ultrasound context
    - Bracketed quantity cues (``[Length]``, ``[Diameter]``, ...)
    - LOINC-ish family phrases (e.g., long bone diaphysis)
    - Qualitative presence/morphology phrasing
    - Ratio variants (literal and spelled out)
    - Curated hard hints

    Parameters
    ----------
    user_term : str
        Original user-provided term.
    normalized : str
        Normalized representation.

    Returns
    -------
    list of str
        Unique variant strings (order preserved).
    """
    terms: List[str] = [normalized, f"{normalized} ultrasound", f"{normalized} US"]
    terms += [
        f"{normalized} [Length]",
        f"{normalized} [Diameter]",
        f"{normalized} [Circumference]",
        f"{normalized} [Mass]",
        f"{normalized} [Rate]",
    ]
    terms += _loincish_variants_for_family(user_term, normalized)
    terms += _qualitative_variants(user_term)
    terms += _ratio_variants(user_term, normalized)
    terms += _inject_hard_hints(user_term)

    # Deduplicate preserving order
    seen: set = set()
    uniq: List[str] = []
    for t in terms:
        if t not in seen:
            uniq.append(t)
            seen.add(t)
    return uniq


# ----------------------------
# PROPERTY expectations
# ----------------------------


def expected_properties_for_intent(
    user_term: str, normalized: str
) -> Optional[Set[str]]:
    """Map intent to expected LOINC PROPERTY codes.

    PROPERTY examples
    -----------------
    ``Len`` (length), ``Diam`` (diameter), ``Circ`` (circumference), ``Mass``,
    ``Rate``, ``Rto`` (ratio).

    Parameters
    ----------
    user_term : str
        Raw term.
    normalized : str
        Normalized term.

    Returns
    -------
    set of str or None
        Expected PROPERTY values, or ``None`` if unconstrained.
    """
    text = f"{user_term} {normalized}".lower()
    if "circumference" in text:
        return {"Circ"}
    if "diameter" in text or "bpd" in text:
        return {"Diam", "Len"}
    if any(
        k in text
        for k in [
            "length",
            "femur",
            "humerus",
            "radius",
            "ulna",
            "tibia",
            "fibula",
            "long bone",
        ]
    ):
        return {"Len"}
    if "heart rate" in text or "fhr" in text or "rate" in text:
        return {"Rate"}
    if any(k in text for k in ["estimated fetal weight", "efw", "weight", "mass"]):
        return {"Mass"}
    if "cisterna magna" in text:
        return {"Diam", "Len"}
    if "cerebellum" in text:
        return {"Diam", "Len", "Circ"}
    if " ratio " in text or "/" in strip_units(user_term):
        return {"Rto"}
    return None


# ----------------------------
# FHIR calls
# ----------------------------


def _fallback_variant(text_filter: str) -> Optional[str]:
    """Loosen a filter by removing bracketed tokens and punctuation.

    Parameters
    ----------
    text_filter : str
        Original filter string.

    Returns
    -------
    str or None
        Fallback filter string if meaningfully different; else ``None``.
    """
    v2 = re.sub(r"\[[^\]]+\]", "", text_filter)
    v2 = re.sub(r"[(),.:;]", " ", v2)
    v2 = re.sub(r"\s+", " ", v2).strip()
    if len(v2) >= 3 and v2 != text_filter:
        return v2
    return None


def expand_valueset_candidates(
    session: requests.Session, text_filter: str, count: int
) -> List[Dict[str, Any]]:
    """Call `$expand` on the all-LOINC ValueSet with a text filter.

    Implements a lightweight fallback if the first call returns nothing.

    Parameters
    ----------
    session : requests.Session
        Authenticated session.
    text_filter : str
        Filter string passed to ``$expand?filter=``.
    count : int
        Server-side count hint.

    Returns
    -------
    list of dict
        Raw ``contains`` entries from the ValueSet expansion (may be empty).
    """
    params = {
        "url": VALUESET_ALL_LOINC,
        "filter": text_filter,
        "count": count,
        "_format": "json",
    }
    data = http_get_json(session, f"{FHIR_BASE}/ValueSet/$expand", params)
    contains = (data.get("expansion", {}) or {}).get("contains", []) or []
    if not contains:
        fb = _fallback_variant(text_filter)
        if fb:
            params["filter"] = fb
            data = http_get_json(session, f"{FHIR_BASE}/ValueSet/$expand", params)
            contains = (data.get("expansion", {}) or {}).get("contains", []) or []
    return contains


def _merge_property_param(out: Dict[str, Any], prop_param: Dict[str, Any]) -> None:
    """Fold a single ``property`` parameter into ``out['properties']``.

    Parameters
    ----------
    out : dict
        Accumulator being built by :func:`_parse_lookup_parameters`. Will gain/extend
        a ``properties`` sub-dict in the shape ``{<LOINC property code>: <value>}``.
    prop_param : dict
        One entry from the FHIR ``parameter`` array whose ``name`` is ``"property"``.
        Its ``part`` array contains the property ``code`` and a single value in one of
        ``valueString`` / ``valueCode`` / ``valueUri``.

    Notes
    -----
    - Missing/partial ``part`` contents are ignored safely.
    - No return value; this mutates ``out`` in place.
    """
    prop_code: Optional[str] = None
    prop_value: Optional[str] = None
    for sub in prop_param.get("part") or []:
        subname = sub.get("name")
        if subname == "code":
            prop_code = sub.get("valueCode")
        elif subname in ("valueString", "valueCode", "valueUri"):
            prop_value = (
                sub.get("valueString") or sub.get("valueCode") or sub.get("valueUri")
            )
    if prop_code and prop_value is not None:
        out.setdefault("properties", {})[prop_code] = prop_value


def _parse_lookup_parameters(data: Dict[str, Any]) -> Dict[str, Any]:
    """Parse the CodeSystem/$lookup ``parameter`` array into a flat dict.

    Extracted fields
    ----------------
    - Top-level strings: ``display``, ``version``, ``status``, ``abstract``, ``definition``
    - Collapsed LOINC properties under ``properties`` as ``{code: value}``

    Parameters
    ----------
    data : dict
        Raw Parameters resource returned by ``$lookup``.

    Returns
    -------
    dict
        Dictionary with flattened top-level keys and a ``properties`` sub-dict.

    Notes
    -----
    - This version is factored to keep cyclomatic complexity low (ruff C901),
      while preserving prior behavior and output shape.
    """
    out: Dict[str, Any] = {}
    keymap = {
        "name": "display",  # Some servers send "name" instead of "display"
        "display": "display",
        "version": "version",
        "status": "status",
        "abstract": "abstract",
        "definition": "definition",
    }

    for p in data.get("parameter") or []:
        name = p.get("name")
        if name == "property":
            _merge_property_param(out, p)
            continue

        target = keymap.get(name)
        if target:
            val = p.get("valueString") or p.get("valueCode") or p.get("valueUri")
            if val is not None:
                out[target] = val

    return out


def lookup_loinc_details(session: requests.Session, code: str) -> Dict[str, Any]:
    """Call `$lookup` for a LOINC code and return a flattened detail dict.

    Parameters
    ----------
    session : requests.Session
        Authenticated session.
    code : str
        LOINC code (e.g., ``"11824-6"``).

    Returns
    -------
    dict
        Flattened detail dictionary including select properties.

    Raises
    ------
    requests.HTTPError
        If the LOINC server returns a non-OK response.
    """
    params = {"system": "http://loinc.org", "code": code, "_format": "json"}
    data = http_get_json(session, f"{FHIR_BASE}/CodeSystem/$lookup", params)
    out: Dict[str, Any] = {"code": code}
    out.update(_parse_lookup_parameters(data))
    return out


# ----------------------------
# Flags & scoring
# ----------------------------


def is_loinc_part(code: Optional[str]) -> bool:
    """True if the code looks like a LOINC Part (``LP...``)."""
    return bool(code) and code.startswith("LP")


def is_loinc_answer_list(code: Optional[str]) -> bool:
    """True if the code looks like a LOINC Answer List/Answer (``LA...``)."""
    return bool(code) and code.startswith("LA")


def is_deprecated(status: Optional[str], display: Optional[str]) -> bool:
    """Detect deprecated records from status or display text."""
    if status and status.upper() == "DEPRECATED":
        return True
    return bool(display) and "deprecated" in display.lower()


def looks_derived_or_methodized(display: Optional[str]) -> bool:
    """Heuristics for methodized/derived content we may want to down-rank.

    Returns
    -------
    bool
        ``True`` if display suggests *derived* or *methodized* content.
    """
    if not display:
        return False
    d = display.lower()
    if "estimated from" in d:
        return True
    if " by " in d and " method" in d:
        return True
    if any(x in d for x in ["hadlock", "jeanty", "merz", "ott", "goldstein"]):
        return True
    if any(x in d for x in ["amniocentesis", "lmp", "menstrual period"]):
        return True
    return False


def mentions_percentile(display: Optional[str]) -> bool:
    """True if display mentions percentile."""
    return bool(display) and "percentile" in display.lower()


def mentions_laterality(display: Optional[str]) -> bool:
    """True if display hints at sidedness (left/right)."""
    if not display:
        return False
    d = display.lower()
    return any(
        w in d for w in [" left ", " right ", " left]", " right]", "left ", "right "]
    )


def _norm_str(x: Optional[str]) -> str:
    """Normalize simple code/property strings for robust comparisons."""
    return (x or "").strip()


def compute_feature_flags(
    details: Dict[str, Any], numeric_expected: bool, expected_props: Optional[Set[str]]
) -> Dict[str, Any]:
    """Compute classification flags used for filtering and scoring.

    Parameters
    ----------
    details : dict
        Flattened output from :func:`lookup_loinc_details`.
    numeric_expected : bool
        Whether the user term implies a numeric (Qn) measurement.
    expected_props : set of str or None
        Acceptable PROPERTY values, or ``None`` if unconstrained.

    Returns
    -------
    dict
        {
          "is_part", "is_answer_list", "is_deprecated", "is_derived",
          "is_percentile", "has_laterality", "property_match", "scale_match"
        }
    """
    disp = details.get("display") or ""
    status = details.get("status")
    props = details.get("properties", {}) if "properties" in details else {}
    code = details.get("code") or ""

    flag_is_part = is_loinc_part(code)
    flag_is_answer = is_loinc_answer_list(code)
    flag_is_depr = is_deprecated(status, disp)
    flag_is_derived = looks_derived_or_methodized(disp)
    flag_is_percentile = mentions_percentile(disp)
    flag_has_laterality = mentions_laterality(disp)

    scale = _norm_str(props.get("SCALE_TYP"))
    prop_kind = _norm_str(props.get("PROPERTY"))
    expected_norm: Set[str] = {(_norm_str(p)) for p in (expected_props or set())}

    prop_match = (not expected_norm) or (prop_kind in expected_norm)
    scale_match = (not numeric_expected) or (scale.upper() == "QN")

    return {
        "is_part": flag_is_part,
        "is_answer_list": flag_is_answer,
        "is_deprecated": flag_is_depr,
        "is_derived": flag_is_derived,
        "is_percentile": flag_is_percentile,
        "has_laterality": flag_has_laterality,
        "property_match": prop_match,
        "scale_match": scale_match,
    }


def prefilter_score_for_variant(display: str, variant_text: str) -> int:
    """Loose pre-score for raw candidates based on text overlap and fetal/US context.

    Parameters
    ----------
    display : str
        Candidate display text returned by ``$expand``.
    variant_text : str
        The filter/variant string used to obtain the candidate.

    Returns
    -------
    int
        Higher is better; used only for lightweight trimming before ``$lookup``.
    """
    disp = (display or "").lower()
    s = 0
    for w in SOFT_CONTEXT_TOKENS:
        if w in disp:
            s += 2
    for tok in variant_text.lower().split():
        if tok and tok in disp:
            s += 1
    return s


def _score_property_scale(
    props: Dict[str, Any], expected_props: Optional[Set[str]], numeric_expected: bool
) -> int:
    """Score component for PROPERTY and SCALE_TYP alignment."""
    s = 0
    if expected_props and props.get("PROPERTY") in expected_props:
        s += 12
    if numeric_expected and props.get("SCALE_TYP") == "Qn":
        s += 6
    if props.get("SCALE_TYP") == "Qn":
        s += 2
    return s


def _score_context(props: Dict[str, Any], disp: str) -> int:
    """Score component for ultrasound/fetal contextual alignment."""
    s = 0
    cls = (props.get("CLASS") or "").lower()
    system = (props.get("SYSTEM") or "").lower()
    sys_core = (props.get("system-core") or "").lower()
    super_sys = (props.get("super-system") or "").lower()

    if "us" in cls or "ob" in cls or "ultrasound" in disp or " us" in disp:
        s += 3
    if (
        "fetal" in disp
        or "fetus" in disp
        or super_sys == "fetus"
        or "fetus" in sys_core
        or "fetus" in system
    ):
        s += 2
    return s


def _score_penalties(
    flags: Dict[str, Any],
    allow_derived: bool,
    allow_percentile: bool,
    want_laterality: bool,
) -> int:
    """Score component applying penalties for derived/percentile/laterality when undesired."""
    s = 0
    if not allow_derived and flags["is_derived"]:
        s -= 10
    elif allow_derived and flags["is_derived"]:
        s -= 4

    if not allow_percentile and flags["is_percentile"]:
        s -= 8
    elif allow_percentile and flags["is_percentile"]:
        s -= 2

    if flags["has_laterality"] and not want_laterality:
        s -= 3
    return s


def _score_misc(
    details: Dict[str, Any],
    disp: str,
    normalized: str,
    numeric_expected: bool,
    props: Dict[str, Any],
) -> int:
    """Score component for narrative/ordinal penalties, token overlap, and safety nudges."""
    s = 0
    if numeric_expected and (
        (props.get("SCALE_TYP") in {"Ord", "Nar"}) or "narrative" in disp
    ):
        s -= 4

    for tok in normalized.lower().split():
        if tok and tok in disp:
            s += 1

    code = details.get("code") or ""
    if not code.startswith(("LP", "LA")):
        s += 1
    if not is_deprecated(details.get("status") or "", details.get("display")):
        s += 1
    return s


def final_score(
    details: Dict[str, Any],
    user_term: str,
    normalized: str,
    flags: Dict[str, Any],
    allow_derived: bool,
    allow_percentile: bool,
    numeric_expected: bool,
    expected_props: Optional[Set[str]],
) -> int:
    """Compute the final ranking score for a candidate.

    Scoring summary
    ---------------
    + PROPERTY/SCALE alignment (dominant if provided)
    + ultrasound/fetal context
    + token overlap with normalized term
    + small nudge against Parts/Answers and deprecated
    - derived/methodized unless allowed
    - percentile unless allowed
    - laterality when not requested
    - narrative/ordinal when numeric expected

    Parameters
    ----------
    details : dict
        Candidate details from :func:`lookup_loinc_details`.
    user_term : str
        Original user term.
    normalized : str
        Normalized term.
    flags : dict
        Feature flags from :func:`compute_feature_flags`.
    allow_derived : bool
        Whether to allow derived/methodized entries.
    allow_percentile : bool
        Whether to allow percentile entries.
    numeric_expected : bool
        Whether a numeric observation is expected.
    expected_props : set of str or None
        Preferred PROPERTY values.

    Returns
    -------
    int
        Higher scores rank earlier.
    """
    disp = (details.get("display") or "").lower()
    props = details.get("properties", {}) if "properties" in details else {}
    want_laterality = any(w in user_term.lower() for w in ["left", "right"])

    s = 0
    s += _score_property_scale(props, expected_props, numeric_expected)
    s += _score_context(props, disp)
    s += _score_penalties(flags, allow_derived, allow_percentile, want_laterality)
    s += _score_misc(details, disp, normalized, numeric_expected, props)
    return s


# ----------------------------
# Filtering & ranking
# ----------------------------


def _sort_key_for_rank(x: Dict[str, Any]) -> Tuple[int, int, int, int, int]:
    """Tie-breaker key to prefer property/scale matches and non-derived/percentile.

    Returns
    -------
    tuple of int
        ``(score, property_match, scale_match, not_derived, not_percentile)`` as ints.
    """
    f = x["_flags"]
    sec = (
        1 if f.get("property_match") else 0,
        1 if f.get("scale_match") else 0,
        0 if f.get("is_derived") else 1,
        0 if f.get("is_percentile") else 1,
    )
    return (x["_score"],) + sec


def filter_and_rank(
    enriched: List[Dict[str, Any]],
    user_term: str,
    normalized: str,
    numeric_expected: bool,
    expected_props: Optional[Set[str]],
    strict: bool,
    allow_derived: bool,
    allow_percentile: bool,
) -> List[Dict[str, Any]]:
    """Apply strict/relaxed filtering and compute scores, returning sorted candidates.

    Strict mode
    -----------
    - Exclude Parts/AnswerLists/Deprecated
    - Exclude derived unless ``allow_derived``
    - Exclude percentile unless ``allow_percentile``
    - If numeric expected: require ``SCALE_TYP=Qn`` and PROPERTY in expected set

    Relaxed mode
    ------------
    - Still excludes Parts/AnswerLists/Deprecated
    - Allows derived/percentile (penalized by score)
    - Drops hard requirements on PROPERTY/SCALE

    Parameters
    ----------
    enriched : list of dict
        Enriched candidates from ``$lookup``.
    user_term : str
        Original user term.
    normalized : str
        Normalized user term.
    numeric_expected : bool
        Whether numeric observation is expected.
    expected_props : set of str or None
        Expected PROPERTY set if any.
    strict : bool
        Whether to apply strict gating rules.
    allow_derived : bool
        Allow derived/methodized entries in strict mode.
    allow_percentile : bool
        Allow percentile entries in strict mode.

    Returns
    -------
    list of dict
        Sorted candidates with ``_flags``, ``_score``, and ``_stage`` annotations.
    """
    accepted: List[Dict[str, Any]] = []
    for d in enriched:
        flags = compute_feature_flags(d, numeric_expected, expected_props)

        if flags["is_part"] or flags["is_answer_list"] or flags["is_deprecated"]:
            continue

        if strict:
            if not allow_derived and flags["is_derived"]:
                continue
            if not allow_percentile and flags["is_percentile"]:
                continue
            if numeric_expected:
                if not flags["scale_match"]:
                    continue
                if expected_props and not flags["property_match"]:
                    continue

        s = final_score(
            details=d,
            user_term=user_term,
            normalized=normalized,
            flags=flags,
            allow_derived=allow_derived,
            allow_percentile=allow_percentile,
            numeric_expected=numeric_expected,
            expected_props=expected_props,
        )
        row = dict(d)
        row["_flags"] = flags
        row["_score"] = s
        row["_stage"] = "strict" if strict else "relaxed"
        accepted.append(row)

    accepted.sort(key=_sort_key_for_rank, reverse=True)
    return accepted


# ----------------------------
# Candidate collection & enrichment
# ----------------------------


def _update_raw_candidates_with_variant(
    raw: Dict[str, Dict[str, Any]], contains: List[Dict[str, Any]], cap: int
) -> None:
    """Merge ``contains`` results into ``raw`` using a simple prefilter score ordering.

    Parameters
    ----------
    raw : dict
        Accumulator mapping code → minimal info.
    contains : list of dict
        Raw results from ``$expand``.
    cap : int
        Per-variant cap after lightweight trimming.
    """
    contains = sorted(
        contains,
        key=lambda c: prefilter_score_for_variant(
            c.get("display", ""), c.get("display", "")
        ),
        reverse=True,
    )
    for c in contains[:cap]:
        code = c.get("code")
        if not code or code in raw:
            continue
        raw[code] = {"code": code, "display": c.get("display")}
        if len(raw) >= GLOBAL_CANDIDATE_CAP:
            break


def gather_raw_candidates(
    session: requests.Session,
    variants: List[str],
    count_per_variant: int,
    per_variant_cap: int,
) -> Dict[str, Dict[str, Any]]:
    """Collect a de-duplicated pool of raw candidate codes via $expand across variants.

    Parameters
    ----------
    session : requests.Session
        Authenticated session.
    variants : list of str
        Text variants to pass to ``$expand``.
    count_per_variant : int
        Server-side candidate count hint per variant.
    per_variant_cap : int
        Local per-variant trim cap after prefilter scoring.

    Returns
    -------
    dict
        Mapping from LOINC code → minimal info from ``$expand``.
    """
    raw_candidates: Dict[str, Dict[str, Any]] = {}
    for v in variants:
        try:
            contains = expand_valueset_candidates(session, v, count=count_per_variant)
        except (requests.RequestException, ValueError) as exc:
            logging.warning(
                "Skipping variant due to request error for '%s': %s", v, exc
            )
            continue
        _update_raw_candidates_with_variant(
            raw_candidates, contains, max(DEFAULT_MAX_VARIANT_RESULTS, per_variant_cap)
        )
        if len(raw_candidates) >= GLOBAL_CANDIDATE_CAP:
            break
    return raw_candidates


def enrich_candidates(
    session: requests.Session,
    raw_candidates: Dict[str, Dict[str, Any]],
    sleep_between_lookups: float,
) -> List[Dict[str, Any]]:
    """Call `$lookup` for each raw candidate and collect enriched detail dicts.

    Parameters
    ----------
    session : requests.Session
        Authenticated session.
    raw_candidates : dict
        Mapping from code → minimal info (from ``$expand``).
    sleep_between_lookups : float
        Delay in seconds between successive ``$lookup`` calls.

    Returns
    -------
    list of dict
        Enriched candidate detail dicts.
    """
    enriched: List[Dict[str, Any]] = []
    for code in raw_candidates:
        time.sleep(sleep_between_lookups)
        try:
            d = lookup_loinc_details(session, code)
            enriched.append(d)
        except (requests.RequestException, ValueError, KeyError) as exc:
            logging.warning("Lookup failed for code %s: %s", code, exc)
            continue
    return enriched


# ----------------------------
# Row builders (CSV-ready)
# ----------------------------


def _empty_best_row(user_term: str, normalized: str, error: str) -> Dict[str, Any]:
    """Build a placeholder 'best row' when no matches were accepted.

    Parameters
    ----------
    user_term : str
        original term
    normalized : str
        normalized term
    error : str
        error description

    Returns
    -------
    dict
        Row with all fields present, ``None`` where unknown, and an ``error`` message.
    """
    return {
        "search_term": user_term,
        "normalized_query": normalized,
        "rank": None,
        "loinc_code": None,
        "display": None,
        "definition": None,
        "status": None,
        "class": None,
        "system": None,
        "property": None,
        "time": None,
        "method": None,
        "scale": None,
        "example_units": None,
        "system_core": None,
        "super_system": None,
        "time_core": None,
        "time_modifier": None,
        "analyte": None,
        "analyte_core": None,
        "analyte_suffix": None,
        "analyte_numerator": None,
        "analyte_divisor": None,
        "analyte_divisor_suffix": None,
        "category": None,
        "search_terms": None,
        "display_name": None,
        "is_part": None,
        "is_answer_list": None,
        "is_deprecated": None,
        "is_derived": None,
        "is_percentile": None,
        "has_laterality": None,
        "property_match": None,
        "scale_match": None,
        "stage": None,
        "score": None,
        "error": error,
    }


def _best_rows_from_ranked(
    chosen: List[Dict[str, Any]], user_term: str, normalized: str
) -> List[Dict[str, Any]]:
    """Convert ranked candidates into compact, CSV-ready rows.

    Parameters
    ----------
    chosen : list of dict
        Ranked candidates selected by strict/relaxed passes.
    user_term : str
        Original user term.
    normalized : str
        Normalized user term.

    Returns
    -------
    list of dict
        Rows with flattened properties, flags, stage, score, and rank.
    """
    best_rows: List[Dict[str, Any]] = []
    for rank_idx, row in enumerate(chosen, start=1):
        props = row.get("properties", {}) if row else {}
        flags = row.get("_flags", {})
        best_rows.append(
            {
                "search_term": user_term,
                "normalized_query": normalized,
                "rank": float(rank_idx),
                "loinc_code": row.get("code"),
                "display": row.get("display"),
                "definition": row.get("definition"),
                "status": row.get("status"),
                "class": props.get("CLASS"),
                "system": props.get("SYSTEM"),
                "property": props.get("PROPERTY"),
                "time": props.get("TIME_ASPCT"),
                "method": props.get("METHOD_TYP"),
                "scale": props.get("SCALE_TYP"),
                "example_units": props.get("EXAMPLE_UCUM_UNITS"),
                "system_core": props.get("system-core"),
                "super_system": props.get("super-system"),
                "time_core": props.get("time-core"),
                "time_modifier": props.get("time-modifier"),
                "analyte": props.get("analyte"),
                "analyte_core": props.get("analyte-core"),
                "analyte_suffix": props.get("analyte-suffix"),
                "analyte_numerator": props.get("analyte-numerator"),
                "analyte_divisor": props.get("analyte-divisor"),
                "analyte_divisor_suffix": props.get("analyte-divisor-suffix"),
                "category": props.get("category"),
                "search_terms": props.get("search"),
                "display_name": props.get("DisplayName"),
                "is_part": flags.get("is_part"),
                "is_answer_list": flags.get("is_answer_list"),
                "is_deprecated": flags.get("is_deprecated"),
                "is_derived": flags.get("is_derived"),
                "is_percentile": flags.get("is_percentile"),
                "has_laterality": flags.get("has_laterality"),
                "property_match": flags.get("property_match"),
                "scale_match": flags.get("scale_match"),
                "stage": row.get("_stage"),
                "score": row.get("_score"),
                "error": None,
            }
        )
    return best_rows


def _all_candidates_rows(
    enriched: List[Dict[str, Any]],
    user_term: str,
    normalized: str,
    numeric_expected: bool,
    expected_props: Optional[Set[str]],
    allow_derived: bool,
    allow_percentile: bool,
) -> List[Dict[str, Any]]:
    """Emit verbose rows for *all* enriched candidates (audit trail).

    Parameters
    ----------
    enriched : list of dict
        Enriched candidates.
    user_term : str
        Original term.
    normalized : str
        Normalized term.
    numeric_expected : bool
        Whether numeric is expected.
    expected_props : set of str or None
        Expected PROPERTY values.
    allow_derived : bool
        Whether derived entries are allowed in strict scoring.
    allow_percentile : bool
        Whether percentile entries are allowed in strict scoring.

    Returns
    -------
    list of dict
        Verbose rows including strict and relaxed scores.
    """
    all_candidates: List[Dict[str, Any]] = []
    for d in enriched:
        flags = compute_feature_flags(d, numeric_expected, expected_props)
        score_strict = final_score(
            d,
            user_term,
            normalized,
            flags,
            allow_derived,
            allow_percentile,
            numeric_expected,
            expected_props,
        )
        score_relaxed = final_score(
            d, user_term, normalized, flags, True, True, False, expected_props
        )
        props = d.get("properties", {}) if d else {}
        all_candidates.append(
            {
                "search_term": user_term,
                "normalized_query": normalized,
                "loinc_code": d.get("code"),
                "display": d.get("display"),
                "status": d.get("status"),
                "definition": d.get("definition"),
                "class": props.get("CLASS"),
                "system": props.get("SYSTEM"),
                "property": props.get("PROPERTY"),
                "time": props.get("TIME_ASPCT"),
                "method": props.get("METHOD_TYP"),
                "scale": props.get("SCALE_TYP"),
                "example_units": props.get("EXAMPLE_UCUM_UNITS"),
                "system_core": props.get("system-core"),
                "super_system": props.get("super-system"),
                "time_core": props.get("time-core"),
                "time_modifier": props.get("time-modifier"),
                "analyte": props.get("analyte"),
                "analyte_core": props.get("analyte-core"),
                "analyte_suffix": props.get("analyte-suffix"),
                "analyte_numerator": props.get("analyte-numerator"),
                "analyte_divisor": props.get("analyte-divisor"),
                "analyte_divisor_suffix": props.get("analyte-divisor-suffix"),
                "category": props.get("category"),
                "search_terms": props.get("search"),
                "display_name": props.get("DisplayName"),
                "is_part": flags["is_part"],
                "is_answer_list": flags["is_answer_list"],
                "is_deprecated": flags["is_deprecated"],
                "is_derived": flags["is_derived"],
                "is_percentile": flags["is_percentile"],
                "has_laterality": flags["has_laterality"],
                "property_match": flags["property_match"],
                "scale_match": flags["scale_match"],
                "score_strict": score_strict,
                "score_relaxed": score_relaxed,
            }
        )
    return all_candidates


# ----------------------------
# Strict/relaxed strategy helpers
# ----------------------------


def _allowances_from_term(user_term: str) -> Tuple[bool, bool]:
    """Derive allowances (derived, percentile) directly from the text term.

    Rules
    -----
    - If term contains ``estimated`` or ``efw`` → allow derived.
    - If term contains ``percentile`` → allow percentile.

    Parameters
    ----------
    user_term : str
        Original term.

    Returns
    -------
    (bool, bool)
        ``(allow_derived, allow_percentile)``.
    """
    lu = user_term.lower()
    allow_derived = (
        ("estimated" in lu) or ("efw" in lu) or ("estimated fetal weight" in lu)
    )
    allow_percentile = "percentile" in lu
    return allow_derived, allow_percentile


def _strict_then_relaxed(
    enriched: List[Dict[str, Any]],
    user_term: str,
    normalized: str,
    numeric_expected: bool,
    expected_props: Optional[Set[str]],
    top_k: int,
    allow_derived: bool,
    allow_percentile: bool,
) -> List[Dict[str, Any]]:
    """Pick top-``k`` with a strict pass, then fill remaining slots from a relaxed pass.

    Parameters
    ----------
    enriched : list of dict
        Enriched candidates.
    user_term : str
        Original term.
    normalized : str
        Normalized term.
    numeric_expected : bool
        Whether numeric is expected.
    expected_props : set of str or None
        Expected PROPERTY values.
    top_k : int
        Number of rows to return.
    allow_derived : bool
        Allow derived in strict mode.
    allow_percentile : bool
        Allow percentile in strict mode.

    Returns
    -------
    list of dict
        Up to ``top_k`` ranked candidates.
    """
    chosen = filter_and_rank(
        enriched,
        user_term,
        normalized,
        numeric_expected=numeric_expected,
        expected_props=expected_props,
        strict=True,
        allow_derived=allow_derived,
        allow_percentile=allow_percentile,
    )[:top_k]

    if len(chosen) < top_k:
        relaxed_ranked = filter_and_rank(
            enriched,
            user_term,
            normalized,
            numeric_expected=numeric_expected,
            expected_props=expected_props,
            strict=False,
            allow_derived=True,
            allow_percentile=True,
        )
        have = {c["code"] for c in chosen}
        for row in relaxed_ranked:
            if len(chosen) >= top_k:
                break
            if row["code"] not in have:
                chosen.append(row)
    return chosen


# ----------------------------
# Term processing
# ----------------------------


def process_one_term(
    session: requests.Session,
    user_term: str,
    count_per_variant: int,
    top_k: int,
    sleep_between_lookups: float,
) -> Dict[str, Any]:
    """Full pipeline for a single term: normalize → variants → expand → lookup → rank.

    Parameters
    ----------
    session : requests.Session
        Authenticated LOINC FHIR session.
    user_term : str
        Free-text term to map to LOINC.
    count_per_variant : int
        Server-side ``$expand`` count hint per variant.
    top_k : int
        Number of results to keep for the 'best' table.
    sleep_between_lookups : float
        Pacing delay between ``$lookup`` requests (seconds).

    Returns
    -------
    dict
        ``{"best_rows": [...], "all_candidates": [...], "normalized": "..."}``
    """
    normalized = normalize_user_term(user_term)
    variants = build_text_variants(user_term, normalized)
    numeric_expected = user_term_implies_numeric(user_term, normalized)
    expected_props = expected_properties_for_intent(user_term, normalized)
    allow_derived, allow_percentile = _allowances_from_term(user_term)

    per_variant_cap = min(count_per_variant, max(200, top_k * 8))

    raw_candidates = gather_raw_candidates(
        session=session,
        variants=variants,
        count_per_variant=count_per_variant,
        per_variant_cap=per_variant_cap,
    )

    if not raw_candidates:
        return {
            "best_rows": [
                _empty_best_row(user_term, normalized, "No candidates returned")
            ],
            "all_candidates": [],
            "normalized": normalized,
        }

    enriched = enrich_candidates(session, raw_candidates, sleep_between_lookups)

    chosen = _strict_then_relaxed(
        enriched=enriched,
        user_term=user_term,
        normalized=normalized,
        numeric_expected=numeric_expected,
        expected_props=expected_props,
        top_k=top_k,
        allow_derived=allow_derived,
        allow_percentile=allow_percentile,
    )

    best_rows = _best_rows_from_ranked(chosen, user_term, normalized)
    if not best_rows:
        best_rows = [
            _empty_best_row(
                user_term, normalized, "No acceptable matches after filtering"
            )
        ]

    all_candidates = _all_candidates_rows(
        enriched=enriched,
        user_term=user_term,
        normalized=normalized,
        numeric_expected=numeric_expected,
        expected_props=expected_props,
        allow_derived=allow_derived,
        allow_percentile=allow_percentile,
    )

    return {
        "best_rows": best_rows,
        "all_candidates": all_candidates,
        "normalized": normalized,
    }


# ----------------------------
# Console preview & warnings
# ----------------------------


def _pick_best_per_term(
    all_best_rows: List[Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """From possibly multiple rows per term, pick the top-ranked one.

    Parameters
    ----------
    all_best_rows : list of dict
        Rows from all terms (may contain multiple ranks per term).

    Returns
    -------
    dict
        Mapping of term → best row.
    """
    best_by_term: Dict[str, Dict[str, Any]] = {}
    for r in all_best_rows:
        key = r["search_term"]
        current_rank = best_by_term.get(key, {}).get("rank") or math.inf
        if r["rank"] is not None and (r["rank"] < current_rank):
            best_by_term[key] = r
        elif key not in best_by_term:
            best_by_term[key] = r
    return best_by_term


def _warning_bits_for_row(term: str, r: Dict[str, Any]) -> List[str]:
    """Collect warning badges for a chosen row (e.g., derived, percentile, scale mismatch)."""
    warn_bits: List[str] = []
    if r.get("is_part"):
        warn_bits.append("LOINC Part (not reportable)")
    if r.get("is_answer_list"):
        warn_bits.append("Answer list (not a measurement)")
    if r.get("is_deprecated"):
        warn_bits.append("DEPRECATED")
    if r.get("is_derived"):
        warn_bits.append("derived/methodized")
    if r.get("is_percentile"):
        warn_bits.append("percentile")

    num_hint = any(u in term.lower() for u in ["(cm)", "(mm)", "(g)", "(bpm)"]) or any(
        k in term.lower()
        for k in ["length", "diameter", "circumference", "weight", "mass", "rate"]
    )
    if num_hint and (not r.get("property_match") or not r.get("scale_match")):
        warn_bits.append("NOT Qn / PROPERTY mismatch for numeric intent")
    return warn_bits


def print_console_preview(all_best_rows: List[Dict[str, Any]]) -> None:
    """Log a human-readable preview (top pick per term) using ``logging``.

    Parameters
    ----------
    all_best_rows : list of dict
        All top-k rows across terms.
    """
    logging.info("=== LOINC lookup preview (best per term) ===")
    best_by_term = _pick_best_per_term(all_best_rows)

    printed_warnings = []
    for term, r in best_by_term.items():
        if r.get("error"):
            logging.info("- %s: ERROR — %s", term, r["error"])
            continue
        code = r.get("loinc_code")
        disp = r.get("display")
        logging.info("- %s: %s — %s", term, code, disp)
        warn_bits = _warning_bits_for_row(term, r)
        if warn_bits:
            printed_warnings.append((term, warn_bits))

    if printed_warnings:
        logging.warning("WARNING: Some 'best' picks have caveats:")
        for term, bits in printed_warnings:
            logging.warning("  • %s: %s", term, "; ".join(bits))
        logging.warning(
            "  (See CSV columns is_derived/is_percentile/scale_match/property_match/stage/score.)"
        )
    logging.info("(See CSV for all matches and properties.)")


# ----------------------------
# CSV writing
# ----------------------------


def write_csv_rows(rows: List[Dict[str, Any]], out_path: str) -> None:
    """Write a list of dict rows to CSV.

    Parameters
    ----------
    rows : list of dict
        Records.
    out_path : str
        Destination file path.

    Notes
    -----
    Uses pandas if available; falls back to the standard library CSV writer.
    """
    if not rows:
        return
    if pd is not None:
        pd.DataFrame(rows).to_csv(out_path, index=False)
    else:
        import csv

        fieldnames = list(rows[0].keys())
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)


# ----------------------------
# Runner
# ----------------------------


def run(
    terms: List[str],
    out_csv: str,
    count_per_variant: int,
    top_k: int,
    sleep_sec: float,
    creds_path: Optional[str],
    save_all_candidates: bool,
    all_candidates_out: str,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Run the full pipeline over one or more terms and write CSV outputs.

    Parameters
    ----------
    terms : list of str
        Search terms; if empty, uses :data:`DEFAULT_TERMS`.
    out_csv : str
        Path for top-k results CSV.
    count_per_variant : int
        ``$expand`` server-side count per variant.
    top_k : int
        Number of top matches to keep per term.
    sleep_sec : float
        Delay between ``$lookup`` calls in seconds.
    creds_path : str or None
        Optional credentials file (two lines).
    save_all_candidates : bool
        If True, write a CSV of all enriched candidates for auditing.
    all_candidates_out : str
        Output path for the all-candidates CSV.

    Returns
    -------
    (list of dict, list of dict)
        Tuple of ``(all_best_rows, all_all_candidates)``.
    """
    user, pw, source = resolve_credentials(creds_path)
    logging.info("[auth] Using credentials from %s", source)
    session = make_requests_session(user, pw)
    try:
        assert_auth_ok(session)
    except (requests.HTTPError, RuntimeError) as exc:
        logging.error("Authentication check failed: %s", exc)
        raise

    all_best_rows: List[Dict[str, Any]] = []
    all_all_candidates: List[Dict[str, Any]] = []

    for term in terms or DEFAULT_TERMS:
        result = process_one_term(session, term, count_per_variant, top_k, sleep_sec)
        all_best_rows.extend(result["best_rows"])
        if save_all_candidates:
            all_all_candidates.extend(result["all_candidates"])

    print_console_preview(all_best_rows)

    if all_best_rows:
        write_csv_rows(all_best_rows, out_csv)
        logging.info("Saved top-%d results to: %s", top_k, out_csv)
    else:
        logging.info("No results to save.")

    if save_all_candidates:
        if all_all_candidates:
            write_csv_rows(all_all_candidates, all_candidates_out)
            logging.info("Saved ALL enriched candidates to: %s", all_candidates_out)
        else:
            logging.info("No 'all candidates' to save (none enriched).")

    return all_best_rows, all_all_candidates

