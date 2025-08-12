#!/usr/bin/env python3
"""
loinc_lookup.py
================

Purpose
-------
Map user-facing ultrasound/OB/GYN terms (e.g., "BPD (cm)", "FHR (bpm)") to
**LOINC observation codes** using the official LOINC FHIR Terminology Service.

This version focuses on:
  • Getting *plain, reportable measurements* (not Parts or AnswerLists)
  • Avoiding deprecated / methodized / derived entries unless you *explicitly* ask
  • Enforcing numeric expectations (PROPERTY + SCALE) when the term implies units
  • Returning *broader* results and *clear flags* so reviewers don't miss issues
  • **Improved recall** for ratio/qualitative terms via curated variants & fallbacks
  • **Deterministic ranking** with secondary ties broken by numeric/Qn preference

Authentication (priority order)
-------------------------------
1) `--creds /path/to/file` (two lines: username, password)   [recommended]
2) Environment variables: `LOINC_USER` and `LOINC_PASS`
3) Local file `./loinc_creds.txt` (two lines: username, password)

Outputs
-------
1) Console preview: best pick per term (+ a warning banner if it’s derived/percentile/etc.)
2) CSV (default `loinc_lookup_results.csv`): the **top-k** per term with rich flags
3) OPTIONAL CSV of **all enriched candidates** per term: `--save-all-candidates`
   (default path `loinc_lookup_all_candidates.csv` unless overridden)

Usage Examples
--------------
# Use defaults and save top-5 per term
python loinc_lookup.py --creds loinc_creds.txt --top-k 5

# Keep top-10, pull up to 100 candidates per text-variant; also save ALL candidates
python loinc_lookup.py --creds loinc_creds.txt --top-k 10 --count 100 --save-all-candidates

# Custom terms; still leverage the gating/scoring
python loinc_lookup.py --creds loinc_creds.txt "BPD (cm)" "HC (cm)" "FHR (bpm)" "EFW (g)"

Notes
-----
- We *do not* hard-limit to prenatal/fetal, but we **softly prefer** fetal/US context.
- If the “best” match is methodized/derived/percentile or not Qn for numeric terms,
  the console prints a **WARNING** and the CSV flags it. This way nobody misses it.
"""

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


# ----------------------------
# FHIR server configuration
# ----------------------------

FHIR_BASE = "https://fhir.loinc.org"
# LOINC's implicit ValueSet for all codes (NOT the older /vs/loinc)
VALUESET_ALL_LOINC = "http://loinc.org/vs"
DEFAULT_SLEEP_SEC = 0.2            # politeness throttle between $lookup calls

# >>> Recall & breadth controls (raised caps; still bounded for runtime sanity)
DEFAULT_MAX_VARIANT_RESULTS = 200  # per-variant trim cap after lightweight scoring (was 20)
GLOBAL_CANDIDATE_CAP = 1200        # guardrail across variants (was 400)


# ----------------------------
# Synonyms and expectations
# ----------------------------

# Basic clinical abbreviations → plain words for search expansion
ABBREV_TO_WORDS = {"BPD": "biparietal diameter", "HC": "head circumference", "AC": "abdominal circumference", "FL": "femur length", "EFW": "estimated fetal weight", "FHR": "fetal heart rate"}

# Soft-preference tokens (we boost but never hard-restrict)
SOFT_CONTEXT_TOKENS = ["fetal", "fetus", "ultrasound", "us", "pregnancy", "obstetric"]

# Hand-curated hints for tricky displays that the filter text alone often misses.
# These give $expand exactly the strings that appear in LOINC displays.
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
DEFAULT_TERMS = ["Gestational Age", "BPD (cm)", "HC (cm)", "AC (cm)", "Femur (cm)", "Cerebellum (cm)", "Cisterna Magna (mm)", "Humerus (cm)", "Radius (cm)", "Ulna (cm)", "Tibia (cm)", "Fibula (cm)", "HC/AC", "FL/AC", "FL/BPD", "EFW (g)", "EFW Percentile", "FHR (bpm)", "Presentation", "Placenta Appearance", "Heart Abnormal", "Head Abnormal", "Face/Neck Abnormal", "Spine Abnormal", "Genitalia Normal"]

# When a term implies a numeric observation, we expect PROPERTY and SCALE_TYP.
# PROPERTY expectations by intent (a *set* so we can accept multiple)
def expected_properties_for_intent(user_term: str, normalized: str) -> Optional[set]:
    """
    Decide expected PROPERTY set based on term intent.

    Parameters
    ----------
    user_term : str
        Original input term (may contain units or abbreviations).
    normalized : str
        Normalized term (abbreviations expanded, ratios verbalized).

    Returns
    -------
    set or None
        A set of acceptable PROPERTY values (e.g., {"Circ"}), or None if not gated.
    """
    text = f"{user_term} {normalized}".lower()
    if "circumference" in text:
        return {"Circ"}
    if "diameter" in text or "bpd" in text:
        # BPD sometimes shows as Length (Len) instead of Diameter (Diam)
        return {"Diam", "Len"}
    if any(k in text for k in ["length", "femur", "humerus", "radius", "ulna", "tibia", "fibula", "long bone"]):
        return {"Len"}
    if "heart rate" in text or "fhr" in text or "rate" in text:
        return {"Rate"}
    if "estimated fetal weight" in text or "efw" in text or "weight" in text or "mass" in text:
        return {"Mass"}
    if "cisterna magna" in text:
        return {"Diam", "Len"}
    if "cerebellum" in text:
        return {"Diam", "Len", "Circ"}
    if " ratio " in text or "/" in strip_units(user_term):
        return {"Rto"}
    # Gestational age: not gating on PROPERTY
    return None


# ----------------------------
# Credentials & HTTP
# ----------------------------

def _tidy(s: str) -> str:
    return (s or "").strip().strip("\r").strip("\n")


def read_creds_from_file(path: str) -> Optional[Tuple[str, str]]:
    """Two-line text file: line1=username (email), line2=password."""
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
    """
    Find LOINC credentials using CLI, env, or local file.

    Returns
    -------
    (user, password, source_desc)
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
    s = requests.Session()
    s.auth = HTTPBasicAuth(user, pw)
    s.headers.update({
        "Accept": "application/fhir+json",
        "User-Agent": "loinc-lookup/2.1 (+local usage)"
    })
    return s


def http_get_json(session: requests.Session, url: str, params: Optional[Dict[str, Any]] = None,
                  max_retries: int = 3, sleep_base: float = 0.5) -> Dict[str, Any]:
    """
    Robust GET wrapper with controlled retries.

    Notes
    -----
    - 401 → clear error with curl hint
    - 429/5xx → exponential backoff
    - else → raise for status, return JSON
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
    _ = http_get_json(session, f"{FHIR_BASE}/CodeSystem",
                      params={"url": "http://loinc.org", "_format": "json"})


# ----------------------------
# Query normalization
# ----------------------------

def strip_units(user_term: str) -> str:
    """Remove trailing '(unit)' to get a clean base expression."""
    t = user_term.strip()
    if "(" in t and ")" in t and t.index("(") < t.index(")"):
        return t[:t.index("(")].strip()
    return t


def normalize_user_term(user_term: str) -> str:
    """
    Expand abbreviations; rewrite ratios into words so free text matching works.
      e.g., "HC/AC" → "head circumference ratio abdominal circumference"
    """
    base = strip_units(user_term)
    upper_no_space = base.upper().replace(" ", "")
    for abbr, words in ABBREV_TO_WORDS.items():
        if upper_no_space == abbr or base.upper() == abbr:
            return words
    if "/" in base:  # ratio shorthand
        a, b = [p.strip() for p in base.split("/", 1)]
        a = ABBREV_TO_WORDS.get(a.upper(), a)
        b = ABBREV_TO_WORDS.get(b.upper(), b)
        return f"{a} ratio {b}"
    return base


def user_term_implies_numeric(user_term: str, normalized: str) -> bool:
    """
    If the term implies a numeric value, we require Qn and PROPERTY in strict mode.

    Signals
    -------
    explicit units; words like length/diameter/circumference/weight/rate.
    """
    t = f"{user_term} {normalized}".lower()
    if any(u in user_term.lower() for u in ["(cm)", "(mm)", "(g)", "(bpm)"]):
        return True
    for k in ["length", "diameter", "circumference", "weight", "mass", "rate"]:
        if k in t:
            return True
    return False


def is_ratio_intent(user_term: str, normalized: str) -> bool:
    return (" ratio " in f"{user_term} {normalized}".lower()) or ("/" in strip_units(user_term))


# ----------------------------
# Phrase variants (critical for bringing in the "right" LOINC candidates)
# ----------------------------

def build_text_variants(user_term: str, normalized: str) -> List[str]:
    """
    Construct a diverse set of LOINC-style phrases to pass to ValueSet/$expand.

    Strategy
    --------
    - Start with normalized text + ultrasound variants
    - Add bracketed unit/quantity hints ([Length], [Diameter], [Circumference], [Mass], [Rate])
    - Inject common LOINC wording patterns (e.g., "Head Diameter.biparietal")
    - Add qualitative presence/morphology phrasing for non-numeric intents
    - **Add literal ratio forms (A/B) and "derived by US" phrasings**
    - **Inject hand-curated HARD_HINTS for tricky, high-value targets**
    """
    terms: List[str] = []

    # Base + ultrasound context (soft)
    terms.append(normalized)
    terms.extend([f"{normalized} ultrasound", f"{normalized} US"])

    # Bracketed style hints seen in LOINC displays
    terms.extend([
        f"{normalized} [Length]",
        f"{normalized} [Diameter]",
        f"{normalized} [Circumference]",
        f"{normalized} [Mass]",
        f"{normalized} [Rate]",
    ])

    ln = normalized.lower()
    lu = user_term.lower()

    # --- Specific measurement families (LOINC-ish phrasing) ---

    # BPD
    if "biparietal diameter" in ln or "bpd" in lu:
        terms.extend([
            "Head Diameter.biparietal [Length] fetus US",
            "Head biparietal diameter [Length] fetus US",
            "Fetal head biparietal diameter [Length] US",
        ])

    # HC
    if "head circumference" in ln or "hc" in lu:
        terms.extend([
            "Head [Circumference] fetus US",
        ])

    # AC
    if "abdominal circumference" in ln or "ac" in lu:
        terms.extend([
            "Abdomen [Circumference] fetus US",
            "Abdominal circumference fetus US",
        ])

    # Long bones
    if any(b in ln or b in lu for b in ["femur", "humerus", "radius", "ulna", "tibia", "fibula"]):
        bone = None
        for b in ["femur", "humerus", "radius", "ulna", "tibia", "fibula"]:
            if b in ln or b in lu:
                bone = b
                break
        if bone:
            terms.extend([
                f"{bone.title()} diaphysis fetus [Length] US",
                f"Fetal {bone} length [Length] US",
            ])

    # Cerebellum / Cisterna magna
    if "cerebellum" in ln:
        terms.extend([
            "Cerebellum fetus [Diameter] US",
            "Cerebellum Diameter transverse fetus US",
        ])
    if "cisterna magna" in ln:
        terms.extend([
            "Cisterna magna fetus [Diameter] US",
            "Fetal cisterna magna sagittal diameter US",
        ])

    # --- Ratios (spell out both ways + literal) ---
    if is_ratio_intent(user_term, normalized):
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
        if a and b:
            # spelled-out variations
            terms.extend([
                f"{a}/{b}",
                f"{a}/{b} ratio",
                f"{a} to {b} ratio",
                f"ratio of {a} to {b}",
                f"{a} divided by {b}",
                f"{a} over {b}",
                f"{a.title()} / {b.title()} derived by US",
                f"{a.title()} / {b.title()} derived by ultrasound",
            ])
            # keep the literal user token too (e.g., "HC/AC")
            raw_lit = strip_units(user_term)
            terms.append(raw_lit)
            terms.append(raw_lit.upper())

    # EFW – focus on mass, not percentile/panels
    if "estimated fetal weight" in ln or "efw" in lu:
        terms.extend([
            "Estimated fetal weight [Mass] US",
            "Fetal body weight [Mass] US",
            "Fetus body weight [Mass] US",
            "EFW [Mass] US",
        ])

    # FHR – include US and auscultation variants
    if "fetal heart rate" in ln or "fhr" in lu:
        terms.extend([
            "Fetal heart rate [Rate] US",
            "Fetal heart rate US",
            "Fetal heart rate by auscultation",
            "Fetal heart rate mean 10 minutes",
            "Fetal heart rate reactivity",
        ])

    # Qualitative presence/morphology variants (beyond HARD_HINTS; keep generic)
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
    for k, variants in qual_map.items():
        if k in lu:
            terms.extend(variants)

    # --- Inject hand-curated HARD_HINTS ---
    key = strip_units(user_term).lower()
    key_compact = key.replace(" ", "")
    terms.extend(HARD_HINTS.get(key, []))
    terms.extend(HARD_HINTS.get(key_compact, []))

    # Deduplicate while preserving order
    seen = set()
    uniq: List[str] = []
    for t in terms:
        if t not in seen:
            uniq.append(t)
            seen.add(t)
    return uniq


# ----------------------------
# FHIR: search and enrichment
# ----------------------------

def _fallback_variant(text_filter: str) -> Optional[str]:
    """
    Loosen a variant string if the first $expand returns nothing.

    The fallback removes bracketed tokens and punctuation (except slashes),
    collapses whitespace, and returns None if no meaningful change occurred.
    """
    v2 = re.sub(r"\[[^\]]+\]", "", text_filter)
    v2 = re.sub(r"[(),.:;]", " ", v2)
    v2 = re.sub(r"\s+", " ", v2).strip()
    if len(v2) >= 3 and v2 != text_filter:
        return v2
    return None


def expand_valueset_candidates(session: requests.Session, text_filter: str, count: int) -> List[Dict[str, Any]]:
    """
    Call ValueSet/$expand and return the raw 'contains' entries.

    Implements a lightweight fallback if the first call returns nothing.
    """
    params = {"url": VALUESET_ALL_LOINC, "filter": text_filter, "count": count, "_format": "json"}
    data = http_get_json(session, f"{FHIR_BASE}/ValueSet/$expand", params)
    contains = (data.get("expansion", {}) or {}).get("contains", []) or []
    if not contains:
        fb = _fallback_variant(text_filter)
        if fb:
            params["filter"] = fb
            data = http_get_json(session, f"{FHIR_BASE}/ValueSet/$expand", params)
            contains = (data.get("expansion", {}) or {}).get("contains", []) or []
    return contains


def lookup_loinc_details(session: requests.Session, code: str) -> Dict[str, Any]:
    """
    Flatten CodeSystem/$lookup (Parameters) to a dict with properties.

    Returns
    -------
    dict
        Keys include: "code", "display", "definition", "status", "properties".
    """
    params = {"system": "http://loinc.org", "code": code, "_format": "json"}
    data = http_get_json(session, f"{FHIR_BASE}/CodeSystem/$lookup", params)
    out: Dict[str, Any] = {"code": code}
    for p in data.get("parameter", []):
        name = p.get("name")
        val = p.get("valueString") or p.get("valueCode") or p.get("valueUri")
        if name in ("name", "display"):
            out["display"] = val
        elif name == "version":
            out["version"] = val
        elif name == "status":
            out["status"] = val
        elif name == "abstract":
            out["abstract"] = val
        elif name == "definition":
            out["definition"] = val
        elif name == "property":
            prop_code = None
            prop_value = None
            for sub in p.get("part", []):
                if sub.get("name") == "code":
                    prop_code = sub.get("valueCode")
                elif sub.get("name") in ("valueString", "valueCode", "valueUri"):
                    prop_value = sub.get("valueString") or sub.get("valueCode") or sub.get("valueUri")
            if prop_code and prop_value:
                out.setdefault("properties", {})[prop_code] = prop_value
    return out


# ----------------------------
# Classification flags & acceptance checks
# ----------------------------

def is_loinc_part(code: Optional[str]) -> bool:
    """LOINC Part codes start with 'LP' and are *not* reportable observations."""
    return bool(code) and code.startswith("LP")


def is_loinc_answer_list(code: Optional[str]) -> bool:
    """LOINC Answer List/Answer codes start with 'LA' and are *not* observations."""
    return bool(code) and code.startswith("LA")


def is_deprecated(status: Optional[str], display: Optional[str]) -> bool:
    if status and status.upper() == "DEPRECATED":
        return True
    return bool(display) and "deprecated" in display.lower()


def looks_derived_or_methodized(display: Optional[str]) -> bool:
    """
    Heuristics for "derived/method" wording we want to avoid unless requested.
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
    return bool(display) and "percentile" in display.lower()


def mentions_laterality(display: Optional[str]) -> bool:
    if not display:
        return False
    d = display.lower()
    return any(w in d for w in [" left ", " right ", " left]", " right]", "left ", "right "])


def _norm_str(x: Optional[str]) -> str:
    """
    Normalize simple code/property strings for robust comparisons.
    """
    return (x or "").strip()


def compute_feature_flags(
    details: Dict[str, Any],
    numeric_expected: bool,
    expected_props: Optional[set],
) -> Dict[str, Any]:
    """
    Derive useful booleans/fields from $lookup 'details' so reviewers can see
    exactly *why* a candidate was (not) accepted or penalized.

    Returns
    -------
    dict
        booleans for is_part/is_answer_list/is_deprecated/is_derived/is_percentile,
        laterality, and property/scale matches.
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
    expected_norm = {(_norm_str(p)) for p in (expected_props or set())}

    # Acceptance checks for strict stage when numeric is expected
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


# ----------------------------
# Scoring and filtering
# ----------------------------

def prefilter_score_for_variant(display: str, variant_text: str) -> int:
    """Lightweight score to trim each variant's raw results before $lookup."""
    disp = (display or "").lower()
    s = 0
    for w in SOFT_CONTEXT_TOKENS:
        if w in disp:
            s += 2
    for tok in variant_text.lower().split():
        if tok and tok in disp:
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
    expected_props: Optional[set],
) -> int:
    """
    Deterministic feature-based scoring:

    Positive
    --------
    + expected PROPERTY (big)
    + SCALE_TYP=Qn for numeric
    + ultrasound/fetal context
    + minor token overlaps
    + small universal nudge for Qn

    Negative
    --------
    - derived/method unless allowed
    - percentile unless allowed
    - laterality unless requested
    - narrative/ordinal signal when numeric expected
    """
    disp = (details.get("display") or "").lower()
    props = details.get("properties", {}) if "properties" in details else {}
    cls = (props.get("CLASS") or "").lower()
    system = (props.get("SYSTEM") or "").lower()
    sys_core = (props.get("system-core") or "").lower()
    super_sys = (props.get("super-system") or "").lower()
    status = details.get("status") or ""

    s = 0

    # PROPERTY/SCALE
    if expected_props and props.get("PROPERTY") in expected_props:
        s += 12  # dominant signal for intent
    if numeric_expected and props.get("SCALE_TYP") == "Qn":
        s += 6

    # Universal nudge toward Qn (helps FHR numeric forms beat reactivity)
    if props.get("SCALE_TYP") == "Qn":
        s += 2

    # Context (soft)
    if "us" in cls or "ob" in cls or "ultrasound" in disp or " us" in disp:
        s += 3
    if "fetal" in disp or "fetus" in disp or super_sys == "fetus" or "fetus" in sys_core or "fetus" in system:
        s += 2

    # Penalties
    if not allow_derived and flags["is_derived"]:
        s -= 10
    if allow_derived and flags["is_derived"]:
        s -= 4

    if not allow_percentile and flags["is_percentile"]:
        s -= 8
    if allow_percentile and flags["is_percentile"]:
        s -= 2

    want_laterality = any(w in user_term.lower() for w in ["left", "right"])
    if flags["has_laterality"] and not want_laterality:
        s -= 3

    # Narrative/ordinal penalty if numeric intent
    if numeric_expected and ((props.get("SCALE_TYP") in {"Ord", "Nar"}) or "narrative" in disp):
        s -= 4

    # Token overlap
    for tok in normalized.lower().split():
        if tok and tok in disp:
            s += 1

    # Safety nudge
    code = details.get("code") or ""
    if not code.startswith(("LP", "LA")):
        s += 1
    if not is_deprecated(status, details.get("display")):
        s += 1

    return s


def filter_and_rank(
    enriched: List[Dict[str, Any]],
    user_term: str,
    normalized: str,
    numeric_expected: bool,
    expected_props: Optional[set],
    strict: bool,
    allow_derived: bool,
    allow_percentile: bool,
) -> List[Dict[str, Any]]:
    """
    STRICT:
      - Exclude Parts/AnswerLists/Deprecated
      - Exclude derived/method unless allow_derived
      - Exclude percentile unless allow_percentile
      - If numeric_expected: require SCALE_TYP=Qn AND PROPERTY in expected set (if provided)
    RELAXED:
      - Still exclude Parts/AnswerLists/Deprecated
      - Allow derived/percentile (but score penalizes them)
      - Drop the hard requirement on PROPERTY/SCALE (but scoring still prefers them)
    """
    accepted: List[Dict[str, Any]] = []

    for d in enriched:
        flags = compute_feature_flags(d, numeric_expected, expected_props)

        # Hard drops common to both modes
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
                # If we know what PROPERTY we want, require it strictly
                if expected_props and not flags["property_match"]:
                    continue

        # Score and annotate
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
        row = dict(d)  # shallow copy
        row["_flags"] = flags
        row["_score"] = s
        row["_stage"] = "strict" if strict else "relaxed"
        accepted.append(row)

    # Sort by main score with deterministic tie-breakers favoring numeric/Qn
    def _key(x: Dict[str, Any]):
        f = x["_flags"]
        # prefer property/scale matches when numeric; prefer not-derived, not-percentile
        sec = (
            1 if f.get("property_match") else 0,
            1 if f.get("scale_match") else 0,
            0 if f.get("is_derived") else 1,
            0 if f.get("is_percentile") else 1,
        )
        return (x["_score"],) + sec

    accepted.sort(key=_key, reverse=True)
    return accepted


# ----------------------------
# Term processing (search → enrich → strict → relaxed)
# ----------------------------

def process_one_term(
    session: requests.Session,
    user_term: str,
    count_per_variant: int,
    top_k: int,
    sleep_between_lookups: float,
) -> Dict[str, Any]:
    """
    Run the full pipeline for a single term.

    Returns
    -------
    dict
        {
          "best_rows": [top-k dicts ready for CSV],
          "all_candidates": [all enriched dicts + flags + scores],
          "normalized": "<normalized text>"
        }
    """
    normalized = normalize_user_term(user_term)
    variants = build_text_variants(user_term, normalized)
    numeric_expected = user_term_implies_numeric(user_term, normalized)
    expected_props = expected_properties_for_intent(user_term, normalized)

    # Heuristic: allowing derived/percentile only if the user plainly asks for it
    lu = user_term.lower()
    allow_derived = ("estimated" in lu) or ("efw" in lu) or ("estimated fetal weight" in lu)
    allow_percentile = ("percentile" in lu)

    # Gather raw candidates across variants (dedup codes)
    raw_candidates: Dict[str, Dict[str, Any]] = {}

    # Per-variant trim scales with caller count and top_k
    per_variant_cap = min(count_per_variant, max(200, top_k * 8))

    for v in variants:
        try:
            contains = expand_valueset_candidates(session, v, count=count_per_variant)
        except Exception:
            continue
        # lightweight trimming per variant
        contains = sorted(contains, key=lambda c: prefilter_score_for_variant(c.get("display", ""), v), reverse=True)
        for c in contains[:max(DEFAULT_MAX_VARIANT_RESULTS, per_variant_cap)]:
            code = c.get("code")
            if not code or code in raw_candidates:
                continue
            raw_candidates[code] = {"code": code, "display": c.get("display")}

        if len(raw_candidates) >= GLOBAL_CANDIDATE_CAP:
            break

    if not raw_candidates:
        return {
            "best_rows": [{
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
                "error": "No candidates returned"
            }],
            "all_candidates": [],
            "normalized": normalized,
        }

    # Enrich with $lookup
    enriched: List[Dict[str, Any]] = []
    for code in raw_candidates:
        time.sleep(sleep_between_lookups)  # politeness
        try:
            d = lookup_loinc_details(session, code)
            enriched.append(d)
        except Exception:
            # skip individual lookup errors; keep going
            continue

    # STRICT filter & rank
    strict_ranked = filter_and_rank(
        enriched, user_term, normalized,
        numeric_expected=numeric_expected, expected_props=expected_props,
        strict=True, allow_derived=allow_derived, allow_percentile=allow_percentile
    )
    chosen = strict_ranked[:top_k]

    # RELAXED fallback (still excludes LP/LA/deprecated; allows derived/percentile)
    if len(chosen) < top_k:
        relaxed_ranked = filter_and_rank(
            enriched, user_term, normalized,
            numeric_expected=numeric_expected, expected_props=expected_props,
            strict=False, allow_derived=True, allow_percentile=True
        )
        # keep adding until top_k
        for row in relaxed_ranked:
            if len(chosen) >= top_k:
                break
            # Avoid duplicates if any code was already chosen
            if row["code"] not in {c["code"] for c in chosen}:
                chosen.append(row)

    # Prepare unified "all candidates" (strict + relaxed scores/flags)
    all_candidates: List[Dict[str, Any]] = []
    for d in enriched:
        flags = compute_feature_flags(d, numeric_expected, expected_props)
        score_strict = final_score(
            details=d,
            user_term=user_term,
            normalized=normalized,
            flags=flags,
            allow_derived=allow_derived,
            allow_percentile=allow_percentile,
            numeric_expected=numeric_expected,
            expected_props=expected_props,
        )
        # Also compute a relaxed score for comparison
        score_relaxed = final_score(
            details=d,
            user_term=user_term,
            normalized=normalized,
            flags=flags,
            allow_derived=True,
            allow_percentile=True,
            numeric_expected=False,        # relaxed ignores numeric gate in score
            expected_props=expected_props,
        )
        props = d.get("properties", {}) if d else {}
        all_candidates.append({
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
        })

    # Convert chosen rows to output schema (top-k only)
    best_rows: List[Dict[str, Any]] = []
    for rank_idx, row in enumerate(chosen, start=1):
        props = row.get("properties", {}) if row else {}
        flags = row.get("_flags", {})
        best_rows.append({
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
            # Flags & scoring
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
            "error": None
        })

    # If *no* acceptable matches after strict+relaxed
    if not best_rows:
        best_rows.append({
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
            "error": "No acceptable matches after filtering"
        })

    return {"best_rows": best_rows, "all_candidates": all_candidates, "normalized": normalized}


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
    """
    Execute the lookup for all terms and write CSVs.

    Parameters
    ----------
    terms : list of str
        Input terms. If empty, DEFAULT_TERMS are used.
    out_csv : str
        Output CSV path for top-k results.
    count_per_variant : int
        Server-side candidate count per text variant for $expand.
    top_k : int
        How many top matches to keep per term.
    sleep_sec : float
        Seconds to sleep between $lookup calls.
    creds_path : str or None
        Path to two-line creds file.
    save_all_candidates : bool
        Whether to also save a CSV of all enriched candidates.
    all_candidates_out : str
        Output CSV path for the "all candidates" file.

    Returns
    -------
    (best_rows_all_terms, all_candidates_all_terms)
    """
    user, pw, source = resolve_credentials(creds_path)
    print(f"[auth] Using credentials from {source}")
    session = make_requests_session(user, pw)
    try:
        assert_auth_ok(session)
    except Exception:
        print("Authentication check failed.\n")
        raise

    all_best_rows: List[Dict[str, Any]] = []
    all_all_candidates: List[Dict[str, Any]] = []

    for term in (terms or DEFAULT_TERMS):
        result = process_one_term(session, term, count_per_variant, top_k, sleep_sec)
        best_rows = result["best_rows"]
        all_cands = result["all_candidates"]

        all_best_rows.extend(best_rows)
        if save_all_candidates:
            all_all_candidates.extend(all_cands)

    # -------- Console preview (best per term) --------
    print("\n=== LOINC lookup preview (best per term) ===")
    # pick the lowest rank (1 is best) per search_term
    best_by_term: Dict[str, Dict[str, Any]] = {}
    for r in all_best_rows:
        key = r["search_term"]
        current_rank = (best_by_term.get(key, {}).get("rank") or math.inf)
        if r["rank"] is not None and (r["rank"] < current_rank):
            best_by_term[key] = r
        elif key not in best_by_term:
            best_by_term[key] = r

    # Print and flag suspicious bests
    printed_warnings = []
    for term, r in best_by_term.items():
        if r.get("error"):
            print(f"- {term}: ERROR — {r['error']}")
            continue
        code = r.get("loinc_code")
        disp = r.get("display")
        print(f"- {term}: {code} — {disp}")

        # Warning banner when best pick is questionable for clinical reporting
        warn_bits = []
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
        # If numeric intent but property/scale mismatch, flag loudly
        num_hint = any(u in term.lower() for u in ["(cm)", "(mm)", "(g)", "(bpm)"]) or \
                   any(k in term.lower() for k in ["length", "diameter", "circumference", "weight", "mass", "rate"])
        if num_hint and (not r.get("property_match") or not r.get("scale_match")):
            warn_bits.append("NOT Qn / PROPERTY mismatch for numeric intent")

        if warn_bits:
            printed_warnings.append((term, warn_bits))

    if printed_warnings:
        print("\nWARNING: Some 'best' picks have caveats:")
        for term, bits in printed_warnings:
            print(f"  • {term}: " + "; ".join(bits))
        print("  (See CSV columns is_derived/is_percentile/scale_match/property_match/stage/score.)")

    print("\n(See CSV for all matches and properties.)\n")

    # -------- Write CSVs --------
    if not all_best_rows:
        print("No results to save.")
        return all_best_rows, all_all_candidates

    if pd is not None:
        pd.DataFrame(all_best_rows).to_csv(out_csv, index=False)
    else:
        import csv
        fieldnames = list(all_best_rows[0].keys())
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader(); w.writerows(all_best_rows)
    print(f"Saved top-{top_k} results to: {out_csv}")

    if save_all_candidates:
        if all_all_candidates:
            if pd is not None:
                pd.DataFrame(all_all_candidates).to_csv(all_candidates_out, index=False)
            else:
                import csv
                fieldnames = list(all_all_candidates[0].keys())
                with open(all_candidates_out, "w", newline="", encoding="utf-8") as f:
                    w = csv.DictWriter(f, fieldnames=fieldnames)
                    w.writeheader(); w.writerows(all_all_candidates)
            print(f"Saved ALL enriched candidates to: {all_candidates_out}")
        else:
            print("No 'all candidates' to save (none enriched).")

    return all_best_rows, all_all_candidates


# ----------------------------
# CLI
# ----------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Look up LOINC codes/labels/definitions via the official LOINC FHIR Terminology Service."
    )
    p.add_argument("terms", nargs="*", help="Search terms. If omitted, uses an obstetric ultrasound default list.")
    p.add_argument("--out", default="loinc_lookup_results.csv",
                   help="Output CSV path for top-k results (default: loinc_lookup_results.csv)")
    p.add_argument("--count", type=int, default=50,
                   help="Server-side candidate count per text variant for $expand (default: 50)")
    p.add_argument("--top-k", type=int, default=5,
                   help="How many top matches to keep per term (default: 5)")
    p.add_argument("--sleep", type=float, default=DEFAULT_SLEEP_SEC,
                   help=f"Seconds to sleep between $lookup calls (default: {DEFAULT_SLEEP_SEC})")
    p.add_argument("--creds", default=None,
                   help=f"Path to creds file (two lines: username, password). If omitted, uses env or ./loinc_creds.txt")
    p.add_argument("--save-all-candidates", action="store_true",
                   help="Also save a CSV with ALL enriched candidates and their flags/scores (audit trail).")
    p.add_argument("--all-out", default="loinc_lookup_all_candidates.csv",
                   help="Output CSV path for ALL enriched candidates (default: loinc_lookup_all_candidates.csv)")
    return p.parse_args()


def main():
    args = parse_args()
    terms = args.terms if args.terms else DEFAULT_TERMS
    run(
        terms=terms,
        out_csv=args.out,
        count_per_variant=args.count,
        top_k=args.top_k,
        sleep_sec=args.sleep,
        creds_path=args.creds,
        save_all_candidates=args.save_all_candidates,
        all_candidates_out=args.all_out,
    )


if __name__ == "__main__":
    main()
