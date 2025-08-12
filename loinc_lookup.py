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

For Clinicians
--------------
- Results are labeled and flagged clearly:
    is_part, is_answer_list, is_deprecated, is_derived, is_percentile,
    has_laterality, property_match, scale_match, stage (strict|relaxed), score.
- Numeric terms (cm/mm/g/bpm; “length/diameter/circumference/weight/rate”)
  prefer **SCALE_TYP=Qn** and the expected **PROPERTY** (e.g., Circ, Diam, Len, Mass, Rate).
- Fallback (RELAXED) still excludes Parts/AnswerLists/deprecated but allows
  methodized/derived/percentile if nothing else is available—**and flags them.**

For Bioinformaticians
---------------------
- Candidates come from `ValueSet/$expand` on the implicit LOINC ValueSet
  `http://loinc.org/vs` using multiple **LOINC-style phrase variants**.
- Details from `CodeSystem/$lookup` populate properties like PROPERTY, SCALE_TYP,
  SYSTEM, METHOD_TYP, and more. We score + filter using these properties.
- We provide an **optional all-candidates CSV** (`--save-all-candidates`)
  containing every enriched candidate with scores and gating flags for audit.

For Computer Scientists
-----------------------
- Clear separation of concerns (auth, HTTP, candidates, enrichment, gating, ranking).
- Idempotent GET with exponential backoff on 429/5xx; fast-fail on 401 with hints.
- Deterministic ranking with explicit feature-based scoring and **stable tie-breaks**.
- CSV outputs are schema-stable; columns are documented in code.

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
