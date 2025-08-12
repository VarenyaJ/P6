#!/usr/bin/env python3
"""
Filtered LOINC lookup for fetal biometry headers (OB ultrasound–focused).

This script searches a local LOINC 2.80 `LoincTableCore.csv` for codes related to
column headers in an Excel sheet (default: 'fetal_biometry'). It returns the requested
LOINC metadata while filtering to fetal/OB ultrasound context to avoid noisy matches
(e.g., adult EKG heart-rate).

Core features
-------------
1) Header term expansion:
   - Normalizes column headers (lowercases, removes units in parentheses).
   - Splits into tokens and adds helpful aliases (e.g., "bpd" → "biparietal diameter").

2) Targeted search fields:
   - Searches within COMPONENT, LONG_COMMON_NAME, SHORTNAME using a single
     word-bounded, non-capturing regex (no "match group" warnings).

3) Fetal/OB context filtering (default ON):
   - Keeps rows that are strongly OB/fetal (e.g., CLASS in {OB.US, PANEL.OB.US},
     SYSTEM mentions fetal/placental structures, METHOD_TYP includes US, or the text
     contains "fetal"/"fetus").

4) Scoring and trimming:
   - Ranks candidates by fetal context strength, header-term presence, and
     property hints (e.g., Len for BPD/HC/AC; Mass for EFW; NRat for FHR).
   - Keeps top N rows per header (default 50).

Outputs
-------
CSV with columns:
  Header, LOINC_NUM, COMPONENT, PROPERTY, TIME_ASPCT, SYSTEM, SCALE_TYP,
  METHOD_TYP, CLASS, CLASSTYPE, LONG_COMMON_NAME, SHORTNAME
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set

import pandas as pd


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


def normalize_header_to_term(header_text: object) -> str:
    """
    Normalize a column header into a search-friendly base term.

    - Converts to lowercase.
    - Removes units or hints in parentheses (e.g., "(cm)").
    - Collapses multiple spaces.

    Parameters
    ----------
    header_text : object
        The header name (usually a string). Non-strings return "".

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


def collect_search_terms_for_header(header_text: str) -> List[str]:
    """
    Build a list of search terms for a given header using:
      - the normalized base term,
      - tokenized words from the base term,
      - helpful aliases for the base term and its first token.

    Parameters
    ----------
    header_text : str
        The column header from the Excel sheet.

    Returns
    -------
    List[str]
        Unique, non-empty search terms.
    """
    base = normalize_header_to_term(header_text)
    terms: Set[str] = set()
    if base:
        terms.add(base)
        # token-level terms (e.g., "cisterna magna" → {"cisterna", "magna"})
        for token in re.split(r"[\/\-\s]+", base):
            token = token.strip()
            if token:
                terms.add(token)

        # alias expansion on base and on first token
        if base in HEADER_ALIASES:
            terms.update(HEADER_ALIASES[base])
        first_token = base.split()[0] if base else ""
        if first_token in HEADER_ALIASES:
            terms.update(HEADER_ALIASES[first_token])

    return sorted(t for t in terms if t)


def build_search_regex(terms: Iterable[str]) -> str:
    """
    Construct a single, word-bounded, non-capturing regex that matches
    any of the given terms (case-insensitive via pandas).

    Parameters
    ----------
    terms : Iterable[str]
        The search terms.

    Returns
    -------
    str
        A regex string suitable for `Series.str.contains(..., regex=True)`.
    """
    escaped_terms = [re.escape(t) for t in terms if t]
    if not escaped_terms:
        # A regex that never matches (used when no terms were generated)
        return r"(?!x)x"
    return r"\b(?:%s)\b" % "|".join(escaped_terms)


def row_has_ob_fetal_context(row: pd.Series) -> bool:
    """
    Heuristic: does a LOINC row clearly indicate OB/fetal ultrasound context?

    Signals (any is sufficient):
      - CLASS is OB.US or PANEL.OB.US
      - SYSTEM contains fetal/placental structures or related anatomy
      - METHOD_TYP contains 'US'
      - Text contains 'fetal'/'fetus' in COMPONENT/LONG_COMMON_NAME/SHORTNAME

    Parameters
    ----------
    row : pd.Series
        A single row from the LOINC table.

    Returns
    -------
    bool
        True if the row is likely OB/fetal; False otherwise.
    """
    loinc_class = (row.get("CLASS") or "").upper()
    system_text = (row.get("SYSTEM") or "").lower()
    component_text = (row.get("COMPONENT") or "").lower()
    long_name_text = (row.get("LONG_COMMON_NAME") or "").lower()
    short_name_text = (row.get("SHORTNAME") or "").lower()
    method_text = (row.get("METHOD_TYP") or "").upper()

    contains_fetal_word = (
        "fetal" in component_text
        or "fetal" in long_name_text
        or "fetal" in short_name_text
    )
    class_is_ob = loinc_class in OB_ULTRASOUND_CLASSES
    method_is_ultrasound = "US" in method_text
    system_is_ob = any(keyword in system_text for keyword in OB_SYSTEM_KEYWORDS)
    system_mentions_fetus = "fetus" in system_text

    return (
        contains_fetal_word
        or class_is_ob
        or (method_is_ultrasound and system_is_ob)
        or system_mentions_fetus
    )


def score_candidate_row(
    row: pd.Series,
    header_terms_lower: Set[str],
    allowed_properties: Set[str],
) -> int:
    """
    Assign a relevance score to a candidate LOINC row.

    Higher scores reflect stronger fetal/OB context and closer relation to the header.

    Scoring signals
    ---------------
    +3 if SYSTEM contains 'fetus'
    +3 if any of COMPONENT/LONG_COMMON_NAME/SHORTNAME contains 'fetal'
    +2 if CLASS is OB.US or PANEL.OB.US
    +2 if METHOD_TYP contains 'US'
    +2 per appearance of any header term in COMPONENT/LONG_COMMON_NAME/SHORTNAME
    +1 if PROPERTY matches a header-specific hint (e.g., Len/Mass/NRat)

    Parameters
    ----------
    row : pd.Series
        The LOINC row being scored.
    header_terms_lower : Set[str]
        Header-derived terms, all lowercased.
    allowed_properties : Set[str]
        PROPERTY values considered especially relevant for this header.

    Returns
    -------
    int
        The computed score (non-negative).
    """
    score = 0

    component_text = (row.get("COMPONENT") or "").lower()
    long_name_text = (row.get("LONG_COMMON_NAME") or "").lower()
    short_name_text = (row.get("SHORTNAME") or "").lower()
    system_text = (row.get("SYSTEM") or "").lower()
    loinc_class = (row.get("CLASS") or "").upper()
    method_text = (row.get("METHOD_TYP") or "").upper()
    property_value = (row.get("PROPERTY") or "")

    # Context signals
    if "fetus" in system_text:
        score += 3
    if "fetal" in component_text or "fetal" in long_name_text or "fetal" in short_name_text:
        score += 3
    if loinc_class in OB_ULTRASOUND_CLASSES:
        score += 2
    if "US" in method_text:
        score += 2

    # Header-term presence
    for term in header_terms_lower:
        if f" {term} " in f" {component_text} " or f" {term} " in f" {long_name_text} " or f" {term} " in f" {short_name_text} ":
            score += 2

    # Property hint (if applicable)
    if allowed_properties and property_value in allowed_properties:
        score += 1

    return score


def main() -> None:
    """
    Parse CLI arguments, perform the LOINC search, and write the output CSV.

    Steps
    -----
    1) Load LOINC CSV with string dtype.
    2) Read the Excel sheet and extract column headers to search for.
    3) For each header:
       - Generate search terms and regex.
       - Search COMPONENT/LONG_COMMON_NAME/SHORTNAME.
       - Optionally filter to OB/fetal context.
       - Score, sort, and cap results per header.
    4) Concatenate and write the final CSV.
    """
    parser = argparse.ArgumentParser(
        description="Search local LOINC for codes related to fetal biometry headers "
                    "(OB ultrasound–focused)."
    )
    parser.add_argument(
        "--loinc-csv",
        required=True,
        help="Path to LoincTableCore.csv (e.g., /Users/you/Downloads/Loinc_2.80/LoincTableCore.csv)",
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
        help="Print per-header diagnostics (raw/filtered/kept counts).",
    )
    args = parser.parse_args()

    loinc_csv_path = Path(args.loinc_csv)
    excel_path = Path(args.excel)
    if not loinc_csv_path.exists():
        raise FileNotFoundError(f"LOINC CSV not found: {loinc_csv_path}")
    if not excel_path.exists():
        raise FileNotFoundError(f"Excel file not found: {excel_path}")

    # Load LOINC with consistent string dtype
    loinc_dataframe = pd.read_csv(loinc_csv_path, dtype=str, low_memory=False).fillna("")
    missing_columns = [c for c in REQUESTED_OUTPUT_COLUMNS if c not in loinc_dataframe.columns]
    if missing_columns:
        raise SystemExit(f"Missing expected columns in LOINC CSV: {missing_columns}")

    # Load headers from Excel sheet
    excel_file = pd.ExcelFile(excel_path)
    if args.sheet not in excel_file.sheet_names:
        raise SystemExit(f"Sheet '{args.sheet}' not found in {excel_path}. "
                         f"Found: {excel_file.sheet_names}")
    fetal_sheet_dataframe = excel_file.parse(args.sheet)
    excel_headers: List[str] = list(fetal_sheet_dataframe.columns)

    fields_to_search = ["COMPONENT", "LONG_COMMON_NAME", "SHORTNAME"]
    all_header_results: List[pd.DataFrame] = []

    for header in excel_headers:
        search_terms = collect_search_terms_for_header(header)
        if not search_terms:
            if args.debug:
                print(f"[DEBUG] Header '{header}': no usable search terms after normalization.")
            continue

        regex = build_search_regex(search_terms)

        # Boolean mask across the target fields (case-insensitive)
        search_mask = pd.Series(False, index=loinc_dataframe.index)
        for field in fields_to_search:
            search_mask = search_mask | loinc_dataframe[field].str.contains(
                regex, case=False, na=False, regex=True
            )

        candidate_rows = loinc_dataframe.loc[search_mask].copy()
        if candidate_rows.empty:
            if args.debug:
                print(f"[DEBUG] Header '{header}': 0 raw matches.")
            continue

        raw_count = len(candidate_rows)

        # Context filter to OB/fetal ultrasound
        if not args.no_context_filter:
            context_mask = candidate_rows.apply(row_has_ob_fetal_context, axis=1)
            candidate_rows = candidate_rows.loc[context_mask]

        filtered_count = len(candidate_rows)
        if candidate_rows.empty:
            if args.debug:
                print(f"[DEBUG] Header '{header}': 0 matches after OB/fetal context filter.")
            continue

        # PROPERTY hints for this header (use both the full normalized header and first token)
        normalized_header = normalize_header_to_term(header)
        hint_keys: Set[str] = {normalized_header}
        if normalized_header:
            hint_keys.add(normalized_header.split()[0])
        allowed_properties_for_header: Set[str] = set()
        for key in hint_keys:
            allowed_properties_for_header |= PROPERTY_HINTS_BY_HEADER.get(key, set())

        # Score, sort, cap
        header_terms_lower = {t.lower() for t in search_terms}
        candidate_rows["__score"] = candidate_rows.apply(
            lambda r: score_candidate_row(r, header_terms_lower, allowed_properties_for_header),
            axis=1,
        )
        candidate_rows = candidate_rows.sort_values(
            ["__score", "CLASS", "SYSTEM", "COMPONENT"],
            ascending=[False, True, True, True],
        )

        if args.max_per_header and args.max_per_header > 0:
            candidate_rows = candidate_rows.head(args.max_per_header)

        kept_count = len(candidate_rows)

        # Finalize columns and add the originating header
        candidate_rows = candidate_rows.drop(columns=["__score"], errors="ignore")
        candidate_rows.insert(0, "Header", header)
        candidate_rows = candidate_rows[["Header"] + REQUESTED_OUTPUT_COLUMNS].drop_duplicates()

        if args.debug:
            print(
                f"[DEBUG] Header '{header}': raw={raw_count}, "
                f"after_context={filtered_count}, kept={kept_count}, "
                f"allowed_props={sorted(allowed_properties_for_header) if allowed_properties_for_header else 'none'}"
            )

        all_header_results.append(candidate_rows)

    if all_header_results:
        output_dataframe = pd.concat(all_header_results, ignore_index=True)
    else:
        output_dataframe = pd.DataFrame(columns=["Header"] + REQUESTED_OUTPUT_COLUMNS)

    output_path = Path(args.out)
    output_dataframe.to_csv(output_path, index=False)
    print(f"Wrote {len(output_dataframe)} rows to {output_path}")


if __name__ == "__main__":
    main()
