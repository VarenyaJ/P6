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



def collect_search_terms_for_header(header_text: str) -> List[str]:
    """
    Build a list of search terms for a given header.

    The list includes:
    - The normalized base term.
    - Tokenized words from the base term.
    - Helpful aliases for the base term and its first token.

    Parameters
    ----------
    header_text : str
        The column header from the Excel sheet.

    Returns
    -------
    list of str
        Unique, non-empty search terms, sorted for determinism.
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
    Construct a word-bounded, non-capturing regex that matches any of the terms.

    Parameters
    ----------
    terms : Iterable[str]
        The search terms.

    Returns
    -------
    str
        Regex string suitable for ``Series.str.contains(..., regex=True)``.
    """
    escaped_terms = [re.escape(t) for t in terms if t]
    if not escaped_terms:
        # A regex that never matches (used when no terms were generated)
        return r"(?!x)x"
    return r"\b(?:%s)\b" % "|".join(escaped_terms)



def row_has_ob_fetal_context(row: pd.Series) -> bool:
    """
    Heuristic: decide whether a LOINC row indicates OB/fetal ultrasound context.

    Signals (any is sufficient)
    ---------------------------
    - CLASS is OB.US or PANEL.OB.US
    - SYSTEM contains fetal/placental structures or related anatomy
    - METHOD_TYP contains 'US'
    - Text contains 'fetal'/'fetus' in COMPONENT/LONG_COMMON_NAME/SHORTNAME

    Parameters
    ----------
    row : pandas.Series
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
    row : pandas.Series
        The LOINC row being scored.
    header_terms_lower : set of str
        Header-derived terms, all lowercased.
    allowed_properties : set of str
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
        if (
            f" {term} " in f" {component_text} "
            or f" {term} " in f" {long_name_text} "
            or f" {term} " in f" {short_name_text} "
        ):
            score += 2

    # Property hint (if applicable)
    if allowed_properties and property_value in allowed_properties:
        score += 1

    return score



def get_property_hints_for_header(header: str) -> Set[str]:
    """
    Retrieve PROPERTY hints for a header using both the full normalized header
    and its first token (if present).

    Parameters
    ----------
    header : str
        Header text from the Excel sheet.

    Returns
    -------
    set of str
        Allowed PROPERTY values for this header (may be empty).
    """
    normalized_header = normalize_header_to_term(header)
    hint_keys: Set[str] = {normalized_header} if normalized_header else set()
    if normalized_header:
        hint_keys.add(normalized_header.split()[0])
    allowed: Set[str] = set()
    for key in hint_keys:
        allowed |= PROPERTY_HINTS_BY_HEADER.get(key, set())
    return allowed



def build_search_mask(
    loinc_dataframe: pd.DataFrame,
    regex: str,
    fields_to_search: Sequence[str],
) -> pd.Series:
    """
    Build a boolean mask over the LOINC dataframe for the given regex and fields.

    Parameters
    ----------
    loinc_dataframe : pandas.DataFrame
        LOINC table.
    regex : str
        Compiled regex string to search for.
    fields_to_search : sequence of str
        Column names to search within.

    Returns
    -------
    pandas.Series
        Boolean mask aligned to ``loinc_dataframe.index``.
    """
    mask = pd.Series(False, index=loinc_dataframe.index)
    for field in fields_to_search:
        mask = mask | loinc_dataframe[field].str.contains(regex, case=False, na=False, regex=True)
    return mask




def process_single_header(
    loinc_dataframe: pd.DataFrame,
    header: str,
    fields_to_search: Sequence[str],
    max_per_header: int,
    use_context_filter: bool,
) -> pd.DataFrame:
    """
    Process one header: generate terms, search, context-filter, score, and trim.

    Parameters
    ----------
    loinc_dataframe : pandas.DataFrame
        The full LOINC table.
    header : str
        Header from the Excel sheet.
    fields_to_search : sequence of str
        LOINC columns to search (e.g., COMPONENT, LONG_COMMON_NAME, SHORTNAME).
    max_per_header : int
        Maximum rows to keep after scoring (``<=0`` means keep all).
    use_context_filter : bool
        If True, apply OB/fetal context filtering.

    Returns
    -------
    pandas.DataFrame
        Result rows for this header with the requested output columns, plus
        a leading ``Header`` column. May be empty.
    """
    search_terms = collect_search_terms_for_header(header)
    if not search_terms:
        LOGGER.debug("[DEBUG] Header '%s': no usable search terms after normalization.", header)
        return pd.DataFrame(columns=["Header"] + REQUESTED_OUTPUT_COLUMNS)

    regex = build_search_regex(search_terms)
    search_mask = build_search_mask(loinc_dataframe, regex, fields_to_search)
    candidate_rows = loinc_dataframe.loc[search_mask].copy()

    if candidate_rows.empty:
        LOGGER.debug("[DEBUG] Header '%s': 0 raw matches.", header)
        return pd.DataFrame(columns=["Header"] + REQUESTED_OUTPUT_COLUMNS)

    raw_count = len(candidate_rows)

    if use_context_filter:
        context_mask = candidate_rows.apply(row_has_ob_fetal_context, axis=1)
        candidate_rows = candidate_rows.loc[context_mask]

    filtered_count = len(candidate_rows)
    if candidate_rows.empty:
        LOGGER.debug(
            "[DEBUG] Header '%s': 0 matches after OB/fetal context filter.", header
        )
        return pd.DataFrame(columns=["Header"] + REQUESTED_OUTPUT_COLUMNS)

    allowed_properties_for_header = get_property_hints_for_header(header)

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

    if max_per_header and max_per_header > 0:
        candidate_rows = candidate_rows.head(max_per_header)

    kept_count = len(candidate_rows)

    # Finalize columns and add the originating header
    candidate_rows = candidate_rows.drop(columns=["__score"], errors="ignore")
    candidate_rows.insert(0, "Header", header)
    candidate_rows = candidate_rows[["Header"] + REQUESTED_OUTPUT_COLUMNS].drop_duplicates()

    LOGGER.debug(
        "[DEBUG] Header '%s': raw=%d, after_context=%d, kept=%d, allowed_props=%s",
        header,
        raw_count,
        filtered_count,
        kept_count,
        sorted(allowed_properties_for_header) if allowed_properties_for_header else "none",
    )

    return candidate_rows




def write_output(output_dataframe: pd.DataFrame, output_path: Path) -> None:
    """
    Write results to CSV and log a short summary.

    Parameters
    ----------
    output_dataframe : pandas.DataFrame
        Final concatenated results with a leading ``Header`` column.
    output_path : Path
        Destination CSV path.
    """
    output_dataframe.to_csv(output_path, index=False)
    LOGGER.info("Wrote %d rows to %s", len(output_dataframe), output_path)

