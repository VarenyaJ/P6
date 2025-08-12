#!/usr/bin/env python3
"""
Search LOINC (local CSV) for codes/terms related to fetal_biometry headers.

Inputs:
  - LOINC CSV (e.g., /Users/varenya/Downloads/Loinc_2.80/LoincTableCore/LoincTableCore.csv)
  - Excel with a 'fetal_biometry' sheet (e.g., Sydney_Python_transformation.xlsx)

Output:
  - CSV of matches with the requested LOINC fields + the source Header.
"""

import argparse
import re
from pathlib import Path

import pandas as pd

REQUESTED_COLS = [
    "LOINC_NUM","COMPONENT","PROPERTY","TIME_ASPCT","SYSTEM",
    "SCALE_TYP","METHOD_TYP","CLASS","CLASSTYPE",
    "LONG_COMMON_NAME","SHORTNAME"
]

# Minimal helpful synonyms for common fetal biometry headers
ALIASES = {
    "bpd": ["bpd", "biparietal diameter"],
    "hc": ["hc", "head circumference"],
    "ac": ["ac", "abdominal circumference"],
    "efw": ["efw", "estimated fetal weight"],
    "fhr": ["fhr", "fetal heart rate"],
    "fl": ["fl", "femur length"],
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
    "occipitofrontal diameter": ["occipitofrontal diameter", "ofd"],
    "head": ["head"],
    "abdomen": ["abdomen", "abdominal"],
}

def normalize_header(h: str) -> str:
    if not isinstance(h, str):
        return ""
    # remove unit hints/parentheticals, keep core term
    base = re.sub(r"\(.*?\)", "", h).strip()
    # collapse multiple spaces & lowercase
    base = re.sub(r"\s+", " ", base).strip().lower()
    return base

def terms_for_header(header: str) -> list[str]:
    base = normalize_header(header)
    terms = set()
    if base:
        terms.add(base)
        # split tokens that might be useful alone (e.g., "Cisterna Magna")
        tokens = re.split(r"[\/\-\s]+", base)
        for t in tokens:
            t = t.strip()
            if t:
                terms.add(t)
        # add known aliases on base and on first token (common abbrev like "bpd")
        terms.update(ALIASES.get(base, []))
        if tokens:
            terms.update(ALIASES.get(tokens[0], []))
    # drop empties and return list
    return [t for t in terms if t]

def compile_pattern(terms: list[str]) -> str:
    # Build a regex that matches any of the terms (escaped)
    escaped = [re.escape(t) for t in terms if t]
    if not escaped:
        return r"^$"  # matches nothing
    return "(" + "|".join(escaped) + ")"

def find_matches_for_header(loinc_df: pd.DataFrame, header: str) -> pd.DataFrame:
    terms = terms_for_header(header)
    if not terms:
        return pd.DataFrame(columns=["Header"] + REQUESTED_COLS)

    pattern = compile_pattern(terms)

    fields_to_search = ["COMPONENT", "LONG_COMMON_NAME", "SHORTNAME"]
    mask = pd.Series(False, index=loinc_df.index)
    for field in fields_to_search:
        if field in loinc_df.columns:
            mask = mask | loinc_df[field].str.contains(pattern, case=False, na=False, regex=True)

    if not mask.any():
        return pd.DataFrame(columns=["Header"] + REQUESTED_COLS)

    matched = loinc_df.loc[mask, REQUESTED_COLS].copy()
    matched.insert(0, "Header", header)
    # drop exact dupes for cleanliness
    matched = matched.drop_duplicates()
    return matched

def main():
    parser = argparse.ArgumentParser(description="Search LOINC for fetal_biometry header matches.")
    parser.add_argument(
        "--loinc-csv",
        default="/Users/varenya/Downloads/Loinc_2.80/LoincTableCore/LoincTableCore.csv",
        help="Path to LoincTableCore.csv",
    )
    parser.add_argument(
        "--excel",
        default="Sydney_Python_transformation.xlsx",
        help="Path to Excel containing fetal_biometry sheet",
    )
    parser.add_argument(
        "--sheet",
        default="fetal_biometry",
        help="Sheet name that has the fetal biometry headers",
    )
    parser.add_argument(
        "--out",
        default="fetal_biometry_loinc_matches.csv",
        help="Output CSV for all matches",
    )
    args = parser.parse_args()

    loinc_csv = Path(args.loinc_csv)
    excel_path = Path(args.excel)

    if not loinc_csv.exists():
        raise FileNotFoundError(f"LOINC CSV not found: {loinc_csv}")
    if not excel_path.exists():
        raise FileNotFoundError(f"Excel file not found: {excel_path}")

    # Load LOINC
    loinc_df = pd.read_csv(
        loinc_csv,
        dtype=str,            # keep everything as strings
        low_memory=False
    ).fillna("")

    # Ensure requested columns exist
    missing_cols = [c for c in REQUESTED_COLS if c not in loinc_df.columns]
    if missing_cols:
        raise ValueError(f"Missing expected columns in LOINC CSV: {missing_cols}")

    # Load fetal_biometry headers
    xls = pd.ExcelFile(excel_path)
    if args.sheet not in xls.sheet_names:
        raise ValueError(f"Sheet '{args.sheet}' not found in {excel_path}. Found: {xls.sheet_names}")
    fetal_df = xls.parse(args.sheet)
    headers = list(fetal_df.columns)

    # Build matches
    all_matches = []
    for h in headers:
        m = find_matches_for_header(loinc_df, h)
        if not m.empty:
            all_matches.append(m)

    if all_matches:
        out_df = pd.concat(all_matches, ignore_index=True)
    else:
        out_df = pd.DataFrame(columns=["Header"] + REQUESTED_COLS)

    out_df.to_csv(args.out, index=False)
    print(f"Wrote {len(out_df)} rows to {args.out}")

if __name__ == "__main__":
    main()
