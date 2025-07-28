#!/usr/bin/env python3
"""
Command‑line interface for P6 toolkit.
Now determines sheet type by presence of multiple required key columns.
"""

import sys
import click
import pandas as pd
from P6.genotype import Genotype
from P6.phenotype import Phenotype

# Columns that need renaming → target dataclass fields
RENAME_MAP = {
    # genotype columns
    "ref":       "reference",
    "alt":       "alternate",
    "gene":      "gene_symbol",
    "start":     "start_position",
    "end":       "end_position",
    "chrom":     "chromosome",

    # phenotype columns
    "hpo":       "hpo_id",
    "timestamp": "date_of_observation",
}

# For any renamed field, the two neighbors it must sit between
EXPECTED_COLUMN_NEIGHBORS = {
    "start_position": ("chromosome", "end_position"),
    "end_position":   ("start_position", "reference"),
    "reference":      ("end_position", "alternate"),
    "alternate":      ("reference", "gene_symbol"),
    "gene_symbol":    ("alternate", "hgvsg"),
}

# Minimal required columns (after renaming) to identify each sheet type
GENOTYPE_KEY_COLUMNS = {
    "contact_email",
    "phasing",
    "chromosome",
    "start_position",
    "end_position",
    "reference",
    "alternate",
    "gene_symbol",
    "hgvsg",
    "hgvsc",
    "hgvsp",
    "zygosity",
    "inheritance",
}

PHENOTYPE_KEY_COLUMNS = {
    "hpo_id",
    "date_of_observation",
    "status",
}


def load_sheets_as_tables(workbook_path: str) -> dict[str, pd.DataFrame]:
    """
    Read each worksheet into a DataFrame:
      - first row = header
      - first column = index
      - normalize all headers to snake_case lowercase
      - apply renames from RENAME_MAP
    """
    excel = pd.ExcelFile(workbook_path, engine="openpyxl")
    tables: dict[str, pd.DataFrame] = {}

    for sheet_name in excel.sheet_names:
        df = pd.read_excel(
            excel,
            sheet_name=sheet_name,
            header=0,
            index_col=0,
            engine="openpyxl",
        )

        # CLEAN & NORMALIZE headers:
        df.columns = (
            df.columns
              .str.strip()
              .str.replace(r"\s*\(.*?\)", "", regex=True)  # drop any "(…)"
              .str.replace(r"\s+", "_",    regex=True)    # spaces → underscore
              .str.replace(":",         "",    regex=False)  # drop colons
              .str.lower()
        )

        # apply specific renames (e.g. "ref" → "reference")
        df = df.rename(
            columns={orig: target for orig, target in RENAME_MAP.items() if orig in df.columns}
        )

        tables[sheet_name] = df

    return tables


def verify_column_order(df: pd.DataFrame, field: str) -> None:
    """
    Ensure `field` sits between its expected neighbors.
    """
    before, after = EXPECTED_COLUMN_NEIGHBORS[field]
    cols = list(df.columns)
    idx = cols.index(field)
    if idx == 0 or idx == len(cols) - 1:
        raise ValueError(f"Field {field!r} at column index {idx} cannot satisfy neighbor check.")
    if cols[idx - 1] != before or cols[idx + 1] != after:
        raise ValueError(
            f"Field {field!r} should be between {before!r} and {after!r}, "
            f"found neighbors {cols[idx - 1]!r} and {cols[idx + 1]!r}."
        )


@click.group()
def main():
    """P6 CLI: parse genomic & phenotypic sheets into typed Python objects."""
    pass


@main.command(name="parse-excel")
@click.argument("excel_file", type=click.Path(exists=True))
def parse_excel(excel_file: str):
    """
    Read each sheet, check column order, then:
      - Identify as a Genotype sheet if ALL GENOTYPE_KEY_COLUMNS are present.
      - Identify as a Phenotype sheet if ALL PHENOTYPE_KEY_COLUMNS are present.
      - Otherwise skip.
    Instantiate objects accordingly.
    """
    all_sheets = load_sheets_as_tables(excel_file)
    genotype_records = []
    phenotype_records = []

    for sheet_name, df in all_sheets.items():
        columns = set(df.columns)
        is_genotype_sheet = GENOTYPE_KEY_COLUMNS.issubset(columns)
        is_phenotype_sheet = PHENOTYPE_KEY_COLUMNS.issubset(columns)

        if is_genotype_sheet == is_phenotype_sheet:
            click.echo(
                f"⚠  Skipping {sheet_name!r}: cannot unambiguously classify as genotype or phenotype",
                err=True,
            )
            continue

        id_column = "genotype_patient_ID" if is_genotype_sheet else "phenotype_patient_ID"
        working = df.reset_index().rename(columns={"index": id_column})

        # Verify ordering of any renamed columns
        for field in EXPECTED_COLUMN_NEIGHBORS:
            if field in working.columns:
                try:
                    verify_column_order(working, field)
                except ValueError as e:
                    click.echo(f"❌ Sheet {sheet_name!r}: {e}", err=True)
                    sys.exit(1)

        if is_genotype_sheet:
            for row_index, row_data in working.iterrows():
                genotype_records.append(Genotype(
                    genotype_patient_ID=str(row_data["genotype_patient_ID"]),
                    contact_email=row_data["contact_email"],
                    phasing=bool(row_data["phasing"]),
                    chromosome=row_data["chromosome"],
                    start_position=int(row_data["start_position"]),
                    end_position=int(row_data["end_position"]),
                    reference=row_data["reference"],
                    alternate=row_data["alternate"],
                    gene_symbol=row_data["gene_symbol"],
                    hgvsg=row_data["hgvsg"],
                    hgvsc=row_data["hgvsc"],
                    hgvsp=row_data["hgvsp"],
                    zygosity=row_data["zygosity"],
                    inheritance=row_data["inheritance"],
                ))
        else:
            for row_index, row_data in working.iterrows():
                phenotype_records.append(Phenotype(
                    phenotype_patient_ID=str(row_data["phenotype_patient_ID"]),
                    HPO_ID=row_data["hpo_id"],
                    date_of_observation=row_data["date_of_observation"],
                    status=bool(row_data["status"]),
                ))

    click.echo(f"Created {len(genotype_records)} Genotype objects")
    click.echo(f"Created {len(phenotype_records)} Phenotype objects")


if __name__ == "__main__":
    main()
