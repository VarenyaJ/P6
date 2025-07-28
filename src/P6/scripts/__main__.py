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
    "ref": "reference",
    "alt": "alternate",
    "gene": "gene_symbol",
    "start": "start_position",
    "end": "end_position",
}

# For any renamed field, the two neighbors it must sit between
EXPECTED_COLUMN_NEIGHBORS = {
    "start_position": ("chromosome", "end_position"),
    "end_position": ("start_position", "reference"),
    "reference": ("end_position", "alternate"),
    "alternate": ("reference", "gene_symbol"),
    "gene_symbol": ("alternate", "hgvsg"),
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
    "HPO_ID",
    "date_of_observation",
    "status",
}


def load_sheets_as_tables(workbook_path: str) -> dict[str, pd.DataFrame]:
    """
    Read each worksheet into a DataFrame:
      - first row = header
      - first column = index
      - apply renames from RENAME_MAP
    """
    excel = pd.ExcelFile(workbook_path, engine="openpyxl")
    tables = {}
    for sheet_name in excel.sheet_names:
        df = pd.read_excel(
            excel,
            sheet_name=sheet_name,
            header=0,
            index_col=0,
            engine="openpyxl",
        )
        df = df.rename(columns={k: v for k, v in RENAME_MAP.items() if k in df.columns})
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
            f"found neighbors {cols[idx-1]!r} and {cols[idx+1]!r}."
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
        # Determine sheet type by the full set of required columns
        columns = set(df.columns)
        is_genotype_sheet = GENOTYPE_KEY_COLUMNS.issubset(columns)
        is_phenotype_sheet = PHENOTYPE_KEY_COLUMNS.issubset(columns)

        if is_genotype_sheet == is_phenotype_sheet:
            # Either both True (unlikely) or both False → ambiguous
            click.echo(f"⚠️  Skipping {sheet_name!r}: cannot unambiguously classify as genotype or phenotype", err=True)
            continue

        # Reset index into the appropriate patient_ID field
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

        # Build domain objects
        if is_genotype_sheet:
            for _, row in working.iterrows():
                genotype_records.append(Genotype(
                    genotype_patient_ID=str(row["genotype_patient_ID"]),
                    contact_email=row["contact_email"],
                    phasing=bool(row["phasing"]),
                    chromosome=row["chromosome"],
                    start_position=int(row["start_position"]),
                    end_position=int(row["end_position"]),
                    reference=row["reference"],
                    alternate=row["alternate"],
                    gene_symbol=row["gene_symbol"],
                    hgvsg=row["hgvsg"],
                    hgvsc=row["hgvsc"],
                    hgvsp=row["hgvsp"],
                    zygosity=row["zygosity"],
                    inheritance=row["inheritance"],
                ))
        else:
            for _, row in working.iterrows():
                phenotype_records.append(Phenotype(
                    phenotype_patient_ID=str(row["phenotype_patient_ID"]),
                    HPO_ID=row["HPO_ID"],
                    date_of_observation=row["date_of_observation"],
                    status=bool(row["status"]),
                ))

    click.echo(f"Created {len(genotype_records)} Genotype objects")
    click.echo(f"Created {len(phenotype_records)} Phenotype objects")


if __name__ == "__main__":
    main()
