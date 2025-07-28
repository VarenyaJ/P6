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

# Map raw zygosity abbreviations to allowed dataclass zygosity values
ZYGOSITY_MAP = {
    "het":     "heterozygous",
    "hom":     "homozygous",
    "comphet": "compound_heterozygosity",
    "hemi":    "hemizygous",
    "mosaic":  "mosaic",
}

# Map raw inheritance abbreviations to allowed dataclass inheritance values
INHERITANCE_MAP = {
    "unknown":   "unknown",
    "inherited": "inherited",
    "denovo":    "de_novo_mutation",
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

        # rename the former-index column (whatever its original header was)
        id_column = "genotype_patient_ID" if is_genotype_sheet else "phenotype_patient_ID"
        working = df.reset_index()
        original_index_col = working.columns[0]
        working = working.rename(columns={original_index_col: id_column})

        # Verify ordering of any renamed columns
        for field in EXPECTED_COLUMN_NEIGHBORS:
            if field in working.columns:
                try:
                    verify_column_order(working, field)
                except ValueError as e:
                    click.echo(f"❌ Sheet {sheet_name!r}: {e}", err=True)
                    sys.exit(1)

        if is_genotype_sheet:
            for _, row in working.iterrows():
                # handle slash‑separated zygosity and inheritance
                zyg_list = [z.strip().lower() for z in str(row["zygosity"]).split("/")]
                inh_list = [i.strip().lower() for i in str(row["inheritance"]).split("/")]

                for z_code, i_code in zip(zyg_list, inh_list):
                    if z_code not in ZYGOSITY_MAP:
                        click.echo(f"❌ Sheet {sheet_name!r}: Unrecognized zygosity code {z_code!r}", err=True)
                        sys.exit(1)
                    if i_code not in INHERITANCE_MAP:
                        click.echo(f"❌ Sheet {sheet_name!r}: Unrecognized inheritance code {i_code!r}", err=True)
                        sys.exit(1)

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
                        zygosity=ZYGOSITY_MAP[z_code],
                        inheritance=INHERITANCE_MAP[i_code],
                    ))
        else:
            for _, row in working.iterrows():
                # --- normalize phenotype fields into valid strings ---
                raw_hpo = row["hpo_id"]
                hpo_str = str(raw_hpo).strip()
                if hpo_str.lower().startswith("hp:"):
                    digits = hpo_str[3:]
                    hpo_id = f"HP:{digits.zfill(7)}"
                else:
                    hpo_id = hpo_str.zfill(7)

                raw_date = row["date_of_observation"]
                date_str = str(raw_date).strip()
                if not date_str.upper().startswith("T"):
                    date_str = f"T{date_str}"

                phenotype_records.append(Phenotype(
                    phenotype_patient_ID=str(row["phenotype_patient_ID"]),
                    HPO_ID=hpo_id,
                    date_of_observation=date_str,
                    status=bool(row["status"]),
                ))

    click.echo(f"Created {len(genotype_records)} Genotype objects")
    click.echo(f"Created {len(phenotype_records)} Phenotype objects")


if __name__ == "__main__":
    main()
