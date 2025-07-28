#!/usr/bin/env python3
"""
P6 CLI: parse one or more genomic/phenotypic Excel workbooks into objects,
with optional verbose logging and timestamped logfile output.
"""

import sys
import logging
import click
import pandas as pd
from P6.genotype import Genotype
from P6.phenotype import Phenotype

# --- Constants --------------------------------------------------------------

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

EXPECTED_COLUMN_NEIGHBORS = {
    "start_position": ("chromosome", "end_position"),
    "end_position":   ("start_position", "reference"),
    "reference":      ("end_position", "alternate"),
    "alternate":      ("reference", "gene_symbol"),
    "gene_symbol":    ("alternate", "hgvsg"),
}

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

ZYGOSITY_MAP = {
    "het":     "heterozygous",
    "hom":     "homozygous",
    "comphet": "compound_heterozygosity",
    "hemi":    "hemizygous",
    "mosaic":  "mosaic",
}

INHERITANCE_MAP = {
    "unknown":   "unknown",
    "inherited": "inherited",
    "denovo":    "de_novo_mutation",
}


def load_sheets_as_tables(workbook_path: str) -> dict[str, pd.DataFrame]:
    """
    Load each sheet into a DataFrame, normalize headers and apply RENAME_MAP.
    """
    excel = pd.ExcelFile(workbook_path, engine="openpyxl")
    sheet_tables: dict[str, pd.DataFrame] = {}

    for sheet_name in excel.sheet_names:
        df = pd.read_excel(
            excel,
            sheet_name=sheet_name,
            header=0,
            index_col=0,
            engine="openpyxl",
        )
        df.columns = (
            df.columns
              .str.strip()
              .str.replace(r"\s*\(.*?\)", "", regex=True)
              .str.replace(r"\s+", "_",          regex=True)
              .str.replace(":",           "",   regex=False)
              .str.lower()
        )
        df = df.rename(
            columns={orig: new for orig, new in RENAME_MAP.items() if orig in df.columns}
        )
        sheet_tables[sheet_name] = df

    return sheet_tables


def verify_column_order(df: pd.DataFrame, field: str) -> None:
    """
    Ensure that `field` appears directly between its expected neighbors.
    """
    before, after = EXPECTED_COLUMN_NEIGHBORS[field]
    cols = list(df.columns)
    idx = cols.index(field)
    if idx == 0 or idx == len(cols) - 1:
        raise ValueError(f"Field {field!r} at index {idx} cannot satisfy neighbor check")
    if cols[idx - 1] != before or cols[idx + 1] != after:
        raise ValueError(
            f"Field {field!r} should sit between {before!r} and {after!r}, "
            f"found {cols[idx-1]!r} / {cols[idx+1]!r}"
        )


@click.group()
def main():
    """Entry point for the P6 CLI."""
    pass


@main.command(name="parse-excel")
@click.argument(
    "workbook_paths",
    nargs=-1,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--verbose-logging",
    is_flag=True,
    help="Also emit debug logs to stderr",
)
@click.option(
    "--log-file-path",
    type=click.Path(dir_okay=False, writable=True),
    help="Append timestamped logs to this file",
)
def parse_excel(
    workbook_paths: tuple[str, ...],
    verbose_logging: bool,
    log_file_path: str,
):
    """
    Parse one or more Excel workbooks; classify each sheet as GENOTYPE or PHENOTYPE,
    instantiate the corresponding objects, and summarize counts.
    """
    if not workbook_paths:
        click.echo("❌  No input files specified.", err=True)
        sys.exit(1)

    # — configure logging —
    handlers: list[logging.Handler] = []
    if log_file_path:
        handlers.append(logging.FileHandler(log_file_path, mode="a", encoding="utf-8"))
    if verbose_logging:
        handlers.append(logging.StreamHandler(sys.stderr))
    if handlers:
        logging.basicConfig(
            level=logging.DEBUG if verbose_logging else logging.INFO,
            format="%(asctime)s %(levelname)s %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            handlers=handlers,
        )

    for workbook_path in workbook_paths:
        logging.info(f"Beginning parse of '{workbook_path}'")
        try:
            sheet_tables = load_sheets_as_tables(workbook_path)
        except Exception as e:
            logging.error(f"Failed to read '{workbook_path}': {e}")
            continue

        logging.debug(f"Loaded sheets: {list(sheet_tables.keys())}")
        count_genotype_records = 0
        count_phenotype_records = 0

        for sheet_name, sheet_df in sheet_tables.items():
            cols = set(sheet_df.columns)
            is_genotype = GENOTYPE_KEY_COLUMNS.issubset(cols)
            is_phenotype = PHENOTYPE_KEY_COLUMNS.issubset(cols)

            if is_genotype == is_phenotype:
                logging.warning(f"Skipping sheet {sheet_name!r}: ambiguous classification")
                continue

            working_df = sheet_df.reset_index()
            orig_index = working_df.columns[0]
            new_index_col = (
                "genotype_patient_ID" if is_genotype else "phenotype_patient_ID"
            )
            working_df = working_df.rename(columns={orig_index: new_index_col})

            # verify any renamed columns are in correct order
            for fld in EXPECTED_COLUMN_NEIGHBORS:
                if fld in working_df.columns:
                    try:
                        verify_column_order(working_df, fld)
                    except ValueError as err:
                        logging.error(f"Sheet {sheet_name!r}: {err}")
                        sys.exit(1)

            if is_genotype:
                logging.debug(f"Sheet {sheet_name!r} → GENOTYPE")
                for _, row in working_df.iterrows():
                    zyg_codes = [z.strip().lower() for z in str(row["zygosity"]).split("/")]
                    inh_codes = [i.strip().lower() for i in str(row["inheritance"]).split("/")]
                    for zc, ic in zip(zyg_codes, inh_codes):
                        if zc not in ZYGOSITY_MAP:
                            logging.error(f"Unknown zygosity code {zc!r} in {sheet_name!r}")
                            sys.exit(1)
                        if ic not in INHERITANCE_MAP:
                            logging.error(f"Unknown inheritance code {ic!r} in {sheet_name!r}")
                            sys.exit(1)

                        Genotype(
                            genotype_patient_ID=str(row[new_index_col]),
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
                            zygosity=ZYGOSITY_MAP[zc],
                            inheritance=INHERITANCE_MAP[ic],
                        )
                        count_genotype_records += 1

            else:
                logging.debug(f"Sheet {sheet_name!r} → PHENOTYPE")
                for _, row in working_df.iterrows():
                    Phenotype(
                        phenotype_patient_ID=str(row[new_index_col]),
                        HPO_ID=row["hpo_id"],
                        date_of_observation=row["date_of_observation"],
                        status=bool(row["status"]),
                    )
                    count_phenotype_records += 1

        click.echo(f"{workbook_path}: Created {count_genotype_records} Genotype objects")
        click.echo(f"{workbook_path}: Created {count_phenotype_records} Phenotype objects")


if __name__ == "__main__":
    main()
