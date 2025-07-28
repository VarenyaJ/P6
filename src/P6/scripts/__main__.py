"""
Command‑line interface for P6 toolkit.
"""

import click
from openpyxl import load_workbook
import pandas as pd

from P6.genotype import Genotype
from P6.phenotype import Phenotype


def load_named_tables(workbook_path: str) -> dict[str, pd.DataFrame]:
    """
    Discover every named table in an Excel workbook and return
    a mapping of table_name -> DataFrame whose columns match
    exactly the attribute names of your domain classes.
    """
    workbook = load_workbook(filename=workbook_path, data_only=True)
    tables: dict[str, pd.DataFrame] = {}

    for worksheet in workbook.worksheets:
        for table in worksheet.tables.values():
            df = pd.read_excel(
                workbook_path,
                sheet_name=worksheet.title,
                header=0,
                index_col=0,
                usecols=table.ref,
                engine="openpyxl",
            )

            # rename only where Excel headers don't already match our attributes
            df = df.rename(columns={
                "ref": "reference",
                "alt": "alternate",
                "gene": "gene_symbol",
                "start": "start_position",
                "end": "end_position",
            })

            tables[table.name] = df

    return tables


@click.group()
def main():
    """P6 CLI: parse genomic & phenotypic tables into Python objects."""
    pass


@main.command()
@click.argument("excel_file", type=click.Path(exists=True))
def parse_excel(excel_file: str):
    """
    Parse every named table in the Excel file, build Genotype and
    Phenotype objects, and report counts.
    """
    named_tables = load_named_tables(excel_file)

    genotype_records: list[Genotype] = []
    phenotype_records: list[Phenotype] = []

    for table_name, df in named_tables.items():
        # Turn the index into the patient‐ID column
        df = df.reset_index().rename(columns={"index": "patient_id"})

        if "phasing" in df.columns:
            # ── Genotype table ──
            for _, row in df.iterrows():
                genotype_records.append(
                    Genotype(
                        genotype_patient_ID=str(row["patient_id"]),
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
                    )
                )

        elif "HPO_ID" in df.columns:
            # ── Phenotype table ──
            for _, row in df.iterrows():
                phenotype_records.append(
                    Phenotype(
                        phenotype_patient_ID=str(row["patient_id"]),
                        HPO_ID=row["HPO_ID"],
                        date_of_observation=row["date_of_observation"],
                        status=bool(row["status"]),
                    )
                )

        else:
            click.echo(f"Skipping table {table_name!r}: unknown schema")

    click.echo(f"Created {len(genotype_records)} Genotype objects")
    click.echo(f"Created {len(phenotype_records)} Phenotype objects")


if __name__ == "__main__":
    main()
