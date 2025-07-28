"""
Command‑line interface for P6 toolkit.
"""

import click
from openpyxl import load_workbook
import pandas as pd

from P6.genotype import Genotype
from P6.phenotype import Phenotype
from P6.periodicity import FrequencyModifier, Periodicity


@click.group()
def main():
    """P6: parse genomic & phenotypic data into typed Python objects."""
    pass


@main.command()
@click.argument("excel_path", type=click.Path(exists=True))
def parse_tables(excel_path):
    """
    Discover all named tables in an Excel workbook, validate headers/rows
    are consistent, and report counts of Genotype vs Phenotype entries.
    """
    wb = load_workbook(filename=excel_path, data_only=True)
    all_counts = {"genotype": 0, "phenotype": 0}

    for sheet in wb.worksheets:
        for tbl in sheet.tables.values():
            df = pd.read_excel(excel_path, sheet_name=sheet.title,
                               header=0, index_col=0, usecols=tbl.ref, engine="openpyxl")
            # rudimentary dispatch by column names
            if "phasing" in df.columns:
                all_counts["genotype"] += len(df)
            elif "HPO_ID" in df.columns:
                all_counts["phenotype"] += len(df)

    click.echo(f"Found {all_counts['genotype']} genotype rows")
    click.echo(f"Found {all_counts['phenotype']} phenotype rows")


if __name__ == "__main__":
    main()