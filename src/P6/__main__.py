"""
Command‑line interface for P6 toolkit.
Now determines sheet type by presence of multiple required key columns,
and normalizes numeric HPO IDs and timestamps.
"""

import pathlib
import sys
import typing
import click

import hpotk

from stairval.notepad import create_notepad

from .loader import load_sheets_as_tables
from .mapper import DefaultMapper



@click.group()
def main():
    """P6: Peter's Parse and Processing of Prenatal Particulars via Pandas."""
    pass


@main.command(name="download")
@click.option("--d", default="data", type=click.Path(exists=True))
@click.option("--hpo-version", default=None, type=typing.Optional[str])
def download(d: str, hpo_version: typing.Optional[str],):
    # TODO: download an HPO
    pass

@main.command(name="parse-excel")
@click.argument("excel_file", type=click.Path(exists=True))
@click.option("--d", default="data", type=click.Path(exists=True))
def parse_excel(excel_file: str, d: str):
    """
    Read each sheet, check column order, then:
      - Identify as a Genotype sheet if ALL GENOTYPE_KEY_COLUMNS are present.
      - Identify as a Phenotype sheet if ALL PHENOTYPE_KEY_COLUMNS are present.
      - Otherwise skip.
    Instantiate objects accordingly, normalizing HPO IDs & timestamps.
    """
    datadir = pathlib.Path(d)
    assert datadir.exists() and datadir.is_dir(), "Data directory must exist"

    fpath_hpo = datadir.joinpath("hp.json")
    assert fpath_hpo.is_file(), "HPO file must exist"
    hpo = hpotk.load_minimal_ontology(str(fpath_hpo))
    mapper = DefaultMapper(hpo)
    
    all_sheets = load_sheets_as_tables(excel_file)
    notepad = create_notepad("phenopackets")
    pps = mapper.apply_mapping(all_sheets, notepad)
    
    assert not notepad.has_errors_or_warnings(include_subsections=True)
    # TODO: write phenopackets to a folder


if __name__ == "__main__":
    main()
