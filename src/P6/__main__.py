"""
Command‑line interface for P6 toolkit.
Now determines sheet type by presence of multiple required key columns,
and normalizes numeric HPO IDs and timestamps.
"""
import click
import hpotk
import pathlib
import requests
import sys
import typing

from datetime import datetime
from google.protobuf.json_format import MessageToJson
from stairval.notepad import create_notepad
from phenopackets.schema.v2.phenopackets_pb2 import Phenopacket

from .loader import load_sheets_as_tables
from .mapper import DefaultMapper


@click.group()
def main():
    """P6: Peter's Parse and Processing of Prenatal Particulars via Pandas."""
    pass


@main.command(name="download")
@click.option("-d", "--data-path", "data_dir", default="tests/data", type=click.Path(exists=True), help="where to save HPO JSON (default: tests/data)")
@click.option("-v", "--hpo-version", default=None, type=str, help="exact HPO release tag (e.g. 2025-03-03 or v2025-03-03)")
def download(data_dir: str, hpo_version: typing.Optional[str]):
    # TODO: download an HPO
    """
    Download a specific or the latest HPO JSON release into the tests/data/ folder.
    """
    datadir = pathlib.Path(data_dir)
    datadir.mkdir(parents=True, exist_ok=True)
    # figure out which tag to download
    if hpo_version:
        tag = hpo_version if hpo_version.startswith("v") else f"v{hpo_version}"
    else:
        # GitHub latest‐release API
        resp = requests.get(
            "https://api.github.com/repos/obophenotype/human-phenotype-ontology/releases/latest"
        )
        resp.raise_for_status()
        tag = resp.json()["tag_name"]
    url = (
        f"https://github.com/obophenotype/human-phenotype-ontology/"
        f"releases/download/{tag}/hp.json"
    )
    click.echo(f"Downloading HPO release {tag} …")
    resp = requests.get(url)
    resp.raise_for_status()

    out = datadir / "hp.json"
    with open(out, "wb") as f:
        f.write(resp.content)

    click.echo(f"Saved HPO JSON to {out}")
    pass


@main.command(name="parse-excel")
@click.option("-e", "--excel-path", "excel_file", required=True, type=click.Path(exists=True, dir_okay=False), help="path to the Excel workbook")
@click.option("-hpo", "--custom-hpo", "hpo_path", type=click.Path(exists=True, dir_okay=False), help="path to a custom HPO JSON file (defaults to tests/data/hp.json)")
def parse_excel(excel_file: str, hpo_path: typing.Optional[str] = None):
    """
    Read each sheet, check column order, then:
      - Identify as a Genotype sheet if ALL GENOTYPE_KEY_COLUMNS are present.
      - Identify as a Phenotype sheet if ALL PHENOTYPE_KEY_COLUMNS are present.
      - Otherwise skip.
    Instantiate objects accordingly, normalizing HPO IDs & timestamps.
    """
    # pick HPO JSON: either custom or default
    if hpo_path:
        hpo_file = pathlib.Path(hpo_path)
    else:
        hpo_file = hpo_file = pathlib.Path("tests/data") / "hp.json"
    if not hpo_file.is_file():
        click.echo(f"Error: HPO file not found at {hpo_file}", err=True)
        sys.exit(1)

    # load ontology & mapper
    hpo = hpotk.load_minimal_ontology(str(hpo_file))
    mapper = DefaultMapper(hpo)

    # read all sheets
    all_sheets = load_sheets_as_tables(excel_file)
    notepad = create_notepad("phenopackets")

    # get back two separate lists, map to domain records
    genotype_records, phenotype_records = mapper.apply_mapping(all_sheets, notepad)

    # if there were errors, show them and exit non‑zero
    if notepad.has_errors(include_subsections=True):
        click.echo("Errors found in mapping:")
        for err in notepad.errors():
            click.echo(f"- {err}")
        # sys.exit(1)
        # no exit—always continue to print the counts

    # show any warnings but keep going
    if notepad.has_warnings(include_subsections=True):
        click.echo("Warnings found in mapping:")
        for w in notepad.warnings():
            click.echo(f"- {w}")

    # pps = mapper.apply_mapping(all_sheets, notepad)
    # assert not notepad.has_errors_or_warnings(include_subsections=True)
    # TODO: write phenopackets to a folder
    # click.echo(f"Created {len(pps)} Phenotype objects")

    # write genotype and phenotype records out as JSON
    # use YYYY-MM-DD_HH-MM-SS for human-readable timestamps
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    base = pathlib.Path.cwd() / "phenopacket-from-excel" / timestamp / "phenopackets"
    base.mkdir(parents=True, exist_ok=True)

    total = 0
    for rec in genotype_records + phenotype_records:
        # one protobuffer per record
        pkt = Phenopacket()
        # every phenopacket needs at least an ID
        if hasattr(rec, "genotype_patient_ID"):
            pkt.id = f"genotype-{rec.genotype_patient_ID}"
        else:
            pkt.id = f"phenotype-{rec.phenotype_patient_ID}"
        # TODO: Look if we can populate other fields here, e.g. observations, units, etc.

        fn = base / f"{pkt.id}.json"
        with open(fn, "w", encoding="utf-8") as out_f:
            # serialize to JSON text
            out_f.write(MessageToJson(pkt))
        total += 1
    click.echo(f"Wrote {total} phenopacket files to {base}")
    click.echo(f"Created {len(genotype_records)} Genotype objects")
    click.echo(f"Created {len(phenotype_records)} Phenotype objects")


if __name__ == "__main__":
    main()
