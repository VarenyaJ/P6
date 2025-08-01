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

from collections import defaultdict
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
@click.option(
    "-d",
    "--data-path",
    "data_dir",
    default="tests/data",
    type=click.Path(exists=True),
    help="where to save HPO JSON (default: tests/data)",
)
@click.option(
    "-v",
    "--hpo-version",
    default=None,
    type=str,
    help="exact HPO release tag (e.g. 2025-03-03 or v2025-03-03)",
)
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
@click.option(
    "-e",
    "--excel-path",
    "excel_file",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="path to the Excel workbook",
)
@click.option(
    "-hpo",
    "--custom-hpo",
    "hpo_path",
    type=click.Path(exists=True, dir_okay=False),
    help="path to a custom HPO JSON file (defaults to tests/data/hp.json)",
)
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
        hpo_file = pathlib.Path("tests/data") / "hp.json"
    if not hpo_file.is_file():
        click.echo(f"Error: HPO file not found at {hpo_file}", err=True)
        sys.exit(1)

    # load ontology & mapper
    ontology = hpotk.load_minimal_ontology(str(hpo_file))
    mapper = DefaultMapper(ontology)

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
    # Group by Patient
    # 1) Group all genotype and phenotype records by patient ID
    records_by_patient: dict[str, dict[str, list]] = defaultdict(
        lambda: {"genotypes": [], "phenotypes": []}
    )
    for genotype in genotype_records:
        patient_id = genotype.genotype_patient_ID
        records_by_patient[patient_id]["genotypes"].append(genotype)

    for phenotype in phenotype_records:
        patient_id = phenotype.phenotype_patient_ID
        records_by_patient[patient_id]["phenotypes"].append(phenotype)

    # 2) Prepare a timestamped output directory
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    output_dir = pathlib.Path.cwd() / "phenopacket-from-excel" / timestamp / "phenopackets"
    output_dir.mkdir(parents=True, exist_ok=True)

    # 3) Build and write one Phenopacket per patient
    for patient_id, patient_data in records_by_patient.items():
        phenopacket = Phenopacket()
        phenopacket.id = patient_id
        phenopacket.subject.id = patient_id

        # TODO: revisit with Daniel to also set feature.onset and feature.resolution using TimeElement

        # 3a) Add HPO-based phenotypic features
        for phenotype in patient_data["phenotypes"]:
            feature = phenopacket.phenotypic_features.add()
            feature.type.id = phenotype.HPO_ID
            # mark as excluded if status is False
            if not phenotype.status:
                feature.excluded = True

        # 3b) Add variant interpretations from genotype records.
        # Genotypes → Interpretation → Diagnosis → GenomicInterpretation
        for interpretation_index, genotype_record in enumerate(patient_data["genotypes"]):
            # Create a new Interpretation entry (must set id and progress_status)
            interpretation = phenopacket.interpretations.add()
            interpretation.id = f"{patient_id}-interp-{interpretation_index}"  # 0-based index
            interpretation.progress_status = interpretation.ProgressStatus.COMPLETED  # correct enum

            # each Interpretation has a `diagnosis` submessage (no `id` field)
            diagnosis = interpretation.diagnosis


            # inside Diagnosis, add GenomicInterpretation
            genomic_interpretation_entry = diagnosis.genomic_interpretations.add()
            genomic_interpretation_entry.subject_or_biosample_id = patient_id
            genomic_interpretation_entry.interpretation_status = genomic_interpretation_entry.InterpretationStatus.CONTRIBUTORY

            # now fill in the VariationDescriptor
            # TODO: set this up later
            # Omit setting gene_context for now.
            # variation_descriptor = genomic_interpretation_entry.variant_interpretation.variation_descriptor
            # we can also set variation_descriptor.gene_context and variation_descriptor.allelic_state here then serialize out as before
            # variation_descriptor.gene_context.gene_symbol = genotype_record.gene_symbol
            # variation_descriptor.allelic_state = variation_descriptor.AllelicState.Value(genotype_record.zygosity.upper())

            # TODO: when ready, add an Expression.HGVS here
            # Record the HGVS genomic notation as an Expression
            # expr = variation_descriptor.expressions.add()
            # expr.syntax = Phenopacket.Diagnosis.GenomicInterpretation.VariantInterpretation.VariationDescriptor.Expression.HGVS
            # expr.value = genotype_record.hgvsg

        # 3d) Serialize to a single JSON file per patient
        out_path = output_dir / f"{patient_id}.json"
        with open(out_path, "w", encoding="utf-8") as out_f:
            out_f.write(MessageToJson(phenopacket))

    click.echo(f"Wrote {len(records_by_patient)} phenopacket files to {output_dir}")
    click.echo(f"Created {len(genotype_records)} Genotype objects")
    click.echo(f"Created {len(phenotype_records)} Phenotype objects")


if __name__ == "__main__":
    main()
