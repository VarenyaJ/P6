"""
Command‑line interface for P6 toolkit.
Now determines sheet type by presence of multiple required key columns,
and normalizes numeric HPO IDs and timestamps.
"""

import click
import hpotk
import pandas as pd  # Not needed for Pandas_Workaround, i.e. don't call declare or call "_read_sheets" at all, just use `tables = load_sheets_as_tables(excel_file)` which only needs `from .loader import load_sheets_as_tables`
import pathlib
import requests
import sys
import typing

from collections import defaultdict, namedtuple
from datetime import datetime
from google.protobuf.json_format import MessageToJson
from stairval.notepad import create_notepad
from phenopackets.schema.v2.phenopackets_pb2 import Phenopacket
import phenopackets.schema.v2 as pps2

from .loader import load_sheets_as_tables
from .mapper import DefaultMapper

AuditEntry = namedtuple("AuditEntry", ["step", "sheet", "message", "level"])

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
@click.option("--strict-variants/--no-strict-variants", default=False, help=("Treat raw↔HGVS mismatches as errors (default: warn)."))
@click.option("--verbose", is_flag=True, help="Show preprocessing and classification steps")
def parse_excel(excel_file: str, hpo_path: typing.Optional[str] = None, verbose: bool = False, strict_variants: bool = False):
    """
    Read each sheet, check column order, then:
      - Identify as a Genotype sheet if ALL GENOTYPE_KEY_COLUMNS are present.
      - Identify as a Phenotype sheet if ALL PHENOTYPE_KEY_COLUMNS are present.
      - Otherwise skip.
    Instantiate objects accordingly, normalizing HPO IDs & timestamps.
    """
    # 1) Load (or locate) the HPO JSON file
    hpo_file = _locate_hpo_file(hpo_path)

    # 2) Build ontology and mapper
    ontology = _load_ontology(str(hpo_file))
    mapper = DefaultMapper(ontology, strict_variants=strict_variants)

    # 3) Read all sheets into DataFrames
    tables = _read_sheets(excel_file)
    # tables = load_sheets_as_tables(excel_file)  # Just use this for Pandas_Workaround. Don't call declare or call "_read_sheets" at all. Just use `tables = load_sheets_as_tables(excel_file)` which only needs `from .loader import load_sheets_as_tables`
    # TODO: Decide if it is better to implement `Pandas_Workaround` or just use Pandas

    # optionally audit preprocessing
    if verbose:
        for entry in preprocess(tables):
            #click.echo(f"[{entry.level.upper():7}] {entry.step:20} {entry.sheet:15} {entry.message}")
            # indent every line…
            indent = "              "
            line = f"{entry.step:20} {entry.sheet:15} {entry.message}"
            # color by level
            click.echo("")  # blank line before mapping output
            if entry.level == "error":
                colored = click.style(line, fg="red")
            elif entry.level in ("warn", "warning"):
                colored = click.style(line, fg="yellow")
            else:
                colored = click.style(line, fg="cyan")
            click.echo(indent + colored)
        click.echo("")  # a blank line before mapping output

    # 4) Apply mapping to get raw records and collect issues
    notepad = create_notepad("phenopackets")
    genotype_records, phenotype_records = mapper.apply_mapping(tables, notepad)

    # 5) Report any errors or warnings
    _report_issues(notepad)

    # 6) Group results by patient
    records_by_patient = _group_records_by_patient(genotype_records, phenotype_records)

    # 7) Prepare output directory with timestamp
    # Will contain genotype and phenotype records as JSON
    generated_phenopacket_output_dir = _prepare_output_dir()

    # 8) Serialize phenopackets per patient
    _write_phenopackets(records_by_patient, generated_phenopacket_output_dir)

    # 9) Final summary
    click.echo(
        f"Wrote {len(records_by_patient)} phenopacket files to {generated_phenopacket_output_dir}"
    )
    click.echo(f"Created {len(genotype_records)} Genotype objects")
    click.echo(f"Created {len(phenotype_records)} Phenotype objects")


def _locate_hpo_file(hpo_path: typing.Optional[str]) -> pathlib.Path:
    # pick HPO JSON: either custom or default
    if hpo_path:
        hpo_file = pathlib.Path(hpo_path)
    else:
        hpo_file = pathlib.Path("tests/data") / "hp.json"
    if not hpo_file.is_file():
        click.echo(f"Error: HPO file not found at {hpo_file}", err=True)
        sys.exit(1)
    return hpo_file


def _load_ontology(hpo_file: str) -> hpotk.MinimalOntology:
    # load ontology from JSON
    return hpotk.load_minimal_ontology(hpo_file)


# Comment out this function and all uses to just get a Pandas_Workaround by only using the original `tables = load_sheets_as_tables(excel_file)`, which only depends on `from .loader import load_sheets_as_tables`
def _read_sheets(excel_file: str) -> dict[str, pd.DataFrame]:
    # read each worksheet into a DataFrame
    return load_sheets_as_tables(excel_file)


def _report_issues(notepad):
    # if there were errors, show them
    if notepad.has_errors(include_subsections=True):
        click.echo("Errors found in mapping:")
        for err in notepad.errors():
            click.echo(f"- {err}")
    # show any warnings but keep going
    if notepad.has_warnings(include_subsections=True):
        click.echo("Warnings found in mapping:")
        for w in notepad.warnings():
            click.echo(f"- {w}")


def _group_records_by_patient(
    genotype_records: list, phenotype_records: list
) -> dict[str, dict[str, list]]:
    # Group genotype & phenotype records by patient ID
    records = defaultdict(lambda: {"genotype_records": [], "phenotype_records": []})
    for genotype in genotype_records:
        records[genotype.genotype_patient_ID]["genotype_records"].append(genotype)
    for phenotype in phenotype_records:
        records[phenotype.phenotype_patient_ID]["phenotype_records"].append(phenotype)
    return records


def _prepare_output_dir() -> pathlib.Path:
    # use YYYY-MM-DD_HH-MM-SS for human-readable timestamps
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    phenopacket_from_excel_output_dir = (
        pathlib.Path.cwd() / "phenopacket_from_excel" / timestamp / "phenopackets"
    )
    phenopacket_from_excel_output_dir.mkdir(parents=True, exist_ok=True)
    return phenopacket_from_excel_output_dir


def _write_phenopackets(
    records_by_patient: dict[str, dict[str, list]],
    generated_phenopacket_output_dir: pathlib.Path,
):
    # Build and write one Phenopacket per patient
    for patient_id, patient_data in records_by_patient.items():
        phenopacket = Phenopacket()
        phenopacket.id = patient_id
        phenopacket.subject.id = patient_id

        # 3a) Add phenotypic features
        for phenotype in patient_data["phenotype_records"]:
            feature = phenopacket.phenotypic_features.add()
            feature.type.id = phenotype.HPO_ID
            # mark as excluded if status is False
            if not phenotype.status:
                feature.excluded = True

        # 3b) Add genotype interpretations
        # Genotypes → Interpretation → Diagnosis → GenomicInterpretation
        for interpretation_index, genotype_record in enumerate(
            patient_data["genotype_records"]
        ):
            # Create a new Interpretation entry (must set id and progress_status)
            interpretation = phenopacket.interpretations.add()
            interpretation.id = f"{patient_id}-interpretation-{interpretation_index}"
            interpretation.progress_status = interpretation.ProgressStatus.COMPLETED

            # each Interpretation has a `diagnosis` submessage (no `id` field)
            diagnosis = interpretation.diagnosis

            # inside Diagnosis, add GenomicInterpretation
            genomic_interpretation_entry = diagnosis.genomic_interpretations.add()
            genomic_interpretation_entry.subject_or_biosample_id = patient_id
            genomic_interpretation_entry.interpretation_status = (
                genomic_interpretation_entry.InterpretationStatus.CONTRIBUTORY
            )

            # TODO: Revise VariationDescriptor and gene_context later, omit setting gene_context for now.
            # variation_descriptor = genomic_interpretation_entry.variant_interpretation.variation_descriptor
            # we can also set variation_descriptor.gene_context and variation_descriptor.allelic_state here then serialize out as before
            # variation_descriptor.gene_context.gene_symbol = genotype_record.gene_symbol
            # variation_descriptor.allelic_state = variation_descriptor.AllelicState.Value(genotype_record.zygosity.upper())

            # Grab the VariantInterpretation and its descriptor
            variant_interpretation = genomic_interpretation_entry.variant_interpretation
            variation_descriptor = variant_interpretation.variation_descriptor

            # 1) Gene symbol & allelic state
            # 'gene_context' is a message; we need to CopyFrom if setting a message,
            # but for its scalar fields we can still assign directly:
            variation_descriptor.gene_context.symbol = genotype_record.gene_symbol
            variation_descriptor.allelic_state.CopyFrom(
                pps2.OntologyClass(
                    id="GENO:"
                    + genotype_record.zygosity_code,  # or however we decide to construct this later on
                    label=genotype_record.zygosity,
                )
            )

            # 2) HGVS expression
            hgvs_expr = variation_descriptor.expressions.add()
            # Attempt to set the HGVS syntax enum if available; otherwise skip.
            try:
                hgvs_expr.syntax = pps2.VariationDescriptor.Expression.HGVS
            except AttributeError:
                pass
            hgvs_expr.value = genotype_record.hgvsg

            # 3) Genomic location (exact interval) and alleles, if supported
            try:
                loc_ctx = variation_descriptor.location
                # use the nested VariationDescriptor.Location enum
                loc_ctx.interval.interval_type = (
                    pps2.VariationDescriptor.Location.Interval.Type.EXACT
                )
                loc_ctx.interval.start = genotype_record.start_position
                loc_ctx.interval.end = genotype_record.end_position
                loc_ctx.reference_sequence_id = genotype_record.chromosome

                # 4) Reference & alternate alleles
                variation_descriptor.reference = genotype_record.reference
                variation_descriptor.alternate = genotype_record.alternate
            except AttributeError:
                # some protobuffs give trouble when trying to expose location/alleles so just skip
                pass

        # 3d) Serialize to JSON
        generated_phenopacket_output_path = (
            generated_phenopacket_output_dir / f"{patient_id}.json"
        )
        with open(generated_phenopacket_output_path, "w", encoding="utf-8") as out_f:
            out_f.write(MessageToJson(phenopacket))

def preprocess(tables: dict[str, pd.DataFrame]) -> list[AuditEntry]:
    """
    Run lightweight audits on each sheet:
      - header normalization
      - sheet classification
      - variant‐column presence (raw vs HGVS)
    """
    from .mapper import (
        RAW_VARIANT_COLUMNS,
        HGVS_VARIANT_COLUMNS,
        GENOTYPE_BASE_COLUMNS,
        PHENOTYPE_KEY_COLUMNS,
    )

    entries: list[AuditEntry] = []

    # Step 1: header counts
    for name, df in tables.items():
        entries.append(AuditEntry(
            step="normalize-headers",
            sheet=name,
            message=f"{len(df.columns)} cols",
            level="info",
        ))

    # Step 2: classify
    for name, df in tables.items():
        cols = set(df.columns)
        has_raw = RAW_VARIANT_COLUMNS.issubset(cols)
        has_hgvs = bool(HGVS_VARIANT_COLUMNS & cols)
        is_gen = GENOTYPE_BASE_COLUMNS.issubset(cols) and (has_raw or has_hgvs)
        is_pheno = PHENOTYPE_KEY_COLUMNS.issubset(cols)

        kind = "genotype" if is_gen else "phenotype" if is_pheno else "skip"
        entries.append(AuditEntry(
            step="classify-sheet",
            sheet=name,
            message=kind + (f" ({'raw+hgvs' if has_raw and has_hgvs else 'raw' if has_raw else 'hgvs'})"),
            level="info",
        ))

    # Step 3: variant columns
    for name, df in tables.items():
        cols = set(df.columns)
        if GENOTYPE_BASE_COLUMNS.issubset(cols):
            if not (RAW_VARIANT_COLUMNS.issubset(cols) or HGVS_VARIANT_COLUMNS & cols):
                entries.append(AuditEntry(
                    step="variant-check",
                    sheet=name,
                    message="missing raw & HGVS",
                    level="error",
                ))
    return entries


if __name__ == "__main__":
    main()
