import abc

# ruff removed: click
import hpotk
from hpotk.validate import (
    ObsoleteTermIdsValidator,
    PhenotypicAbnormalityValidator,
    AnnotationPropagationValidator,
    ValidationRunner,
)
import pandas as pd
import re
import typing

from collections import defaultdict
from dataclasses import dataclass
from phenopackets.schema.v2.phenopackets_pb2 import Phenopacket
from stairval.notepad import Notepad

from .biosample import BiosampleRecord
from .disease import DiseaseRecord
from .genotype import Genotype
from .measurement import MeasurementRecord
from .phenotype import Phenotype


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

PHENOTYPE_KEY_COLUMNS = {"hpo_id", "date_of_observation", "status"}

# Key columns to identify additional sheets
DISEASE_KEY_COLUMNS = {"disease_term", "disease_onset"}
MEASUREMENT_KEY_COLUMNS = {"measurement_type", "measurement_value", "measurement_unit"}
BIOSAMPLE_KEY_COLUMNS = {"biosample_id", "biosample_type", "collection_date"}


# Map raw zygosity abbreviations to allowed dataclass zygosity values
ZYGOSITY_MAP = {
    "het": "heterozygous",
    "hom": "homozygous",
    "comphet": "compound_heterozygosity",
    "hemi": "hemizygous",
    "mosaic": "mosaic",
}

# Map raw inheritance abbreviations to allowed dataclass inheritance values
INHERITANCE_MAP = {
    "unknown": "unknown",
    "inherited": "inherited",
    "denovo": "de_novo_mutation",
}

# Check Variant column groups: allow either the full raw coordinates or any HGVS notation
RAW_VARIANT_COLUMNS = {
    "chromosome",
    "start_position",
    "end_position",
    "reference",
    "alternate",
}
HGVS_VARIANT_COLUMNS = {"hgvsg", "hgvsc", "hgvsp"}
# minimal base columns to call something a genotype sheet (we bring the index in later)
GENOTYPE_BASE_COLUMNS = {"contact_email", "phasing"}

# Friendly aliases → reduces friction while keeping behavior explicit
KNOWN_SHEET_ALIASES: dict[str, set[str]] = {"genotype": {"genotype", "variants", "variant", "geno"}, "phenotype": {"phenotype", "hpo", "pheno"}, "diseases": {"disease", "diseases"}, "measurements": {"measurement", "measurements", "labs"}, "biosamples": {"biosample", "biosamples", "samples"}}

@dataclass
class TypedTables:
    """
    Explicit, typed access to workbook sheets.
    Any field can be `None`, meaning that the sheet not provided.
    """
    genotype: pd.DataFrame | None
    phenotype: pd.DataFrame | None
    diseases: pd.DataFrame | None
    measurements: pd.DataFrame | None
    biosamples: pd.DataFrame | None

def _choose_named_tables(self, tables: dict[str, pd.DataFrame], notepad: Notepad) -> TypedTables:
    """
    Prefer explicit sheet names (plus common aliases).
    """
    def by_alias(kind: str) -> pd.DataFrame | None:
        aliases = KNOWN_SHEET_ALIASES[kind]
        for sheet_name, df in tables.items():
            if sheet_name.strip().casefold() in aliases:
                return df
        return None

    selected = TypedTables(
        genotype=by_alias("genotype"),
        phenotype=by_alias("phenotype"),
        diseases=by_alias("diseases"),
        measurements=by_alias("measurements"),
        biosamples=by_alias("biosamples"),
    )

    # Hard-minimum: at least genotype or phenotype must exist
    if selected.genotype is None and selected.phenotype is None:
        notepad.add_error("Missing required sheet: either 'genotype' or 'phenotype'.")

    return selected

class TableMapper(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def apply_mapping(
        self, tables: dict[str, pd.DataFrame], notepad: Notepad
    ) -> typing.Sequence[Phenopacket]:
        # return fully-assembled Phenopacket messages, not intermediate parts.
        raise NotImplementedError


class DefaultMapper(TableMapper):
    def __init__(self, hpo: hpotk.MinimalOntology, strict_variants: bool = False):
        """
        - False: raw⇄HGVS mismatches are logged as WARNINGS
        - True : raw⇄HGVS mismatches are logged as ERRORS
        """
        self._hpo = hpo
        self.strict_variants = strict_variants

    def apply_mapping(
        self, tables: dict[str, pd.DataFrame], notepad: Notepad
    ) -> list[Phenopacket]:
        """
        Process:
        1) choose/validate input tables
        2) map rows to domain records
        3) group records per patient
        4) construct Phenopacket per patient
        5) return the list of packets

        """
        # TODO: implement the placeholders I am going to temporarily call
        typed_tables = self._choose_named_tables(tables, notepad)
        genotype_records = self._map_genotype_table(typed_tables.genotype, notepad)
        phenotype_records = self._map_phenotype_table(typed_tables.phenotype, notepad)
        disease_records = self._map_diseases_table(typed_tables.diseases, notepad)
        measurement_records = self._map_measurements_table(
            typed_tables.measurements, notepad
        )
        biosample_records = self._map_biosamples_table(typed_tables.biosamples, notepad)

        grouped = self._group_records_by_patient(
            genotype_records,
            phenotype_records,
            disease_records,
            measurement_records,
            biosample_records,
        )

        packets: list[Phenopacket] = [
            self.construct_phenopacket_for_patient(patient_id, bundle, notepad)
            for patient_id, bundle in grouped.items()
        ]
        return packets

    def _check_hgvs_consistency(
        self, sheet_name: str, df: pd.DataFrame, notepad: Notepad
    ) -> None:
        """
        If both raw coordinates and HGVS notation are present, ensure that the genotype notations match
        """

        pattern = re.compile(
            r"^(?:chr)?(?P<chromosome_name>[^:]+):g\.(?P<mutation_position>\d+)"
            r"(?P<reference_allele>[ACGT]+)>(?P<alternative_allele>[ACGT]+)$",
            re.IGNORECASE,
        )

        for idx, row in df.iterrows():
            hgvs = str(row.get("hgvsg", "")).strip()
            m = pattern.match(hgvs)
            if not m:
                # always treat malformed HGVS as an error
                notepad.add_error(
                    f"Sheet {sheet_name!r}, row {idx}: malformed HGVS g. notation {hgvs!r}"
                )
                continue
            chromosome_name = m.group("chromosome_name")
            mutation_position = int(m.group("mutation_position"))
            reference_allele = m.group("reference_allele")
            alternative_allele = m.group("alternative_allele")

            # compare against raw columns
            mismatch_msg = (
                f"Sheet {sheet_name!r}, row {idx}: HGVS '{hgvs}' disagrees with "
                f"raw ({row['chromosome']}:{row['start_position']}-"
                f"{row['end_position']} {row['reference']}>{row['alternate']})"
            )
            if (
                str(row["chromosome"]) != chromosome_name
                or int(row["start_position"]) != mutation_position
                or int(row["end_position"]) != mutation_position
                or str(row["reference"]) != reference_allele
                or str(row["alternate"]) != alternative_allele
            ):
                if self.strict_variants:
                    notepad.add_error(mismatch_msg)
                else:
                    notepad.add_warning(mismatch_msg)

    def _prepare_sheet(self, df: pd.DataFrame, is_genotype: bool) -> pd.DataFrame:
        """Bring the index into a column and name it appropriately."""
        column_id = "genotype_patient_ID" if is_genotype else "phenotype_patient_ID"
        working = df.reset_index()
        original = working.columns[0]
        return working.rename(columns={original: column_id})

    def _map_genotype(
        self, sheet_name: str, df: pd.DataFrame, notepad: Notepad
    ) -> list[Genotype]:
        records: list[Genotype] = []
        for idx, row in df.iterrows():
            # handle slash‑separated zygosity and inheritance
            list_of_zygosity_types = [
                z.strip().lower() for z in str(row["zygosity"]).split("/")
            ]
            list_of_inheritance_types = [
                i.strip().lower() for i in str(row["inheritance"]).split("/")
            ]
            for zygosity_type, inheritance_type in zip(
                list_of_zygosity_types, list_of_inheritance_types
            ):
                if zygosity_type not in ZYGOSITY_MAP:
                    notepad.add_error(
                        f"Sheet {sheet_name!r}: Unrecognized zygosity code {zygosity_type!r}"
                    )
                if inheritance_type not in INHERITANCE_MAP:
                    notepad.add_error(
                        f"Sheet {sheet_name!r}: Unrecognized inheritance code {inheritance_type!r}"
                    )
                # allow missing/NaN contact_email → substitute dummy
                raw_email = row["contact_email"]
                contact_email = (
                    "unknown@example.com"
                    if pd.isna(raw_email)
                    else str(raw_email).strip()
                )
                kwargs = {
                    "genotype_patient_ID": str(row["genotype_patient_ID"]),
                    "contact_email": contact_email,
                    "phasing": bool(row["phasing"]),
                    "chromosome": str(row["chromosome"]),
                    "start_position": int(row["start_position"]),
                    "end_position": int(row["end_position"]),
                    "reference": str(row["reference"]),
                    "alternate": str(row["alternate"]),
                    "gene_symbol": str(row["gene_symbol"]),
                    "hgvsg": str(row["hgvsg"]),
                    "hgvsc": str(row["hgvsc"]),
                    "hgvsp": str(row["hgvsp"]),
                    "zygosity": ZYGOSITY_MAP[zygosity_type],
                    "inheritance": INHERITANCE_MAP[inheritance_type],
                }
                try:
                    records.append(Genotype(**kwargs))
                except (ValueError, TypeError) as e:
                    notepad.add_error(f"Sheet {sheet_name!r}, row {idx}: {e}")
        return records

    def _map_phenotype(
        self, sheet_name: str, df: pd.DataFrame, notepad: Notepad
    ) -> list[Phenotype]:
        records: list[Phenotype] = []
        # Collect every HPO ID in this sheet, so we can validate propagation later:
        all_ids: list[hpotk.TermId] = []

        for idx, row in df.iterrows():
            # normalize phenotype fields into valid strings
            hpo_cell = str(row["hpo_id"]).strip()
            # Parse optional label and digits
            # extract the last token (it should just be the HPO code), case‑insensitive
            m = re.match(
                r"""
                ^\s*
                (?P<label>.*?)              # optional label
                \s*                         # whitespace
                \(?                         # optional "("
                (?:HP:?)?(?P<digits>\d+)    # digits, with optional "HP"
                \)?                         # optional ")"
                \s*$
                """,
                hpo_cell,
                re.VERBOSE | re.IGNORECASE,
            )
            if not m:
                notepad.add_error(
                    f"Sheet {sheet_name!r}, row {idx}: Cannot parse HPO term+ID from {hpo_cell!r}"
                )
                continue

            raw_label = m.group("label").strip()
            digits = m.group("digits")
            curie = f"HP:{digits.zfill(7)}"
            term_id = hpotk.TermId.from_curie(curie)

            # 1) Normalize the date_of_observation
            # Normalize the date_of_observation
            # if it's numeric, cast to int; else treat as string
            raw_date = row["date_of_observation"]
            if isinstance(raw_date, (int, float)):
                date_str = f"T{int(raw_date)}"
            else:
                s = str(raw_date).strip()
                date_str = s if s.upper().startswith("T") else f"T{s}"

            # 2) Append the Phenotype record
            # Append Phenotype record
            records.append(
                Phenotype(
                    phenotype_patient_ID=str(row["phenotype_patient_ID"]),
                    HPO_ID=curie,
                    date_of_observation=date_str,
                    status=bool(row["status"]),
                )
            )

            # 3) The IDs must exist in the ontology:
            # Validate ID against ontology
            term = self._hpo.get_term(term_id)
            if term is None:
                notepad.add_warning(
                    f"Skipping row {idx} in {sheet_name!r}: HPO ID {curie!r} not found in ontology"
                )
                continue

            # 4) If the term is obsolete, flag it:
            if term.is_obsolete:
                replacements = ", ".join(str(t) for t in term.alt_term_ids)
                notepad.add_warning(
                    f"Sheet {sheet_name!r}, row {idx}: {curie!r} is obsolete; use {replacements}"
                )

            # 5) If they gave a label, check that it matches (case-insensitive):
            if raw_label and raw_label.lower() != term.name.lower():
                notepad.add_warning(
                    f"Sheet {sheet_name!r}, row {idx}: label {raw_label!r} "
                    f"does not match ontology name {term.name!r}"
                )

            # Only now record for batch‐validation
            all_ids.append(term_id)

        # Bulk‐validate all collected IDs
        if all_ids:
            validators = [
                ObsoleteTermIdsValidator(self._hpo),
                PhenotypicAbnormalityValidator(self._hpo),
                AnnotationPropagationValidator(self._hpo),
            ]
            runner = ValidationRunner(validators=validators)
            validation_runner = runner.validate_all(all_ids)
            for issue in validation_runner.results:
                msg = f"Sheet {sheet_name!r}: {issue.message}"
                if issue.level.name == "ERROR":
                    notepad.add_error(msg)
                else:
                    notepad.add_warning(msg)

        return records

    def _map_disease(
        self, sheet_name: str, df: pd.DataFrame, notepad: Notepad
    ) -> list[DiseaseRecord]:
        # TODO: implement row→DiseaseRecord, row→MeasurementRecord conversion, and row→BiosampleRecord conversions
        """
        Map each row in a disease sheet to a DiseaseRecord.
        """
        # TODO: fix as this is not in use now: records: list[DiseaseRecord] = []
        raise NotImplementedError

    def _map_measurement(
        self, sheet_name: str, df: pd.DataFrame, notepad: Notepad
    ) -> list[MeasurementRecord]:
        """
        Map each row in a measurement sheet to a MeasurementRecord.
        """
        # TODO: fix as this is not in use now: records: list[MeasurementRecord] = []
        raise NotImplementedError

    def _map_biosample(
        self, sheet_name: str, df: pd.DataFrame, notepad: Notepad
    ) -> list[BiosampleRecord]:
        """
        Map each row in a biosample sheet to a BiosampleRecord.
        """
        # TODO: fix as this is not in use now: records: list[BiosampleRecord] = []
        raise NotImplementedError
