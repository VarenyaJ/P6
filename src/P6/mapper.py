import abc

# ruff removed: click
import hpotk
from hpotk.validate import ObsoleteTermIdsValidator, PhenotypicAbnormalityValidator, AnnotationPropagationValidator, ValidationRunner
import pandas as pd
import re
import typing

from phenopackets.schema.v2.phenopackets_pb2 import Phenopacket
from stairval.notepad import Notepad

from .genotype import Genotype
from .phenotype import Phenotype


class TableMapper(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def apply_mapping(
        self, tables: dict[str, pd.DataFrame], notepad: Notepad
    ) -> typing.Sequence[Phenopacket]:
        pass


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


class DefaultMapper(TableMapper):
    def __init__(self, hpo: hpotk.MinimalOntology):
        self._hpo = hpo

    def apply_mapping(
        self, tables: dict[str, pd.DataFrame], notepad: Notepad
    ) -> tuple[list[Genotype], list[Phenotype]]:
        genotype_records: list[Genotype] = []
        phenotype_records: list[Phenotype] = []

        for sheet_name, df in tables.items():
            columns = set(df.columns)
            """Send each sheet to the right extractor and collect all records."""
            is_genotype_sheet = GENOTYPE_KEY_COLUMNS.issubset(columns)
            is_phenotype_sheet = PHENOTYPE_KEY_COLUMNS.issubset(columns)

            if is_genotype_sheet == is_phenotype_sheet:
                # ambiguous sheet should give a warning instead of an error
                notepad.add_warning(
                    f"Skipping {sheet_name!r}: cannot unambiguously classify as genotype or phenotype"
                )
                continue

            # rename the former-index column
            working = self._prepare_sheet(df, is_genotype_sheet)

            if is_genotype_sheet:
                genotype_records.extend(
                    self._map_genotype(sheet_name, working, notepad)
                )
            else:
                phenotype_records.extend(
                    self._map_phenotype(sheet_name, working, notepad)
                )

        return genotype_records, phenotype_records

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
                re.VERBOSE | re.IGNORECASE
            )
            if not m:
                notepad.add_error(f"Sheet {sheet_name!r}, row {idx}: Cannot parse HPO term+ID from {hpo_cell!r}")
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
                notepad.add_warning(f"Skipping row {idx} in {sheet_name!r}: HPO ID {curie!r} not found in ontology")
                continue

            # 4) If the term is obsolete, flag it:
            if term.is_obsolete:
                replacements = ", ".join(str(t) for t in term.alt_term_ids)
                notepad.add_warning(f"Sheet {sheet_name!r}, row {idx}: {curie!r} is obsolete; use {replacements}")

            # 5) If they gave a label, check that it matches (case-insensitive):
            if raw_label and raw_label.lower() != term.name.lower():
                notepad.add_warning(
                    f"Sheet {sheet_name!r}, row {idx}: label {raw_label!r} "
                    f"does not match ontology name {term.name!r}")

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