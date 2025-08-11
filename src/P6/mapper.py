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
from typing import Sequence, Dict, List, TypeVar, Tuple

T = TypeVar("T")
RowParseResult = Tuple[List[T], List[hpotk.TermId]]
# gives us one consistent return shape: (parsed_items, aux_ids_for_batch_validation)


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
KNOWN_SHEET_ALIASES: dict[str, set[str]] = {"genotype": {"genotype", "variants", "variant", "geno"},
                                            "phenotype": {"phenotype", "hpo", "pheno"},
                                            "diseases": {"disease", "diseases"},
                                            "measurements": {"measurement", "measurements", "labs"},
                                            "biosamples": {"biosample", "biosamples", "samples"}}


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
        # Map each selected sheet to domain-specific records via the table-level wrappers
        # The wrappers handle index→patient id normalization and any sheet-level checks then delegate to the row mappers.
        typed_tables = self._choose_named_tables(tables, notepad)
        genotype_records = self._map_genotype_table(typed_tables.genotype, notepad)
        phenotype_records = self._map_phenotype_table(typed_tables.phenotype, notepad)
        disease_records = self._map_diseases_table(typed_tables.diseases, notepad)
        measurement_records = self._map_measurements_table(typed_tables.measurements, notepad)
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

    def _prepare_sheet(self, df: pd.DataFrame, is_genotype: bool) -> pd.DataFrame:
        """Bring the index into a column and name it appropriately."""
        column_id = "genotype_patient_ID" if is_genotype else "phenotype_patient_ID"
        working = df.reset_index()
        original = working.columns[0]
        return working.rename(columns={original: column_id})

    @staticmethod
    def _normalize_time_like(value: typing.Any) -> str:
        """
        Phenotype/measurement/biosample timestamps:
        - numeric values are prefixed with 'T' (e.g., 20200101 -> 'T20200101')
        - strings are trimmed; if not already prefixed with 'T', we add it
        - empty/NaN -> empty string
        """
        # Handle None, NaN, NaT, pandas NA, and empty/whitespace-only strings
        if value is None or pd.isna(value) or (isinstance(value, str) and not value.strip()):
            return ""
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            return f"T{int(value)}"
        s = str(value).strip()
        if not s:
            return ""
        return s if s.upper().startswith("T") else f"T{s}"

    @staticmethod
    def _to_bool(value: typing.Any) -> bool:
        """
        Robust boolean parsing:
        - True for: 1, '1', 'true', 't', 'yes', 'y' (case-insensitive)
        - False for: 0, '0', 'false', 'f', 'no', 'n', '', None
        - Fallback: Python truthiness on other values (rare)
        """
        if isinstance(value, bool):
            return value
        if value is None:
            return False
        s = str(value).strip().lower()
        if s in {"1", "true", "t", "yes", "y"}:
            return True
        if s in {"0", "false", "f", "no", "n", ""}:
            return False
        return bool(value)

    @staticmethod
    def parse_genotype_row(row: pd.Series, sheet_name: str, notepad: Notepad) -> RowParseResult[Genotype]:
        """
        Parse a single genotype row into zero or more Genotype dataclass instances.
        Returns ([], []) if validation fails for this row.
        """
        genotypes: list[Genotype] = []

        # handle slash-separated zygosity and inheritance
        list_of_zygosity_types = [zygosity_entry.strip().lower() for zygosity_entry in
                                  str(row.get("zygosity", "")).split("/")]
        list_of_inheritance_types = [inheritance_entry.strip().lower() for inheritance_entry in
                                     str(row.get("inheritance", "")).split("/")]

        # zip will truncate to the shorter of the two, matching the previous behavior
        for zygosity_type, inheritance_type in zip(list_of_zygosity_types, list_of_inheritance_types):
            if zygosity_type not in ZYGOSITY_MAP:
                notepad.add_error(f"Sheet {sheet_name!r}: Unrecognized zygosity code {zygosity_type!r}")
                return [], []  # bail on this row
            if inheritance_type not in INHERITANCE_MAP:
                notepad.add_error(f"Sheet {sheet_name!r}: Unrecognized inheritance code {inheritance_type!r}")
                return [], []  # bail on this row

            # allow missing/NaN contact_email → substitute dummy
            raw_email = row.get("contact_email")
            contact_email = ("unknown@example.com" if pd.isna(raw_email) else str(raw_email).strip())

            try:
                genotypes.append(
                    Genotype(
                        genotype_patient_ID=str(row["genotype_patient_ID"]),
                        contact_email=contact_email,
                        phasing=DefaultMapper._to_bool(row.get("phasing")),
                        chromosome=str(row["chromosome"]),
                        start_position=int(row["start_position"]),
                        end_position=int(row["end_position"]),
                        reference=str(row["reference"]),
                        alternate=str(row["alternate"]),
                        gene_symbol=str(row["gene_symbol"]),
                        hgvsg=str(row["hgvsg"]),
                        hgvsc=str(row["hgvsc"]),
                        hgvsp=str(row["hgvsp"]),
                        zygosity=ZYGOSITY_MAP[zygosity_type],
                        inheritance=INHERITANCE_MAP[inheritance_type]
                    )
                )
            except (ValueError, TypeError) as e:
                notepad.add_error(f"Sheet {sheet_name!r}: {e}")
                return [], []  # treat any construction error as fatal for this row

        return genotypes, []  # no batch IDs for genotypes (yet)

    @staticmethod
    def parse_phenotype_row(row: pd.Series, hpo: hpotk.MinimalOntology, sheet_name: str, notepad: Notepad) -> RowParseResult[Phenotype]:
        """
        Parse a single phenotype row into zero or more Phenotype dataclasses.
        Also return any parsed TermIds so the caller can run batch validators later.
        Returns ([], []) if critical validation fails.
        """
        phenotypes: list[Phenotype] = []
        term_ids: list[hpotk.TermId] = []

        # normalize phenotype fields into valid strings
        hpo_cell = str(row.get("hpo_id", "")).strip()

        # Parse optional label and digits; extract the last token (it should just be the HPO code), case-insensitive
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
            notepad.add_error(f"Sheet {sheet_name!r}: Cannot parse HPO term+ID from {hpo_cell!r}")
            return [], []

        raw_label = m.group("label").strip()
        digits = m.group("digits")
        curie = f"HP:{digits.zfill(7)}"
        term_id = hpotk.TermId.from_curie(curie)

        # 1) Normalize the date_of_observation
        # if it's numeric, cast to int; else treat as string
        date_str = DefaultMapper._normalize_time_like(row.get("date_of_observation"))

        # 2) Append the Phenotype record
        try:
            phenotype = Phenotype(
                phenotype_patient_ID=str(row["phenotype_patient_ID"]),
                HPO_ID=curie,
                date_of_observation=date_str,
                status=DefaultMapper._to_bool(row.get("status")),
            )
        except (ValueError, TypeError) as e:
            notepad.add_error(f"Sheet {sheet_name!r}: {e}")
            return [], []

        phenotypes.append(phenotype)
        term_ids.append(term_id)

        # 3) The IDs must exist in the ontology:
        term = hpo.get_term(term_id)
        if term is None:
            notepad.add_warning(f"Sheet {sheet_name!r}: HPO ID {curie!r} not found in ontology")
        else:
            # 4) If the term is obsolete, flag it:
            if term.is_obsolete:
                replacements = ", ".join(str(t) for t in term.alt_term_ids)
                notepad.add_warning(f"Sheet {sheet_name!r}: {curie!r} is obsolete; use {replacements}")
            # 5) If they gave a label, check that it matches (case-insensitive):
            if raw_label and raw_label.lower() != term.name.lower():
                notepad.add_warning(
                    f"Sheet {sheet_name!r}: label {raw_label!r} does not match ontology name {term.name!r}")

        return phenotypes, term_ids

    def _map_genotype(self, sheet_name: str, df: pd.DataFrame, notepad: Notepad) -> list[Genotype]:
        records: list[Genotype] = []
        for _, row in df.iterrows():
            # Parse this row into zero or more Genotype records
            row_records, _ = self.parse_genotype_row(row, sheet_name, notepad)
            records.extend(row_records)
        return records

    def _map_phenotype(self, sheet_name: str, df: pd.DataFrame, notepad: Notepad) -> list[Phenotype]:
        records: list[Phenotype] = []
        # Collect every HPO ID in this sheet, so we can validate propagation later:
        all_ids: list[hpotk.TermId] = []

        for _, row in df.iterrows():
            row_records, row_term_ids = self.parse_phenotype_row(row, self._hpo, sheet_name, notepad)
            # Each parser returns lists; extend the accumulators.
            records.extend(row_records)
            all_ids.extend(row_term_ids)

        # Bulk-validate all collected IDs
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

    @staticmethod
    def check_hgvs_consistency(item: pd.Series, sheet_name: str, notepad: Notepad, strict: bool) -> None:
        """
        If both raw coordinates and HGVS notation are present, ensure that the genotype notations match
        """
        pattern = re.compile(
            r"^(?:chr)?(?P<chromosome_name>[^:]+):g\.(?P<mutation_position>\d+)"
            r"(?P<reference_allele>[ACGT]+)>(?P<alternative_allele>[ACGT]+)$",
            re.IGNORECASE,
        )
        hgvs = str(item.get("hgvsg", "")).strip()
        m = pattern.match(hgvs)
        if not m:
            notepad.add_error(f"Sheet {sheet_name!r}: malformed HGVS g. notation {hgvs!r}")
            return

        # mismatch = (str(item["chromosome"]) != m.group("chromosome_name") or int(item["start_position"]) != int(m.group("mutation_position")) or int(item["end_position"]) != int(m.group("mutation_position")) or str(item["reference"]) != m.group("reference_allele") or str(item["alternate"]) != m.group("alternative_allele")

        # Normalize cases and optional 'chr' prefix for robust comparison
        chrom_cell = str(item["chromosome"]).strip().lower()
        if chrom_cell.startswith("chr"):
            chrom_cell = chrom_cell[3:]
        chrom_hgvs = m.group("chromosome_name").strip().lower()

        pos_hgvs = int(m.group("mutation_position"))
        ref_cell = str(item["reference"]).strip().upper()
        alt_cell = str(item["alternate"]).strip().upper()
        ref_hgvs = m.group("reference_allele").upper()
        alt_hgvs = m.group("alternative_allele").upper()

        mismatch = (
            chrom_cell != chrom_hgvs
            or int(item["start_position"]) != pos_hgvs
            or int(item["end_position"]) != pos_hgvs
            or ref_cell != ref_hgvs
            or alt_cell != alt_hgvs
        )
        if mismatch:
            msg = (f"Sheet {sheet_name!r}: HGVS '{hgvs}' disagrees with "
                   f"raw ({item['chromosome']}:{item['start_position']}-"
                   f"{item['end_position']} {item['reference']}>{item['alternate']})")
            (notepad.add_error if strict else notepad.add_warning)(msg)

    def _prepare_sheet_for_patient(self, df: pd.DataFrame, patient_id_column: str) -> pd.DataFrame:
        """
        Similar to _prepare_sheet, but used for sheets whose patient identifier column is named 'patient_ID' (diseases, measurements, biosamples).
        This brings the current index into a named column so downstream mappers can access a consistent patient identifier.
        """
        working = df.reset_index()
        original = working.columns[0]
        return working.rename(columns={original: patient_id_column})

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

    # Table-level wrapper mappers
    def _map_genotype_table(self, df: pd.DataFrame | None, notepad: Notepad) -> list[Genotype]:
        """
        Sheet-level wrapper for Genotype rows:
          - normalize index to 'genotype_patient_ID'
          - require all key genotype columns (raw + HGVS)
          - optionally check HGVS vs raw coordinate consistency
          - delegate row conversion to _map_genotype
        """
        if df is None:
            return []
        working = self._prepare_sheet(df, is_genotype=True)

        have = set(working.columns)
        # Row parser and Genotype dataclass expect these columns to exist:
        missing = sorted(GENOTYPE_KEY_COLUMNS - have)
        if missing:
            notepad.add_error(f"Sheet 'genotype': missing required columns: {missing}")
            return []
        # Cross-check HGVS vs raw coordinates for every row (since both are present)
        for _, row in working.iterrows():
            self.check_hgvs_consistency(row, "genotype", notepad, self.strict_variants)

        # Must have the base columns, plus EITHER raw coordinates OR at least one HGVS field.
        #if not GENOTYPE_BASE_COLUMNS.issubset(have):
            #missing = sorted(GENOTYPE_BASE_COLUMNS - have)
            #notepad.add_error(f"Sheet 'genotype': missing required base columns: {missing}")
            #return []
        #if not (RAW_VARIANT_COLUMNS.issubset(have) or (HGVS_VARIANT_COLUMNS & have)):
            #notepad.add_error("Sheet 'genotype': provide either all raw coordinates ", f"{sorted(RAW_VARIANT_COLUMNS)} or at least one HGVS column ", f"from {sorted(HGVS_VARIANT_COLUMNS)}")
            #return []
        ## If BOTH groups are present, cross-check consistency:
        #columns_present = have
        #if RAW_VARIANT_COLUMNS.issubset(columns_present) and (HGVS_VARIANT_COLUMNS & columns_present):
            #for _, row in working.iterrows():
                #self.check_hgvs_consistency(row, "genotype", notepad, self.strict_variants)

        return self._map_genotype("genotype", working, notepad)

    def _map_phenotype_table(self, df: pd.DataFrame | None, notepad: Notepad) -> list[Phenotype]:
        """
        Sheet-level wrapper for Phenotype rows:
          - normalize index to 'phenotype_patient_ID'
          - delegate row conversion to _map_phenotype
        """
        if df is None:
            return []
        working = self._prepare_sheet(df, is_genotype=False)
        missing = PHENOTYPE_KEY_COLUMNS - set(working.columns)
        if missing:
            notepad.add_error(f"Sheet 'phenotype': missing expected columns: {sorted(missing)}")
            return []

        return self._map_phenotype("phenotype", working, notepad)

    def _map_diseases_table(self, df: pd.DataFrame | None, notepad: Notepad) -> list[DiseaseRecord]:
        """
        Sheet-level wrapper for Disease rows:
          - normalize index to 'patient_ID'
          - delegate row conversion to _map_disease
        """
        if df is None:
            return []
        working = self._prepare_sheet_for_patient(df, "patient_ID")
        return self._map_disease("diseases", working, notepad)

    def _map_measurements_table(self, df: pd.DataFrame | None, notepad: Notepad) -> list[MeasurementRecord]:
        """
        Sheet-level wrapper for Measurement rows:
          - normalize index to 'patient_ID'
          - delegate row conversion to _map_measurement
        """
        if df is None:
            return []
        working = self._prepare_sheet_for_patient(df, "patient_ID")
        return self._map_measurement("measurements", working, notepad)

    def _map_biosamples_table(self, df: pd.DataFrame | None, notepad: Notepad) -> list[BiosampleRecord]:
        """
        Sheet-level wrapper for Biosample rows:
          - normalize index to 'patient_ID'
          - delegate row conversion to _map_biosample
        """
        if df is None:
            return []
        working = self._prepare_sheet_for_patient(df, "patient_ID")
        return self._map_biosample("biosamples", working, notepad)

    def _map_disease(self, sheet_name: str, df: pd.DataFrame, notepad: Notepad) -> list[DiseaseRecord]:
        """
        Map each row in a disease sheet to a DiseaseRecord.
        Required columns: patient_ID, disease_term, disease_onset, disease_status.
        Optional column: disease_label.
        """
        records: list[DiseaseRecord] = []
        required_columns = {"patient_ID", "disease_term", "disease_onset", "disease_status"}
        missing = required_columns - set(df.columns)
        if missing:
            notepad.add_error(f"Sheet {sheet_name!r}: missing required columns: {sorted(missing)}")
            return records

        for index, row in df.iterrows():
            try:
                disease_record = DiseaseRecord(
                    patient_ID=str(row["patient_ID"]),
                    disease_term=str(row["disease_term"]).strip(),
                    disease_label=(str(row.get("disease_label", "")).strip() or None),
                    disease_onset=str(row["disease_onset"]).strip(),
                    disease_status=DefaultMapper._to_bool(row.get("disease_status")),
                )
                records.append(disease_record)
            except (ValueError, TypeError) as exception:
                notepad.add_error(f"Sheet {sheet_name!r}, row {index}: {exception}")
        return records

    def _map_measurement(self, sheet_name: str, df: pd.DataFrame, notepad: Notepad) -> list[MeasurementRecord]:
        """
        Map each row in a measurement sheet to a MeasurementRecord.
        Required columns: patient_ID, measurement_type, measurement_value, measurement_unit.
        Optional column: measurement_timestamp (numeric values are prefixed with 'T' for consistency).
        """
        records: list[MeasurementRecord] = []
        required_columns = {"patient_ID", "measurement_type", "measurement_value", "measurement_unit"}
        missing = required_columns - set(df.columns)
        if missing:
            notepad.add_error(f"Sheet {sheet_name!r}: missing required columns: {sorted(missing)}")
            return records

        for index, row in df.iterrows():
            try:
                measurement_timestamp = self._normalize_time_like(row.get("measurement_timestamp")) or None

                measurement_record = MeasurementRecord(
                    patient_ID=str(row["patient_ID"]),
                    measurement_type=str(row["measurement_type"]).strip(),
                    measurement_value=float(row["measurement_value"]),
                    measurement_unit=str(row["measurement_unit"]).strip(),
                    measurement_timestamp=measurement_timestamp,
                )
                records.append(measurement_record)
            except (ValueError, TypeError) as exception:
                notepad.add_error(f"Sheet {sheet_name!r}, row {index}: {exception}")
        return records

    def _map_biosample(self, sheet_name: str, df: pd.DataFrame, notepad: Notepad) -> list[BiosampleRecord]:
        """
        Map each row in a biosample sheet to a BiosampleRecord.
        Required columns: patient_ID, biosample_id, biosample_type, collection_date.
        Numeric collection_date values are prefixed with 'T' for consistency with phenotype dates.
        """
        records: list[BiosampleRecord] = []
        required_columns = {"patient_ID", "biosample_id", "biosample_type", "collection_date"}
        missing = required_columns - set(df.columns)
        if missing:
            notepad.add_error(f"Sheet {sheet_name!r}: missing required columns: {sorted(missing)}")
            return records

        for index, row in df.iterrows():
            try:
                collection_date = self._normalize_time_like(row.get("collection_date")) or ""

                biosample_record = BiosampleRecord(
                    patient_ID=str(row["patient_ID"]),
                    biosample_id=str(row["biosample_id"]).strip(),
                    biosample_type=str(row["biosample_type"]).strip(),
                    collection_date=collection_date,
                )
                records.append(biosample_record)
            except (ValueError, TypeError) as exception:
                notepad.add_error(f"Sheet {sheet_name!r}, row {index}: {exception}")
        return records

    # Grouping and phenopacket construction
    def _group_records_by_patient(self, genotype_records: list[Genotype], phenotype_records: list[Phenotype], disease_records: list[DiseaseRecord], measurement_records: list[MeasurementRecord], biosample_records: list[BiosampleRecord], ) -> dict[str, dict[str, list]]:
        """
        Group all domain records by patient identifier, producing a bundle per patient
        """
        grouped = defaultdict(
            lambda: {
                "genotype_records": [],
                "phenotype_records": [],
                "disease_records": [],
                "measurement_records": [],
                "biosample_records": [],
            }
        )
        for genotype in genotype_records:
            grouped[genotype.genotype_patient_ID]["genotype_records"].append(genotype)
        for phenotype in phenotype_records:
            grouped[phenotype.phenotype_patient_ID]["phenotype_records"].append(phenotype)
        for disease in disease_records:
            grouped[disease.patient_ID]["disease_records"].append(disease)
        for measurement in measurement_records:
            grouped[measurement.patient_ID]["measurement_records"].append(measurement)
        for biosample in biosample_records:
            grouped[biosample.patient_ID]["biosample_records"].append(biosample)
        return grouped

    def construct_phenopacket_for_patient(self, patient_id: str, bundle: dict[str, list],
                                          notepad: Notepad) -> Phenopacket:
        """
        Build a Phenopacket for a single patient using their grouped records.
        Field assignments follow the explicit naming and serialization style.
        """
        from phenopackets.schema.v2.phenopackets_pb2 import Phenopacket
        import phenopackets.schema.v2 as pps2

        phenopacket = Phenopacket()
        phenopacket.id = patient_id
        phenopacket.subject.id = patient_id

        # Phenotypic features
        for phenotype in bundle.get("phenotype_records", []):
            feature = phenopacket.phenotypic_features.add()
            feature.type.id = phenotype.HPO_ID
            if not phenotype.status:
                feature.excluded = True

        # Genotype interpretations (minimal HGVS expression to start)
        for interpretation_index, genotype_record in enumerate(bundle.get("genotype_records", [])):
            interpretation = phenopacket.interpretations.add()
            interpretation.id = f"{patient_id}-interpretation-{interpretation_index}"
            interpretation.progress_status = interpretation.ProgressStatus.COMPLETED

            genomic_interpretation_entry = interpretation.diagnosis.genomic_interpretations.add()
            genomic_interpretation_entry.subject_or_biosample_id = patient_id
            genomic_interpretation_entry.interpretation_status = (
                genomic_interpretation_entry.InterpretationStatus.CONTRIBUTORY
            )

            variant_interpretation = genomic_interpretation_entry.variant_interpretation
            variation_descriptor = variant_interpretation.variation_descriptor

            # Add HGVS expression; set syntax enum when available
            expression = variation_descriptor.expressions.add()
            try:
                expression.syntax = pps2.VariationDescriptor.Expression.HGVS
            except AttributeError:
                pass
            expression.value = genotype_record.hgvsg or ""

            # Optional: attempt to set a subset of location/alleles if supported
            try:
                location_context = variation_descriptor.location
                location_context.interval.interval_type = (
                    pps2.VariationDescriptor.Location.Interval.Type.EXACT
                )
                location_context.interval.start = genotype_record.start_position
                location_context.interval.end = genotype_record.end_position
                location_context.reference_sequence_id = genotype_record.chromosome
                variation_descriptor.reference = genotype_record.reference
                variation_descriptor.alternate = genotype_record.alternate
            except AttributeError:
                # Some library builds do not expose these submessages; try to skip gracefully if we cannot implement this feature.
                pass

        # Optional sections (diseases, measurements, biosamples).
        # Keep assignments minimal and consistent with the earlier CLI code.
        for disease_record in bundle.get("disease_records", []):
            disease_message = phenopacket.diseases.add()
            disease_message.term.id = disease_record.disease_term
            if getattr(disease_record, "disease_label", None):
                disease_message.term.label = disease_record.disease_label
            # Onset/status message wiring can be expanded later as needed.

        for measurement_record in bundle.get("measurement_records", []):
            measurement_message = phenopacket.measurements.add()
            measurement_message.type.id = measurement_record.measurement_type
            # Depending on the installed protobuf version, value/unit/timestamp may require message types;
            # keep this minimal for now (extend down the line).

        for biosample_record in bundle.get("biosample_records", []):
            biosample_message = phenopacket.biosamples.add()
            biosample_message.id = biosample_record.biosample_id
            biosample_message.type.id = biosample_record.biosample_type

        return phenopacket
