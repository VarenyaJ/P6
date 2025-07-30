import abc
import click
import hpotk
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
    ) -> typing.Sequence[Phenopacket]:
        genotype_records = []
        phenotype_records = []

        for sheet_name, df in tables.items():
            columns = set(df.columns)
            is_genotype_sheet = GENOTYPE_KEY_COLUMNS.issubset(columns)
            is_phenotype_sheet = PHENOTYPE_KEY_COLUMNS.issubset(columns)

            if is_genotype_sheet == is_phenotype_sheet:
                # ambiguous sheet should give a warning instead of an error
                notepad.add_warning(
                    f"⚠ Skipping {sheet_name!r}: cannot unambiguously classify as genotype or phenotype"
                )
                continue

            # rename the former-index column
            id_column = (
                "genotype_patient_ID" if is_genotype_sheet else "phenotype_patient_ID"
            )
            working = df.reset_index()
            original_index_col = working.columns[0]
            working = working.rename(columns={original_index_col: id_column})

            # Verify ordering of any renamed columns
            for field in EXPECTED_COLUMN_NEIGHBORS:
                if field in working.columns:
                    pass
                    # try:
                    #     _verify_column_order(working, field)
                    # except ValueError as e:
                    #     click.echo(f"Sheet {sheet_name!r}: {e}", err=True)
                    #     sys.exit(1)

            if is_genotype_sheet:
                for idx, row in working.iterrows():
                    # handle slash‑separated zygosity and inheritance
                    zyg_list = [
                        z.strip().lower() for z in str(row["zygosity"]).split("/")
                    ]
                    inh_list = [
                        i.strip().lower() for i in str(row["inheritance"]).split("/")
                    ]

                    for z_code, i_code in zip(zyg_list, inh_list):
                        if z_code not in ZYGOSITY_MAP:
                            notepad.add_error(
                                f"Sheet {sheet_name!r}: Unrecognized zygosity code {z_code!r}"
                            )

                        if i_code not in INHERITANCE_MAP:
                            notepad.add_error(
                                f"Sheet {sheet_name!r}: Unrecognized inheritance code {i_code!r}"
                            )

                        ####
                        # allow missing/NaN contact_email → substitute dummy
                        raw_email = row["contact_email"]
                        if pd.isna(raw_email):
                            contact_email = "unknown@example.com"
                        else:
                            contact_email = str(raw_email).strip()
                        # ensure all “string” fields are actually strings:
                        kwargs = {
                            "genotype_patient_ID": str(row[id_column]),
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
                            "zygosity": ZYGOSITY_MAP[z_code],
                            "inheritance": INHERITANCE_MAP[i_code],
                        }
                        try:
                            g = Genotype(**kwargs)
                        except (ValueError, TypeError) as e:
                            notepad.add_error(f"Sheet {sheet_name!r}, row {idx}: {e}")
                            continue
                        genotype_records.append(g)
                        ####
            else:
                for idx, row in working.iterrows():
                    # normalize phenotype fields into valid strings
                    raw_hpo = row["hpo_id"]
                    hpo_str = str(raw_hpo).strip()
                    # extract the last token (it should just be the HPO code), case‑insensitive
                    m = re.search(r"(?:hp:)?(\d+)", hpo_str, re.IGNORECASE)
                    if not m:
                        notepad.add_error(
                            f"Sheet {sheet_name!r}, row {idx}: Cannot parse HPO ID from {hpo_str!r}"
                        )
                        continue
                    digits = m.group(1)
                    hpo_id = f"HP:{digits.zfill(7)}"
                    raw_date = row["date_of_observation"]
                    # if it's numeric, cast to int; else treat as string
                    if isinstance(raw_date, (int, float)):
                        date_str = f"T{int(raw_date)}"
                    else:
                        date_str = str(raw_date).strip()
                        if not date_str.upper().startswith("T"):
                            date_str = f"T{date_str}"

                    phenotype_records.append(
                        Phenotype(
                            phenotype_patient_ID=str(row[id_column]),
                            HPO_ID=hpo_id,
                            date_of_observation=date_str,
                            status=bool(row["status"]),
                        )
                    )

        # click.echo(f"Created {len(genotype_records)} Genotype objects")
        # click.echo(f"Created {len(phenotype_records)} Phenotype objects")
        return genotype_records, phenotype_records
