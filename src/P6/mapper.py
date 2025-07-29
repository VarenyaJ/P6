import abc
import typing

import hpotk
import pandas as pd

from phenopackets.schema.v2.phenopackets_pb2 import Phenopacket
from stairval.notepad import Notepad

from .genotype import Genotype
from .phenotype import Phenotype


class TableMapper(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def apply_mapping(
        self, 
        tables: dict[str, pd.DataFrame],
        notepad: Notepad,
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
    def __init__(self, hpo: hpotk.MinimalOntology,):
        self._hpo = hpo

    def apply_mapping(
        self,
        tables: dict[str, pd.DataFrame],
        notepad: Notepad,
    ) -> typing.Sequence[Phenopacket]:
        genotype_records = []
        phenotype_records = []

        for sheet_name, df in tables.items():
            columns = set(df.columns)
            is_genotype_sheet = GENOTYPE_KEY_COLUMNS.issubset(columns)
            is_phenotype_sheet = PHENOTYPE_KEY_COLUMNS.issubset(columns)

            if is_genotype_sheet == is_phenotype_sheet:
                notepad.add_error(
                    f"⚠  Skipping {sheet_name!r}: cannot unambiguously classify as genotype or phenotype",
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
                for _, row in working.iterrows():
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
                                f"Sheet {sheet_name!r}: Unrecognized zygosity code {z_code!r}",
                            )
                            
                        if i_code not in INHERITANCE_MAP:
                            notepad.add_error(
                                f"Sheet {sheet_name!r}: Unrecognized inheritance code {i_code!r}",
                            )

                        genotype_records.append(
                            Genotype(
                                genotype_patient_ID=str(row[id_column]),
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
                                zygosity=ZYGOSITY_MAP[z_code],
                                inheritance=INHERITANCE_MAP[i_code],
                            )
                        )
            else:
                for _, row in working.iterrows():
                    # --- normalize phenotype fields into valid strings ---
                    raw_hpo = row["hpo_id"]
                    hpo_str = str(raw_hpo).strip()
                    if hpo_str.lower().startswith("hp:"):
                        digits = hpo_str[3:]
                        hpo_id = f"HP:{digits.zfill(7)}"
                    else:
                        hpo_id = hpo_str.zfill(7)

                    raw_date = row["date_of_observation"]
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
        return ()
