"""
Genotype domain model.

Defines the Genotype class which encapsulates all relevant fields for a
genomic variant entry, with validation of each attribute.
"""

import re
from dataclasses import dataclass
# from typing import ClassVar, Dict, List, Optional

# Patterns and allowed enums
_VALID_ID = re.compile(r"^[A-Za-z0-9]+$")
_EMAIL_PATTERN = re.compile(r"^[\w\.\+\-]+@[\w\.\-]+\.[A-Za-z]+$")
_ALLOWED_CHROM_ENCODINGS = {"hgvs", "ucsc", "refseq", "ensembl", "ncbi", "ega"}
_ALLOWED_ZYGOSITIES = {
    "compound_heterozygosity",
    "homozygous",
    "heterozygous",
    "hemizygous",
    "mosaic",
}
_ALLOWED_INHERITANCE_MODES = {"unknown", "inherited", "de_novo_mutation"}
# Mapping from the normalized zygosity terms to Genotype Ontology codes
_GENO_ALLELIC_STATE_CODES = {
    "heterozygous": "0000135",
    "homozygous": "0000134",
    "compound_heterozygosity": "0000191",
    "hemizygous": "0000136",
    "mosaic": "0000150",
}


@dataclass
class Genotype:
    """
    Represents a single genomic variant call for a patient.

    Attributes:
        genotype_patient_ID: Unique alphanumeric patient identifier.
        contact_email: Email for follow‑up communications.
        phasing: True if variant is phased, False otherwise.
        chromosome: Chromosome name or encoding (e.g., 'chr16' or 'hgvs').
        start_position: 1‑based start coordinate (nonnegative integer).
        end_position: 1‑based end coordinate (nonnegative integer).
        reference: Reference allele sequence.
        alternate: Alternate allele sequence.
        gene_symbol: Official gene symbol (e.g., “BRCA1”).
        hgvsg: HGVS genomic notation (e.g., “g.100A>T”).
        hgvsc: HGVS coding DNA notation (e.g., “c.200A>T”).
        hgvsp: HGVS protein notation (e.g., “p.Lys67Asn”).
        zygosity: One of the allowed zygosity terms.
        inheritance: One of the allowed inheritance modes.
    """

    genotype_patient_ID: str
    contact_email: str
    phasing: bool
    chromosome: str
    start_position: int
    end_position: int
    reference: str
    alternate: str
    gene_symbol: str
    hgvsg: str
    hgvsc: str
    hgvsp: str
    zygosity: str
    inheritance: str

    def __post_init__(self):
        # Validate patient ID
        if not _VALID_ID.match(self.genotype_patient_ID):
            raise ValueError(f"Invalid patient ID: {self.genotype_patient_ID!r}")

        # Validate email format
        if not _EMAIL_PATTERN.match(self.contact_email):
            raise ValueError(f"Invalid contact email: {self.contact_email!r}")

        # Validate chromosome: allow either a known encoding or real 'chr*' names
        chrom_lower = self.chromosome.lower()
        if not (
            chrom_lower in _ALLOWED_CHROM_ENCODINGS or chrom_lower.startswith("chr")
        ):
            raise ValueError(f"Unrecognized chromosome: {self.chromosome!r}")

        # Validate positions
        for attr in ("start_position", "end_position"):
            val = getattr(self, attr)
            if not isinstance(val, int) or val < 0:
                raise ValueError(f"{attr} must be a non‑negative integer, got {val!r}")

        # Validate allele/gene/HGVS strings
        for attr in (
            "reference",
            "alternate",
            "gene_symbol",
            "hgvsg",
            "hgvsc",
            "hgvsp",
        ):
            val = getattr(self, attr)
            if not isinstance(val, str) or not val.strip():
                raise ValueError(f"{attr} must be a nonempty string")

        # Validate zygosity
        if self.zygosity not in _ALLOWED_ZYGOSITIES:
            raise ValueError(f"Invalid zygosity: {self.zygosity!r}")

        # Validate inheritance
        if self.inheritance not in _ALLOWED_INHERITANCE_MODES:
            raise ValueError(f"Invalid inheritance mode: {self.inheritance!r}")

    @property
    def zygosity_code(self) -> str:
        """
        Returns the numeric portion of the GENO: alleleic_state code
        corresponding to this Genotype's zygosity.
        """
        try:
            return _GENO_ALLELIC_STATE_CODES[self.zygosity]
        except KeyError:
            raise ValueError(f"No GENO code defined for zygosity {self.zygosity!r}")


# Map our human‐readable zygosity → the GA4GH GENO codes
# allelic_state_GENO_zygosity_codes: dict[str, str] = {"heterozygous": "0000135", "homozygous": "0000136", "mosaic": "0000539", "hemizygous": "0000144", "compound_heterozygosity": "0000140"}
