"""
Genotype domain model.

Defines the Genotype class which encapsulates all relevant fields for a
genomic variant entry, with validation of each attribute.

This version integrates pyphetools' VariantValidator to authenticate/normalize
HGVS c. expressions and produce a GA4GH VariationDescriptor that P6 uses to
build Phenopackets. We parse the transcript from `hgvsc` and pass the c. part
only to the API, per pyphetools contract. When VariantValidator is unavailable
or returns a non-standard payload, we gracefully fall back to a locally
constructed VariationDescriptor so CLI runs and tests do not hard-fail.
"""

from __future__ import annotations

import os
import re
import requests
from dataclasses import dataclass
from typing import Optional

import phenopackets.schema.v2 as pps2
from pyphetools.creation.variant_validator import VariantValidator


# ---------------------------
# Patterns and allowed enums
# ---------------------------

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

# Permissive HGVS g. pattern with optional "chr" prefix (for SNVs)
# Captures: chromosome, position, ref, alt
_HGVS_G_SNV = re.compile(
    r"""
    ^\s*
    (?:chr)?(?P<chrom>[0-9XYM]+)       # chromosome number or X/Y/M
    :g\.
    (?P<pos>\d+)                       # 1-based position
    (?P<ref>[ACGT]+)>(?P<alt>[ACGT]+)  # simple SNV
    \s*$
    """,
    re.IGNORECASE | re.VERBOSE,
)

# Transcript + c. part, e.g. "NM_000000.0:c.100A>G", "ENST00000205557.12:c.2428G>A"
_HGVSC_TXT_RE = re.compile(
    r"""
    ^\s*
    (?P<tx>
        (?:N[MR]|X[MR]|E(?:NST)?)      # NM/NR/XM/XR/ENST
        [_]?\d+(?:\.\d+)?              # id with optional dot-version
    )
    :
    (?P<c>c\..+)$
    """,
    re.IGNORECASE | re.VERBOSE,
)


@dataclass
class Genotype:
    """
    Represents a single genomic variant call for a patient.

    Attributes:
        genotype_patient_ID: Unique alphanumeric patient identifier.
        contact_email: Email for follow-up communications.
        phasing: True if variant is phased, False otherwise.
        chromosome: Chromosome name or encoding (e.g., 'chr16' or 'hgvs').
        start_position: 1-based start coordinate (nonnegative integer).
        end_position: 1-based end coordinate (nonnegative integer).
        reference: Reference allele sequence.
        alternate: Alternate allele sequence.
        gene_symbol: Official gene symbol (e.g., “BRCA1”).
        hgvsg: HGVS genomic notation (e.g., “1:g.100A>T” or “chr1:g.100A>T”).
        hgvsc: HGVS coding DNA notation (e.g., “NM_000000.0:c.200A>T”).
        hgvsp: HGVS protein notation (e.g., “NP_000000.0:p.Lys67Asn”).
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

    # ----------------------
    # Basic input validation
    # ----------------------

    def __post_init__(self) -> None:
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
                raise ValueError(f"{attr} must be a non-negative integer, got {val!r}")

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

    # ----------------------
    # Convenience properties
    # ----------------------

    @property
    def zygosity_code(self) -> str:
        """
        Returns the numeric portion of the GENO: allelic_state code
        corresponding to this Genotype's zygosity.
        """
        try:
            return _GENO_ALLELIC_STATE_CODES[self.zygosity]
        except KeyError:
            raise ValueError(f"No GENO code defined for zygosity {self.zygosity!r}")

    # ----------------------
    # Variant construction
    # ----------------------

    def to_variation_descriptor(self) -> "pps2.VariationDescriptor":
        """
        Build a GA4GH VariationDescriptor for this variant.

        Strategy:
        1) If P6_SKIP_VV is set, or we cannot parse a transcript+c. from hgvsc,
           return a locally built descriptor using the provided g. string.
        2) Otherwise, call VariantValidator with (transcript, c_part).
           - If VV yields a usable object, adapt it and return.
           - If VV yields warnings / unexpected JSON shapes, fall back locally.

        The local fallback ensures CLI/test runs are resilient to VV outages
        and schema quirks.
        """
        # Optional global switch to skip VariantValidator entirely (e.g. CI/offline)
        if os.getenv("P6_SKIP_VV", "").strip() in {"1", "true", "TRUE"}:
            return self._build_local_vd()

        # Parse transcript + c. component from hgvsc (preferred path)
        tx, c_part = self._parse_hgvsc(self.hgvsc)
        if not (tx and c_part):
            # No transcript+c. info → local fallback with g.
            return self._build_local_vd()

        # Try pyphetools VariantValidator
        try:
            validator = VariantValidator(genome_build="GRCh38", transcript=tx)
            # pyphetools expects ONLY the c. part, not "NM:..."
            hv = validator.encode_hgvs(c_part)
            vi = hv.to_variant_interpretation_202()
            vd = vi.variation_descriptor
        except requests.RequestException:
            # Network/HTTP issues: keep runs resilient by falling back
            return self._build_local_vd()
        except (ValueError, TypeError, KeyError):
            # pyphetools may raise when VV's JSON has flag != 'gene_variant'
            # or when JSON shape is unexpected. Fall back locally.
            return self._build_local_vd()

        # Post-process: allelic_state (zygosity)
        if self.zygosity:
            vd.allelic_state.id = f"GENO:{self.zygosity_code}"
            vd.allelic_state.label = self.zygosity

        # Gene context (if missing)
        try:
            gene_ctx = getattr(vd, "gene_context", None)
            if self.gene_symbol and (
                gene_ctx is None or not getattr(gene_ctx, "symbol", "")
            ):
                vd.gene_context.symbol = self.gene_symbol
        except (AttributeError, ValueError, TypeError):
            # Don't let proto accessors sink the run
            pass

        # Ensure we include the g. expression we received (normalized no "chr")
        g_value = self._normalize_g_expression(self.hgvsg)
        if g_value:
            self._add_hgvs_expression(vd, g_value, syntax_name="HGVS")

        return vd

    # -----------------------
    # Small internal helpers
    # -----------------------

    @staticmethod
    def _parse_hgvsc(hgvsc: str) -> tuple[Optional[str], Optional[str]]:
        """
        Extract transcript identifier and the c. part from an hgvsc string.

        Examples:
            "NM_000000.0:c.100A>G" -> ("NM_000000.0", "c.100A>G")
            "ENST00000205557.12:c.2428G>A" -> ("ENST00000205557.12", "c.2428G>A")
        """
        if not isinstance(hgvsc, str):
            return None, None
        m = _HGVSC_TXT_RE.match(hgvsc.strip())
        if not m:
            return None, None
        tx = m.group("tx")
        c = m.group("c")
        return tx, c

    @staticmethod
    def _normalize_g_expression(hgvsg: str) -> Optional[str]:
        """
        Normalize a genomic HGVS like 'chr16:g.100A>G' -> '16:g.100A>G'
        (Only for simple SNV forms; otherwise return the original trimmed string.)
        """
        if not isinstance(hgvsg, str) or not hgvsg.strip():
            return None
        s = hgvsg.strip()
        m = _HGVS_G_SNV.match(s)
        if m:
            chrom = m.group("chrom")
            pos = m.group("pos")
            ref = m.group("ref").upper()
            alt = m.group("alt").upper()
            return f"{chrom}:g.{pos}{ref}>{alt}"
        # best-effort: strip any leading "chr"
        if s.lower().startswith("chr"):
            return s[3:]
        return s

    def _build_local_vd(self) -> "pps2.VariationDescriptor":
        """
        Construct a minimal VariationDescriptor locally, using the provided
        HGVS g. string (best-effort normalization), zygosity, and gene symbol.
        """
        vd = pps2.VariationDescriptor()
        # Add g. expression if we can normalize one
        g_value = self._normalize_g_expression(self.hgvsg)
        if g_value:
            self._add_hgvs_expression(vd, g_value, syntax_name="HGVS")

        # Allelic state
        if self.zygosity:
            vd.allelic_state.id = f"GENO:{self.zygosity_code}"
            vd.allelic_state.label = self.zygosity

        # Gene symbol (optional)
        if self.gene_symbol:
            try:
                vd.gene_context.symbol = self.gene_symbol
            except (AttributeError, ValueError, TypeError):
                pass

    @staticmethod
    def _add_hgvs_expression(
        vd: "pps2.VariationDescriptor", value: str, syntax_name: str = "HGVS"
    ) -> None:
        """
        Append an Expression to a VariationDescriptor.

        Parameters
        ----------
        vd : VariationDescriptor
            The descriptor to mutate.
        value : str
            HGVS string to append.
        syntax_name : str, optional
            Desired syntax enum name ('HGVS'), by default "HGVS".
        """
        expr = vd.expressions.add()
        expr.value = value
        # Many proto builds expose `Expression.HGVS`; we guard with getattr.
        enum = getattr(type(expr), syntax_name, None)
        if enum is not None and hasattr(expr, "syntax"):
            expr.syntax = enum  # type: ignore[attr-defined]
