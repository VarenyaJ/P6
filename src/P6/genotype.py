"""
Genotype domain model.

Defines the Genotype class which encapsulates all relevant fields for a
genomic variant entry, with validation of each attribute.

High-level role in P6:
- Mapper parses Excel rows → builds Genotype objects.
- Genotype.to_variation_descriptor() returns a GA4GH VariationDescriptor.
- The CLI assembles VariationDescriptors + HPO features into Phenopackets.

VariantDescriptor construction strategy:
1) Prefer pyphetools' VariantValidator (VV) using the transcript+c. parsed
   from `hgvsc` (e.g., "ENST00000205557.12:c.2428G>A"). If VV is reachable
   and returns a usable object, adapt it directly.
2) If VV is disabled (P6_SKIP_VV=1/true), missing, offline, or returns a
   non-standard payload, fall back to a minimal local descriptor:
   - include a normalized g.HGVS (strip "chr"), gene symbol, and zygosity.
3) Always de-duplicate expressions so we don't add the same g.HGVS twice.


Environment flags
----------------------------------------
P6_SKIP_VV=1           : Force the local fallback path (useful for CI/offline).
P6_ENRICH_GENE_XREFS=1 : If set and vv_lookup is importable, ask VV for HGNC/
                         Ensembl IDs and add them to gene_context where possible.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Optional, Tuple

import requests
import phenopackets.schema.v2 as pps2
from pyphetools.creation.variant_validator import VariantValidator


# ----------------------------------
# Patterns and small constant tables
# ----------------------------------

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

# GENO allelic_state codes mapped from normalized zygosity terms
_GENO_ALLELIC_STATE_CODES = {
    "heterozygous": "0000135",
    "homozygous": "0000134",
    "compound_heterozygosity": "0000191",
    "hemizygous": "0000136",
    "mosaic": "0000150",
}

# Permissive HGVS g. SNV pattern with optional "chr" prefix (captures chrom/pos/ref/alt)
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


# ----------------------
# Core domain data class
# ----------------------


@dataclass
class Genotype:
    """
    Represents a single genomic variant call for a patient.

    Attributes:
        genotype_patient_ID: Unique alphanumeric patient identifier.
        contact_email: Email for follow-up communications.
        phasing: True if variant is phased, False otherwise.
        chromosome: Chromosome name/encoding (e.g., 'chr16' or 'hgvs').
        start_position: 1-based start coordinate (non-negative integer).
        end_position: 1-based end coordinate (non-negative integer).
        reference: Reference allele sequence.
        alternate: Alternate allele sequence.
        gene_symbol: Official HGNC gene symbol (e.g., “BRCA1”).
        hgvsg: Genomic HGVS notation (e.g., “1:g.100A>T” or “chr1:g.100A>T”).
        hgvsc: Coding DNA HGVS notation (e.g., “NM_000000.0:c.200A>T”).
        hgvsp: Protein HGVS notation (e.g., “NP_000000.0:p.Lys67Asn”).
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

    # -----------------------------
    # Input validation on init time
    # -----------------------------

    def __post_init__(self) -> None:
        """Validate basic identifier formats and required string fields."""
        if not _VALID_ID.match(self.genotype_patient_ID):
            raise ValueError(f"Invalid patient ID: {self.genotype_patient_ID!r}")

        if not _EMAIL_PATTERN.match(self.contact_email):
            raise ValueError(f"Invalid contact email: {self.contact_email!r}")

        chrom_lower = self.chromosome.lower()
        if not (
            chrom_lower in _ALLOWED_CHROM_ENCODINGS or chrom_lower.startswith("chr")
        ):
            raise ValueError(f"Unrecognized chromosome: {self.chromosome!r}")

        for attr in ("start_position", "end_position"):
            val = getattr(self, attr)
            if not isinstance(val, int) or val < 0:
                raise ValueError(f"{attr} must be a non-negative integer, got {val!r}")

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

        if self.zygosity not in _ALLOWED_ZYGOSITIES:
            raise ValueError(f"Invalid zygosity: {self.zygosity!r}")

        if self.inheritance not in _ALLOWED_INHERITANCE_MODES:
            raise ValueError(f"Invalid inheritance mode: {self.inheritance!r}")

    # ----------------------
    # Convenience properties
    # ----------------------

    @property
    def zygosity_code(self) -> str:
        """Return the numeric part of the GENO: allelic_state code for this zygosity."""
        try:
            return _GENO_ALLELIC_STATE_CODES[self.zygosity]
        except KeyError as e:
            raise ValueError(
                f"No GENO code defined for zygosity {self.zygosity!r}"
            ) from e

    # --------------------------------------------------------------------------
    # Core responsibility: build a VariationDescriptor (VV path or local fallback)
    # --------------------------------------------------------------------------

    def to_variation_descriptor(self) -> "pps2.VariationDescriptor":
        """
        Build a GA4GH VariationDescriptor for this variant.

        Resolution order:
        1) If P6_SKIP_VV is set → build locally.
        2) Else, parse transcript+c. from hgvsc and attempt VV via pyphetools.
           - If VV returns a usable object: enrich & return.
           - On VV/network/shape issues: build locally.
        All paths de-duplicate expressions to avoid double g.HGVS entries.
        """
        if os.getenv("P6_SKIP_VV", "").strip().lower() in {"1", "true"}:
            vd = self._build_local_descriptor()
            return self._enrich_descriptor_common(vd)

        tx, c_part = self._parse_hgvsc(self.hgvsc)
        if not (tx and c_part):
            vd = self._build_local_descriptor()
            return self._enrich_descriptor_common(vd)

        # Try building via VariantValidator; keep failures graceful.
        try:
            vv = VariantValidator(genome_build="GRCh38", transcript=tx)
            hv = vv.encode_hgvs(c_part)  # pyphetools expects ONLY the c. part
            vi = hv.to_variant_interpretation_202()
            vd = vi.variation_descriptor
        except (
            requests.RequestException,
            ValueError,
            TypeError,
            AttributeError,
            KeyError,
        ):
            vd = self._build_local_descriptor()

        return self._enrich_descriptor_common(vd)

    # -----------------------
    # Internal helper methods
    # -----------------------

    @staticmethod
    def _parse_hgvsc(hgvsc: str) -> Tuple[Optional[str], Optional[str]]:
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
        return m.group("tx"), m.group("c")

    @staticmethod
    def _normalize_g_expression(hgvsg: str) -> Optional[str]:
        """
        Normalize a genomic HGVS like 'chr16:g.100A>G' -> '16:g.100A>G' for simple SNVs.
        For non-SNV or non-matching patterns, return the trimmed original string.
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
        if s.lower().startswith("chr"):
            return s[3:]
        return s

    # ---- Descriptor builders --------------------------------------------------

    def _build_local_descriptor(self) -> "pps2.VariationDescriptor":
        """
        Construct a minimal local VariationDescriptor using normalized g.HGVS,
        gene symbol, and zygosity. Used when VV is unavailable or disabled.
        """
        vd = pps2.VariationDescriptor()

        # Add g. expression if present (local path never has any expressions yet)
        g_value = self._normalize_g_expression(self.hgvsg)
        if g_value:
            self._add_hgvs_expression(vd, g_value, syntax_name="HGVS")

        # Allelic state (GENO)
        if self.zygosity:
            vd.allelic_state.id = f"GENO:{self.zygosity_code}"
            vd.allelic_state.label = self.zygosity

        # Gene context (optional)
        if self.gene_symbol:
            try:
                vd.gene_context.symbol = self.gene_symbol
            except (AttributeError, TypeError):
                # Do not let proto accessors sink the run
                pass

        return vd

    def _enrich_descriptor_common(
        self, vd: "pps2.VariationDescriptor"
    ) -> "pps2.VariationDescriptor":
        """
        Common post-processing for both VV-built and locally-built descriptors:
        - ensure allelic_state matches our zygosity,
        - ensure gene_context has our gene symbol (if missing),
        - ensure a normalized g.HGVS expression is present exactly once.
        """
        # Allelic state
        if self.zygosity:
            vd.allelic_state.id = f"GENO:{self.zygosity_code}"
            vd.allelic_state.label = self.zygosity

        # Gene symbol (do not overwrite a non-empty symbol VV may have provided)
        try:
            gene_ctx = getattr(vd, "gene_context", None)
            if self.gene_symbol and (
                gene_ctx is None or not getattr(gene_ctx, "symbol", "")
            ):
                vd.gene_context.symbol = self.gene_symbol
        except (AttributeError, TypeError):
            pass

        # Add normalized g.HGVS only if not already present (dedupe patch)
        g_value = self._normalize_g_expression(self.hgvsg)
        if g_value:
            self._add_hgvs_expression_if_missing(vd, g_value, syntax_name="HGVS")

        return vd

    # ---- Expression utilities -------------------------------------------------

    @staticmethod
    def _expression_values(vd: "pps2.VariationDescriptor") -> set[str]:
        """Collect all current Expression.value strings for fast membership checks."""
        try:
            return {e.value for e in vd.expressions}
        except (AttributeError, TypeError):
            return set()

    @classmethod
    def _add_hgvs_expression_if_missing(
        cls, vd: "pps2.VariationDescriptor", value: str, *, syntax_name: str = "HGVS"
    ) -> None:
        """
        Add a new HGVS Expression only if an identical value is not already present.
        Keeps VV-returned expressions from being duplicated by local enrichment.
        """
        if not value:
            return
        if value in cls._expression_values(vd):
            return
        cls._add_hgvs_expression(vd, value, syntax_name=syntax_name)

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
        enum = getattr(type(expr), syntax_name, None)
        if enum is not None and hasattr(expr, "syntax"):
            expr.syntax = enum  # type: ignore[attr-defined]
