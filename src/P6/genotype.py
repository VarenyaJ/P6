"""
Genotype domain model and VariationDescriptor construction.

High level (P6 perspective)
---------------------------
This module defines the Genotype record used by P6's mapper. Its main job is
to validate raw genotype fields coming from Excel/CSV, and to convert each row
into a GA4GH Phenopackets v2 `VariationDescriptor` that the CLI writes into
phenopackets.

Two-path strategy for variation building
----------------------------------------
1) Preferred: Use `pyphetools`' VariantValidator (VV) to authenticate a
   transcript+c. HGVS (parsed from `hgvsc`) and obtain a normalized descriptor.
2) Fallback: If VV is unavailable, times out, or returns unexpected JSON,
   construct a minimal local `VariationDescriptor` using the provided genomic
   `hgvsg`, the gene symbol, and the zygosity → GENO allelic_state mapping.

Environment flags (developer ergonomics)
----------------------------------------
P6_SKIP_VV=1           : Force the local fallback path (useful for CI/offline).
P6_ENRICH_GENE_XREFS=1 : If set and vv_lookup is importable, ask VV for HGNC/
                         Ensembl IDs and add them to gene_context where possible.
"""

from __future__ import annotations

# ---------------------------
# Imports
# ---------------------------

import os
import re
from dataclasses import dataclass
from typing import Optional, Tuple

import requests
import phenopackets.schema.v2 as pps2
from pyphetools.creation.variant_validator import VariantValidator

# vv_lookup is optional; code must work without it.
try:
    from .vv_lookup import get_gene_xrefs_vv, VVLookupError  # type: ignore
except ImportError:  # pragma: no cover
    get_gene_xrefs_vv = None  # type: ignore[assignment]

    class VVLookupError(RuntimeError):  # type: ignore[no-redef]
        """Fallback sentinel so isinstance checks compile even if vv_lookup is absent."""

        pass


# ---------------------------
# Validation patterns & enums
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

# zygosity → GENO allelic_state (numeric suffix)
_GENO_ALLELIC_STATE_CODES = {
    "heterozygous": "0000135",
    "homozygous": "0000134",
    "compound_heterozygosity": "0000191",
    "hemizygous": "0000136",
    "mosaic": "0000150",
}

# Permissive HGVS g. SNV regex (optional 'chr' prefix).
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

# Transcript + c. component, e.g., "NM_000000.0:c.100A>G" or "ENST00000205557.12:c.2428G>A"
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


# ==============================================================================
# Data model
# ==============================================================================


@dataclass
class Genotype:
    """
    A single genomic variant call plus minimal context.

    The fields mirror the genotype sheet columns; `to_variation_descriptor()`
    converts this record into a GA4GH VariationDescriptor for phenopackets.
    """

    # Raw columns (as parsed by the mapper)
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

    # --------------------------------------------------------------------------
    # Input validation (ensures downstream code can rely on sane types/shapes)
    # --------------------------------------------------------------------------
    def __post_init__(self) -> None:
        # IDs and email
        if not _VALID_ID.match(self.genotype_patient_ID):
            raise ValueError(f"Invalid patient ID: {self.genotype_patient_ID!r}")
        if not _EMAIL_PATTERN.match(self.contact_email):
            raise ValueError(f"Invalid contact email: {self.contact_email!r}")

        # Chromosome encoding (either known keyword or 'chr*')
        chrom_lower = self.chromosome.lower()
        if not (
            chrom_lower in _ALLOWED_CHROM_ENCODINGS or chrom_lower.startswith("chr")
        ):
            raise ValueError(f"Unrecognized chromosome: {self.chromosome!r}")

        # Positions must be non-negative ints
        for attr in ("start_position", "end_position"):
            val = getattr(self, attr)
            if not isinstance(val, int) or val < 0:
                raise ValueError(f"{attr} must be a non-negative integer, got {val!r}")

        # Non-empty strings for allele/gene/HGVS/protein
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

        # Controlled enums
        if self.zygosity not in _ALLOWED_ZYGOSITIES:
            raise ValueError(f"Invalid zygosity: {self.zygosity!r}")
        if self.inheritance not in _ALLOWED_INHERITANCE_MODES:
            raise ValueError(f"Invalid inheritance mode: {self.inheritance!r}")

    # --------------------------------------------------------------------------
    # Convenience
    # --------------------------------------------------------------------------
    @property
    def zygosity_code(self) -> str:
        """
        Map the human‐readable zygosity to the GENO allelic_state numeric suffix.
        """
        try:
            return _GENO_ALLELIC_STATE_CODES[self.zygosity]
        except KeyError as exc:
            raise ValueError(
                f"No GENO code defined for zygosity {self.zygosity!r}"
            ) from exc

    # --------------------------------------------------------------------------
    # Core responsibility: build a VariationDescriptor (VV path or local fallback)
    # --------------------------------------------------------------------------
    def to_variation_descriptor(self) -> "pps2.VariationDescriptor":
        """
        Build a GA4GH VariationDescriptor for this variant.

        Algorithm
        ---------
        1) If P6_SKIP_VV=1 → return a locally constructed descriptor.
        2) Else try VariantValidator with (transcript, c. part) parsed from hgvsc.
           On success, adapt & enrich; on VV errors, fall back locally.
        """
        if self._should_skip_vv():
            return self._build_local_descriptor()

        tx, c_part = self._parse_hgvsc(self.hgvsc)
        if not (tx and c_part):
            return self._build_local_descriptor()

        vd = self._try_build_descriptor_via_vv(tx, c_part)
        if vd is None:
            vd = self._build_local_descriptor()

        self._enrich_descriptor_common(vd)
        self._maybe_enrich_gene_xrefs(vd)
        return vd

    # --------------------------------------------------------------------------
    # Internal helpers (focused, testable)
    # --------------------------------------------------------------------------
    @staticmethod
    def _parse_hgvsc(hgvsc: str) -> Tuple[Optional[str], Optional[str]]:
        """
        Extract (transcript, c. part) from an hgvsc string.

        Examples
        --------
        "NM_000000.0:c.100A>G"          → ("NM_000000.0", "c.100A>G")
        "ENST00000205557.12:c.2428G>A"  → ("ENST00000205557.12", "c.2428G>A")
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
        Normalize a genomic HGVS like 'chr16:g.100A>G' → '16:g.100A>G'.

        Handles simple SNVs via regex; otherwise returns the trimmed original.
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
        if s.lower().startswith("chr"):  # best-effort strip
            return s[3:]
        return s

    def _build_local_descriptor(self) -> "pps2.VariationDescriptor":
        """
        Construct a minimal VariationDescriptor locally (VV-free).

        Adds (if available):
        - normalized genomic HGVS expression (g.)
        - allelic_state (from zygosity)
        - gene_context.symbol (from gene_symbol)
        """
        vd = pps2.VariationDescriptor()

        # expression: g.
        g_value = self._normalize_g_expression(self.hgvsg)
        if g_value:
            self._add_hgvs_expression(vd, g_value, syntax_name="HGVS")

        # allelic state
        if self.zygosity:
            vd.allelic_state.id = f"GENO:{self.zygosity_code}"
            vd.allelic_state.label = self.zygosity

        # gene symbol (optional)
        if self.gene_symbol:
            try:
                vd.gene_context.symbol = self.gene_symbol
            except AttributeError:
                pass

        return vd

    @staticmethod
    def _add_hgvs_expression(
        vd: "pps2.VariationDescriptor", value: str, *, syntax_name: str = "HGVS"
    ) -> None:
        """
        Append an Expression (value + syntax) to a VariationDescriptor.

        `syntax_name` is looked up defensively because different proto builds
        expose the enum slightly differently.
        """
        expr = vd.expressions.add()
        expr.value = value
        enum = getattr(type(expr), syntax_name, None)
        if enum is not None and hasattr(expr, "syntax"):
            expr.syntax = enum  # type: ignore[attr-defined]

    # ---------------------------
    # Complexity-splitting helpers
    # ---------------------------

    @staticmethod
    def _should_skip_vv() -> bool:
        """Return True if the P6_SKIP_VV flag is set."""
        return os.getenv("P6_SKIP_VV", "").strip() in {"1", "true", "TRUE"}

    def _try_build_descriptor_via_vv(
        self, transcript: str, c_part: str
    ) -> Optional["pps2.VariationDescriptor"]:
        """
        Attempt to build a VariationDescriptor via VariantValidator.

        Returns None on network/schema errors so caller can fall back locally.
        """
        try:
            validator = VariantValidator(genome_build="GRCh38", transcript=transcript)
            hgvs_variant = validator.encode_hgvs(c_part)  # VV call
            var_interp = hgvs_variant.to_variant_interpretation_202()
            vd = var_interp.variation_descriptor
            return vd
        except requests.RequestException:
            return None
        except (ValueError, TypeError, KeyError):
            return None

    def _enrich_descriptor_common(self, vd: "pps2.VariationDescriptor") -> None:
        """
        Add common enrichments that apply to both VV and local descriptors:
        - allelic_state (from zygosity)
        - gene_context.symbol (if missing)
        - normalized g. expression
        """
        if self.zygosity:
            vd.allelic_state.id = f"GENO:{self.zygosity_code}"
            vd.allelic_state.label = self.zygosity

        try:
            gene_ctx = getattr(vd, "gene_context", None)
            if self.gene_symbol and (
                gene_ctx is None or not getattr(gene_ctx, "symbol", "")
            ):
                vd.gene_context.symbol = self.gene_symbol
        except AttributeError:
            pass

        g_value = self._normalize_g_expression(self.hgvsg)
        if g_value:
            self._add_hgvs_expression(vd, g_value, syntax_name="HGVS")

    def _maybe_enrich_gene_xrefs(self, vd: "pps2.VariationDescriptor") -> None:
        """
        Best-effort enrichment of gene_context.id with HGNC via VV gene xref API.

        Controlled by P6_ENRICH_GENE_XREFS; never raises.
        """
        if not (
            os.getenv("P6_ENRICH_GENE_XREFS", "").strip() in {"1", "true", "TRUE"}
            and get_gene_xrefs_vv
            and self.gene_symbol
        ):
            return
        try:
            xrefs = get_gene_xrefs_vv(self.gene_symbol)
            if xrefs.get("hgnc_id") and not getattr(vd.gene_context, "id", ""):
                vd.gene_context.id = xrefs["hgnc_id"]
        except VVLookupError:
            pass
        except (requests.RequestException, ValueError, KeyError, TypeError):
            pass
