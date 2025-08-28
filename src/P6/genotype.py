# src/P6/genotype.py
"""
Genotype domain model (P6).

P6 HIGH-LEVEL: Excel → (DefaultMapper) → domain objects (Genotype, Phenotype, …)
→ Phenopacket protobufs → JSON files.

This module provides the **Genotype** domain object and the logic to
convert a row-level variant into a GA4GH **VariationDescriptor** that
ultimately lands in a Phenopacket "interpretations" block.

INTEGRATION STRATEGY (NETWORK-AWARE):
- Preferred path: use pyphetools' VariantValidator client to normalize
  the c. description (transcript-scoped), then adapt to a VariationDescriptor.
- Resilience: if VariantValidator is disabled or unavailable, fall back to a
  locally constructed VariationDescriptor using our provided g. notation.
- Optional enrichment: when enabled, augment gene_context with HGNC/Ensembl IDs.

ENV FLAGS:
- P6_SKIP_VV=1|true     → skip VariantValidator path, always use local fallback
- P6_ENRICH_GENE_XREFS=1|true → try to enrich gene_context via vv_lookup (HGNC/ENSG)
"""

from __future__ import annotations

# =============================================================================
# [P6 SECTION] Imports & module-wide constants
# -----------------------------------------------------------------------------
# Why: Keep dependencies explicit; define regex patterns and enumerations used
#       across validation and construction paths.
# =============================================================================

import os
import re
import requests  # network exception handling—even if VV is skipped
from dataclasses import dataclass
from typing import Optional, Tuple

import phenopackets.schema.v2 as pps2
from pyphetools.creation.variant_validator import VariantValidator

# Optional: VV gene cross-ref enrichment (behind P6_ENRICH_GENE_XREFS)
try:
    from .vv_lookup import get_gene_xrefs_vv, VVLookupError  # type: ignore
except Exception:  # pragma: no cover
    get_gene_xrefs_vv = None           # type: ignore[assignment]
    VVLookupError = Exception          # type: ignore[assignment]

# --------------------------- Patterns and allowed enums -----------------------

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

# Mapping from normalized zygosity terms to GENO allelic_state codes (numeric tail)
_GENO_ALLELIC_STATE_CODES = {
    "heterozygous": "0000135",
    "homozygous": "0000134",
    "compound_heterozygosity": "0000191",
    "hemizygous": "0000136",
    "mosaic": "0000150",
}

# HGVS g. (SNV) with optional "chr" prefix: captures chromosome, position, ref, alt
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

# Transcript + c. part, e.g., "NM_000000.0:c.100A>G", "ENST00000205557.12:c.2428G>A"
_HGVSC_TXT_RE = re.compile(
    r"""
    ^\s*
    (?P<tx>
        (?:N[MR]|X[MR]|E(?:NST)?)      # NM/NR/XM/XR/ENST (and tolerant of ENST prefix)
        [_]?\d+(?:\.\d+)?              # identifier with optional dot-version
    )
    :
    (?P<c>c\..+)$
    """,
    re.IGNORECASE | re.VERBOSE,
)


# =============================================================================
# [P6 SECTION] Genotype domain object
# -----------------------------------------------------------------------------
# Why: Represents a single row-level variant call; carries raw input fields and
#       provides conversion to a GA4GH VariationDescriptor.
# =============================================================================

@dataclass
class Genotype:
    """
    A single genomic variant call for a patient (row-level abstraction).

    Attributes
    ----------
    genotype_patient_ID : str
        Unique alphanumeric patient identifier (Propagates to Phenopacket.subject.id).
    contact_email : str
        Email used for traceability or follow-up.
    phasing : bool
        Whether this call is phased (mapper-level semantics).
    chromosome : str
        Chromosome name or encoding (e.g., "chr16" or "hgvs" keyword).
    start_position : int
        1-based start coordinate (non-negative).
    end_position : int
        1-based end coordinate (non-negative).
    reference : str
        Reference allele sequence.
    alternate : str
        Alternate allele sequence.
    gene_symbol : str
        HGNC gene symbol (e.g., "BRCA1", "ABCC6").
    hgvsg : str
        Genomic HGVS (e.g., "16:g.16177614C>T" or "chr16:g.16177614C>T").
    hgvsc : str
        Transcript-scoped HGVS c. (e.g., "NM_000000.0:c.100A>G").
    hgvsp : str
        Protein HGVS p. (optional for construction; retained for display/reporting).
    zygosity : str
        One of: heterozygous, homozygous, compound_heterozygosity, hemizygous, mosaic.
    inheritance : str
        One of: unknown, inherited, de_novo_mutation.
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

    # -------------------------------------------------------------------------
    # [P6 SUBSECTION] Input validation (constructor-time)
    # Why: Fail-fast on malformed input to keep downstream construction simple.
    # -------------------------------------------------------------------------
    def __post_init__(self) -> None:
        # Patient ID
        if not _VALID_ID.match(self.genotype_patient_ID):
            raise ValueError(f"Invalid patient ID: {self.genotype_patient_ID!r}")

        # Email
        if not _EMAIL_PATTERN.match(self.contact_email):
            raise ValueError(f"Invalid contact email: {self.contact_email!r}")

        # Chromosome encoding or 'chr*'
        chrom_lower = self.chromosome.lower()
        if not (
            chrom_lower in _ALLOWED_CHROM_ENCODINGS or chrom_lower.startswith("chr")
        ):
            raise ValueError(f"Unrecognized chromosome: {self.chromosome!r}")

        # Coordinates
        for attr_name in ("start_position", "end_position"):
            value = getattr(self, attr_name)
            if not isinstance(value, int) or value < 0:
                raise ValueError(f"{attr_name} must be a non-negative integer, got {value!r}")

        # Required strings
        for attr_name in ("reference", "alternate", "gene_symbol", "hgvsg", "hgvsc", "hgvsp"):
            value = getattr(self, attr_name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{attr_name} must be a nonempty string")

        # Controlled vocabularies
        if self.zygosity not in _ALLOWED_ZYGOSITIES:
            raise ValueError(f"Invalid zygosity: {self.zygosity!r}")
        if self.inheritance not in _ALLOWED_INHERITANCE_MODES:
            raise ValueError(f"Invalid inheritance mode: {self.inheritance!r}")

    # -------------------------------------------------------------------------
    # [P6 SUBSECTION] Convenience properties (mappings & codes)
    # Why: Defer small lookups to accessors for readability in builders.
    # -------------------------------------------------------------------------
    @property
    def zygosity_code(self) -> str:
        """Return numeric portion of the GENO allelic_state code for this zygosity."""
        try:
            return _GENO_ALLELIC_STATE_CODES[self.zygosity]
        except KeyError:
            raise ValueError(f"No GENO code defined for zygosity {self.zygosity!r}")

    # =============================================================================
    # [P6 SECTION] VariationDescriptor construction
    # -----------------------------------------------------------------------------
    # Why: Provide the exact object that the Phenopacket "variantInterpretation"
    #       expects (via mapper). Prefers VariantValidator but falls back locally.
    # =============================================================================
    def to_variation_descriptor(self) -> "pps2.VariationDescriptor":
        """
        Build a GA4GH VariationDescriptor for this variant.

        PREFERRED PATH (networked):
            1) Parse transcript + c. from `hgvsc`.
            2) Call VariantValidator via pyphetools with (transcript, c_part).
            3) Convert result to VariationDescriptor.
        FALLBACK PATH (offline/resilient):
            - If P6_SKIP_VV=1|true or parsing/remote fails, synthesize a minimal
              VariationDescriptor locally using the (normalized) g. expression.

        In both paths, we:
            - attach allelic_state (GENO) from zygosity,
            - include gene_context.symbol,
            - add the g. expression as an Expression (value only; syntax HGVS).

        Optional enrichment (P6_ENRICH_GENE_XREFS=1|true):
            - If available, set gene_context.id to HGNC or ENSEMBL:ENSG…
        """
        # ---- feature flag: skip VV entirely
        if os.getenv("P6_SKIP_VV", "").strip().lower() in {"1", "true"}:
            return self._build_local_variation_descriptor()

        # ---- attempt VariantValidator path
        transcript_id, c_notation = self._parse_hgvsc(self.hgvsc)
        if not (transcript_id and c_notation):
            # no transcript+c. → fallback
            return self._build_local_variation_descriptor()

        try:
            vv_client = VariantValidator(genome_build="GRCh38", transcript=transcript_id)
            # pyphetools expects just the "c.*" portion (not "NM_...:c.*")
            vv_variant = vv_client.encode_hgvs(c_notation)
            variant_interpretation = vv_variant.to_variant_interpretation_202()
            variation_descriptor = variant_interpretation.variation_descriptor
        except requests.RequestException:
            # network/HTTP issue → fallback
            return self._build_local_variation_descriptor()
        except (ValueError, TypeError, KeyError):
            # schema/flag surprises from VV → fallback
            return self._build_local_variation_descriptor()

        # ---- post-process common fields on the VV-derived descriptor
        self._attach_common_fields(variation_descriptor)

        # ---- optional enrichment of gene_context.id (HGNC/ENSEMBL) via VV tools
        if (
            os.getenv("P6_ENRICH_GENE_XREFS", "").strip().lower() in {"1", "true"}
            and self.gene_symbol
            and get_gene_xrefs_vv is not None
        ):
            try:
                xrefs = get_gene_xrefs_vv(
                    self.gene_symbol,
                    genome_build="GRCh38",
                    transcript_set="all",
                    limit_transcripts="mane",
                )
                gene_id_curie = xrefs.get("hgnc_id") or xrefs.get("ensembl_gene")
                if gene_id_curie:
                    if gene_id_curie.startswith("ENSG"):
                        variation_descriptor.gene_context.id = f"ENSEMBL:{gene_id_curie}"
                    else:
                        variation_descriptor.gene_context.id = gene_id_curie
            except VVLookupError:
                pass  # ignore enrichment hiccups
            except Exception:
                pass  # never break the main path

        return variation_descriptor

    # -------------------------------------------------------------------------
    # [P6 SUBSECTION] Shared post-processing for VariationDescriptor
    # Why: DRY—attach allelic_state, gene_context.symbol, and g. expression.
    # -------------------------------------------------------------------------
    def _attach_common_fields(self, variation_descriptor: "pps2.VariationDescriptor") -> None:
        """Attach zygosity, gene symbol, and g. expression to the descriptor."""
        # allelic_state (GENO)
        if self.zygosity:
            variation_descriptor.allelic_state.id = f"GENO:{self.zygosity_code}"
            variation_descriptor.allelic_state.label = self.zygosity

        # gene_context.symbol (keep benign if field missing in proto)
        if self.gene_symbol:
            try:
                # fill only if absent or empty
                gene_ctx = getattr(variation_descriptor, "gene_context", None)
                if gene_ctx is None or not getattr(gene_ctx, "symbol", ""):
                    variation_descriptor.gene_context.symbol = self.gene_symbol
            except Exception:
                pass

        # include g. expression we received (normalized; strip 'chr' if present)
        g_value = self._normalize_hgvs_g(self.hgvsg)
        if g_value:
            self._add_hgvs_expression(variation_descriptor, g_value, syntax_name="HGVS")

    # -------------------------------------------------------------------------
    # [P6 SUBSECTION] Local fallback builder (no VV)
    # Why: Provide a minimal but valid VariationDescriptor independently of VV.
    # -------------------------------------------------------------------------
    def _build_local_variation_descriptor(self) -> "pps2.VariationDescriptor":
        """Create a minimal VariationDescriptor using our g. string and zygosity."""
        vd = pps2.VariationDescriptor()
        self._attach_common_fields(vd)
        return vd

    # -------------------------------------------------------------------------
    # [P6 SUBSECTION] Small internal helpers (pure functions)
    # Why: Keep parsing/normalization isolated for easier testing.
    # -------------------------------------------------------------------------
    @staticmethod
    def _parse_hgvsc(hgvsc: str) -> Tuple[Optional[str], Optional[str]]:
        """
        Extract (transcript_id, c_notation) from `hgvsc`.

        Examples
        --------
        "NM_000000.0:c.100A>G"            → ("NM_000000.0", "c.100A>G")
        "ENST00000205557.12:c.2428G>A"    → ("ENST00000205557.12", "c.2428G>A")
        """
        if not isinstance(hgvsc, str):
            return None, None
        match = _HGVSC_TXT_RE.match(hgvsc.strip())
        if not match:
            return None, None
        return match.group("tx"), match.group("c")

    @staticmethod
    def _normalize_hgvs_g(hgvsg: str) -> Optional[str]:
        """
        Normalize genomic HGVS to a compact form:
            "chr16:g.100A>G" → "16:g.100A>G" (SNV shape only).
        Falls back to returning a 'chr'-stripped value or original string.
        """
        if not isinstance(hgvsg, str) or not hgvsg.strip():
            return None
        s = hgvsg.strip()
        match = _HGVS_G_SNV.match(s)
        if match:
            chrom = match.group("chrom")
            pos = match.group("pos")
            ref = match.group("ref").upper()
            alt = match.group("alt").upper()
            return f"{chrom}:g.{pos}{ref}>{alt}"
        if s.lower().startswith("chr"):
            return s[3:]
        return s

    @staticmethod
    def _add_hgvs_expression(
        variation_descriptor: "pps2.VariationDescriptor",
        value: str,
        *,
        syntax_name: str = "HGVS",
    ) -> None:
        """
        Append an Expression(value, syntax) to a VariationDescriptor.

        Parameters
        ----------
        variation_descriptor : VariationDescriptor
            Descriptor to mutate.
        value : str
            HGVS string to attach (e.g., "16:g.100A>G").
        syntax_name : str, optional
            Enum name to set on `expr.syntax` if available (default: "HGVS").
        """
        expr = variation_descriptor.expressions.add()
        expr.value = value
        # Some proto builds present Expression.HGVS; guard with getattr.
        enum = getattr(type(expr), syntax_name, None)
        if enum is not None and hasattr(expr, "syntax"):
            expr.syntax = enum  # type: ignore[attr-defined]
