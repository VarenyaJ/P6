"""
Genotype domain model.

Defines the Genotype class which encapsulates all relevant fields for a
genomic variant entry, with validation of each attribute.
"""

import re
from dataclasses import dataclass
import phenopackets.schema.v2 as pps2
# import pyphetools
from google.protobuf.any_pb2 import Any

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

# permissive HGVS g. pattern with optional "chr" prefix
# Captures: chromosome, position, ref, alt (SNVs)
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

    # Variation construction

    def to_variation_descriptor(self) -> "pps2.VariationDescriptor":
        """
        Build a GA4GH VariationDescriptor for this genotype.

        - Adds HGVS expressions (g, c, p).
        - Sets gene_context and allelic_state.
        - Stores a manual VRS Allele JSON into the description field for visibility.
        - Adds molecule_context and stable id if possible.
        """
        import json, hashlib

        vd = pps2.VariationDescriptor()

        # ---- gene_context ----
        try:
            vd.gene_context.symbol = self.gene_symbol
        except Exception:
            pass

        # ---- allelic_state ----
        try:
            vd.allelic_state.id = f"GENO:{self.zygosity_code}"
            vd.allelic_state.label = self.zygosity
        except Exception:
            pass

        # ---- expressions ----
        g_expr = self._canonicalize_g_hgvs(self.hgvsg)
        g_expr = self._attempt_pyphetools_normalize(g_expr) or g_expr

        def _add_expr(value: str):
            expr = vd.expressions.add()
            expr.value = value
            if hasattr(type(expr), "HGVS"):
                try:
                    expr.syntax = type(expr).HGVS
                except Exception:
                    pass

        if g_expr.strip():
            _add_expr(g_expr)
        if self.hgvsc.strip():
            _add_expr(self.hgvsc.strip())
        if self.hgvsp.strip():
            _add_expr(self.hgvsp.strip())

        # ---- variation info (serialize as JSON string into description) ----
        try:
            m = _HGVS_G_SNV.match(g_expr)
            if m:
                pos_1b = int(m.group("pos"))
                allele_json = {
                    "@type": "ga4gh.vrs.v1.Allele",
                    "location": {
                        "@type": "SequenceLocation",
                        "sequenceId": self._to_accession(m.group("chrom")),
                        "interval": {
                            "@type": "SequenceInterval",
                            "start": pos_1b - 1,
                            "end": pos_1b,
                        },
                    },
                    "state": {
                        "@type": "LiteralSequenceExpression",
                        "sequence": m.group("alt").upper(),
                    },
                }
            else:
                allele_json = {
                    "@type": "ga4gh.vrs.v1.Allele",
                    "location": {
                        "@type": "SequenceLocation",
                        "sequenceId": self._to_accession(self.chromosome),
                        "interval": {
                            "@type": "SequenceInterval",
                            "start": int(self.start_position) - 1,
                            "end": int(self.end_position),
                        },
                    },
                    "state": {
                        "@type": "LiteralSequenceExpression",
                        "sequence": str(self.alternate).upper(),
                    },
                }

            vd.description = json.dumps(allele_json)
        except Exception:
            pass

        # ---- molecule context ----
        try:
            vd.molecule_context = vd.MoleculeContext.DNA
        except Exception:
            pass

        # ---- stable id ----
        try:
            vd.id = "vd:" + hashlib.sha1(g_expr.encode()).hexdigest()[:16]
        except Exception:
            pass

        return vd


    # --- helper for chromosome → accession mapping ---
    @staticmethod
    def _to_accession(chrom: str) -> str:
        """Map chromosome names to GRCh38 RefSeq accessions; fallback to input string."""
        table = {
            "1": "NC_000001.11",
            "2": "NC_000002.12",
            "3": "NC_000003.12",
            "16": "NC_000016.10",
            "X": "NC_000023.11",
            "Y": "NC_000024.10",
            "MT": "NC_012920.1",
        }
        c = chrom.strip().lower().removeprefix("chr")
        c = "MT" if c in ("m", "mt") else c.upper()
        return table.get(c, chrom)


    # -----------------------
    # Small internal helpers
    # -----------------------

    @staticmethod
    def _canonicalize_g_hgvs(hgvs_g: str) -> str:
        """
        Remove optional 'chr' prefix from genomic HGVS expressions.

        Parameters
        ----------
        hgvs_g : str
            Genomic HGVS string (e.g., 'chr16:g.100A>G' or '16:g.100A>G').

        Returns
        -------
        str
            Canonicalized HGVS (e.g., '16:g.100A>G').
        """
        s = (hgvs_g or "").strip()
        if s.lower().startswith("chr"):
            return s[3:]
        return s

    @staticmethod
    def _add_hgvs_expression(
            vd: "pps2.VariationDescriptor", value: str, syntax_name: str = "HGVS"
    ) -> None:
        """
        Append an Expression to a VariationDescriptor with best-effort syntax.

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
        # Many proto builds expose `Expression.HGVS`; guard with try/except.
        try:
            enum = getattr(type(expr), syntax_name)
            expr.syntax = enum  # type: ignore[attr-defined]
        except Exception:
            # Older builds may not have the enum; value alone still helps.
            pass

    @staticmethod
    def _attempt_pyphetools_normalize(hgvs_g: str) -> str | None:
        """
        Best-effort normalization using pyphetools if a suitable API exists.

        This is intentionally defensive: pyphetools' public surface for direct
        HGVS normalization isn’t guaranteed. If we cannot import a normalizer,
        we quietly return None and keep the original HGVS string.

        Parameters
        ----------
        hgvs_g : str
            Input genomic HGVS expression.

        Returns
        -------
        Optional[str]
            A normalized HGVS string if available, else None.
        """
        try:
            # NOTE:
            # pyphetools primarily exposes normalization inside its TemplateImporter
            # flow (via VariantValidator). If, in a future release, a direct API
            # becomes public (e.g., pyphetools.hgvs.normalize), wire it here.
              # noqa: F401
            # Placeholder: no stable, documented function to call here yet.
            # If we decide to add a thin wrapper later, plug it in:
            #   from pyphetools.hgvs import normalize_g
            #   return normalize_g(hgvs_g)  # hypothetical
            import pyphetools  # noqa: F401
            # Placeholder: no stable, documented function to call here yet.
            return None
        except Exception:
            return None