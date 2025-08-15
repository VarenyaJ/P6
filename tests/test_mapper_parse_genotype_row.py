# tests/test_mapper_parse_genotype_row.py
"""
Unit tests for DefaultMapper.parse_genotype_row.

These focus narrowly on:
- tokenized zygosity/inheritance (e.g., "het/hom" + "inherited/denovo")
- graceful handling of unknown codes (should emit an error and return [])
- missing contact_email (should default to 'unknown@example.com')

We don’t need the ontology here because genotype parsing doesn’t consult it.
"""

import pandas as pd

from P6.mapper import DefaultMapper
from stairval.notepad import create_notepad

# Small helpers for clarity


def make_mapper() -> DefaultMapper:
    """Create a DefaultMapper with a dummy ontology (not used for genotype parsing)."""
    # We can pass None for the ontology because parse_genotype_row doesn't access it.
    return DefaultMapper(hpo=None, strict_variants=False)


def base_row() -> dict:
    """
    A minimal, syntactically valid genotype row.
    Individual tests override zygosity / inheritance (and other fields) as needed.
    """
    return {
        "genotype_patient_ID": "P100",
        "contact_email": "user@example.com",
        "phasing": 1,
        "chromosome": "16",
        "start_position": 100,
        "end_position": 100,
        "reference": "A",
        "alternate": "G",
        "gene_symbol": "GENE1",
        # Allow optional 'chr' prefix; mapper.check_hgvs_consistency tolerates it.
        "hgvsg": "chr16:g.100A>G",
        # Minimal but syntactically valid c./p. strings (keeps row construction happy).
        "hgvsc": "NM_000000.0:c.100A>G",
        "hgvsp": "NP_000000.0:p.(Lys34Glu)",
        # Defaults (most tests override these)
        "zygosity": "het",
        "inheritance": "inherited",
    }


# The actual tests


def test_parse_genotype_row_multi_tokens_creates_two_records():
    """
    zygosity 'het/hom' and inheritance 'inherited/denovo' should produce 2 Genotype objects,
    pairing tokens positionally (het↔inherited, hom↔denovo).
    """
    m = make_mapper()
    note = create_notepad("genotype")
    row = pd.Series(
        {**base_row(), "zygosity": "het/hom", "inheritance": "inherited/denovo"}
    )

    items, aux = m.parse_genotype_row(row, "genotype", note)
