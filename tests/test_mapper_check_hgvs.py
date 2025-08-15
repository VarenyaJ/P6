"""
Tests for DefaultMapper.check_hgvs_consistency.

We check:
- matching with/without 'chr',
- BED-like SNV convention (start=pos-1, end=pos),
- error in strict mode when there's a mismatch.
"""

import pandas as pd
import hpotk
from stairval.notepad import create_notepad
from P6.mapper import DefaultMapper

HPO_PATH = "tests/data/hp.v2024-04-26.json.gz"


def make_mapper(strict=False):
    return DefaultMapper(hpotk.load_minimal_ontology(HPO_PATH), strict_variants=strict)


def test_check_hgvs_consistency_ok_bed_like_and_chr_prefix():
    """
    start=99, end=100 vs HGVS pos=100 is acceptable.
    Accept both 'chr1' and '1'.
    """
    # m = make_mapper()
    note = create_notepad("genotype")
    row = pd.Series(
        {
            "chromosome": "chr1",
            "start_position": 99,
            "end_position": 100,
            "reference": "A",
            "alternate": "G",
            "hgvsg": "1:g.100A>G",
        }
    )
    DefaultMapper.check_hgvs_consistency(row, "genotype", note, strict=False)
    assert not note.has_errors(include_subsections=True)
    assert not note.has_warnings(include_subsections=True)


def test_check_hgvs_consistency_strict_errors_on_mismatch():
    """
    With strict_variants=True we should get an error when positions disagree.
    """
    # m = make_mapper(strict=True)
    note = create_notepad("genotype")
    row = pd.Series(
        {
            "chromosome": "1",
            "start_position": 100,
            "end_position": 100,
            "reference": "A",
            "alternate": "G",
            "hgvsg": "1:g.101A>G",
        }
    )
    DefaultMapper.check_hgvs_consistency(row, "genotype", note, strict=True)
    assert note.has_errors(include_subsections=True)
