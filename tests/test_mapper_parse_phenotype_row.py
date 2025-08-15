"""
Focused tests for DefaultMapper.parse_phenotype_row.

These verify:
- normal parsing of HPO ids (with/without labels),
- handling of the "NAD" placeholder,
- label/name mismatch warnings.
"""

import pandas as pd
import hpotk
from stairval.notepad import create_notepad
from P6.mapper import DefaultMapper

HPO_PATH = "tests/data/hp.v2024-04-26.json.gz"


def make_mapper():
    """Helper that builds a mapper with the test HPO graph."""
    return DefaultMapper(hpotk.load_minimal_ontology(HPO_PATH))


def test_parse_phenotype_row_ok_numeric_date_prefixed():
    """
    Numeric dates should be normalized to 'T<digits>'.
    HPO ID should be zero-padded to 7 digits.
    """
    m = make_mapper()
    note = create_notepad("phenotype")
    row = pd.Series(
        {
            "phenotype_patient_ID": "P1",
            "hpo_id": "HP:510",  # will become HP:0000510
            "date_of_observation": 20200101,
            "status": 1,  # True
        }
    )
    items, ids = m.parse_phenotype_row(row, m._hpo, "phenotype", note)
    assert len(items) == 1
    assert items[0].HPO_ID == "HP:0000510"
    assert items[0].date_of_observation == "T20200101"
    assert not note.has_errors(include_subsections=True)


def test_parse_phenotype_row_nad_skips_with_warning():
    """
    'NAD' indicates 'No Abnormality Detected' and should be skipped with a warning,
    not treated as an error.
    """
    m = make_mapper()
    note = create_notepad("phenotype")
    row = pd.Series(
        {
            "phenotype_patient_ID": "P1",
            "hpo_id": "NAD",
            "date_of_observation": "2020",
            "status": 0,
        }
    )
    items, ids = m.parse_phenotype_row(row, m._hpo, "phenotype", note)
    assert items == []
    assert note.has_warnings(include_subsections=True)
    assert not note.has_errors(include_subsections=True)


def test_parse_phenotype_row_label_mismatch_emits_warning():
    """
    If a user supplies a label that doesn't match the ontology label,
    we flag a warning but still keep the row.
    """
    m = make_mapper()
    note = create_notepad("phenotype")
    row = pd.Series(
        {
            "phenotype_patient_ID": "P1",
            "hpo_id": "Schizophrenia (HP:510)",  # label won't match real HP:0000510 name
            "date_of_observation": "T2020",
            "status": 1,
        }
    )
    m.parse_phenotype_row(row, m._hpo, "phenotype", note)
    assert note.has_warnings(include_subsections=True)
