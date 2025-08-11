"""
Ensure the table-level wrappers enforce required columns and report errors.

Each *_table method should return [] and emit an error if required columns are missing.
"""

import pandas as pd
import hpotk
from stairval.notepad import create_notepad
from P6.mapper import DefaultMapper

HPO_PATH = "tests/data/hp.v2024-04-26.json.gz"


def make_mapper():
    return DefaultMapper(hpotk.load_minimal_ontology(HPO_PATH))


def test_map_genotype_table_missing_required_columns_errors():
    m = make_mapper()
    note = create_notepad("genotype")
    # Intentionally missing many required genotype columns
    df = pd.DataFrame({"contact_email": ["a@b.com"]}).set_index(pd.Index(["P1"]))
    records = m._map_genotype_table(df, note)
    assert records == []
    assert note.has_errors(include_subsections=True)


def test_map_phenotype_table_missing_required_columns_errors():
    m = make_mapper()
    note = create_notepad("phenotype")
    df = pd.DataFrame({"hpo_id": ["HP:1"]}).set_index(pd.Index(["P1"]))
    records = m._map_phenotype_table(df, note)
    assert records == []
    assert note.has_errors(include_subsections=True)


def test_map_diseases_table_missing_required_columns_errors():
    m = make_mapper()
    note = create_notepad("diseases")
    df = pd.DataFrame({"patient_ID": ["P1"], "disease_term": ["MONDO:000"]})
    records = m._map_diseases_table(df, note)
    assert records == []
    assert note.has_errors(include_subsections=True)


def test_map_measurements_table_missing_required_columns_errors():
    m = make_mapper()
    note = create_notepad("measurements")
    df = pd.DataFrame({"patient_ID": ["P1"], "measurement_type": ["X"]})
    records = m._map_measurements_table(df, note)
    assert records == []
    assert note.has_errors(include_subsections=True)


def test_map_biosamples_table_missing_required_columns_errors():
    m = make_mapper()
    note = create_notepad("biosamples")
    df = pd.DataFrame({"patient_ID": ["P1"], "biosample_id": ["B1"]})
    records = m._map_biosamples_table(df, note)
    assert records == []
    assert note.has_errors(include_subsections=True)
