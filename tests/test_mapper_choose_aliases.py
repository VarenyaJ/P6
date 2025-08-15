"""
Tests that the alias-based sheet selection works.

We pass 'variants' (→ genotype), 'hpo' (→ phenotype), and 'labs' (→ measurements)
and ensure they are recognized.
"""

import pandas as pd
import hpotk
from stairval.notepad import create_notepad
from P6.mapper import DefaultMapper

HPO_PATH = "tests/data/hp.v2024-04-26.json.gz"


def test_choose_named_tables_aliases():
    m = DefaultMapper(hpotk.load_minimal_ontology(HPO_PATH))
    note = create_notepad("alias-test")
    tables = {"variants": pd.DataFrame(), "hpo": pd.DataFrame(), "labs": pd.DataFrame()}
    selected = m._choose_named_tables(tables, note)
    assert selected.genotype is not None
    assert selected.phenotype is not None
    assert selected.measurements is not None
