import pytest
import pandas as pd
from P6.__main__ import preprocess
from P6.loader import load_sheets_as_tables

@pytest.fixture
def simple_workbook(tmp_path):
    # build a tiny Excel with one sheet and minimal columns
    df = pd.DataFrame({
        "contact_email": ["a@b.com"],
        "phasing": [True],
        "chromosome": ["chr1"],
        "start_position": [100],
        "end_position": [100],
        "reference": ["A"],
        "alternate": ["T"],
        "gene_symbol": ["GENE"],
        "hgvsg": ["g.100A>T"],
        "hgvsc": [""],
        "hgvsp": [""],
        "zygosity": ["het"],
        "inheritance": ["unknown"],
    }, index=["PAT1"])
    path = tmp_path / "wb.xlsx"
    with pd.ExcelWriter(path, engine="openpyxl") as w:
        df.to_excel(w, sheet_name="geno")
    return str(path)

def test_preprocess_detects_variant_sheet(simple_workbook):
    tables = load_sheets_as_tables(simple_workbook)
    entries = preprocess(tables)
    # we expect at least one classify-sheet entry
    assert any(e.step == "classify-sheet" and e.sheet == "geno" for e in entries)
    # and no variant-check errors since this is valid
    assert not any(e.step == "variant-check" and e.level == "error" for e in entries)
