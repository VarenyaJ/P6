import glob
import os
import pytest
import re
from click.testing import CliRunner
from P6.__main__ import main


@pytest.fixture(scope="session", autouse=True)
def verify_hpo_file():
    # ensure HPO file is in place
    data_dir = os.path.join(os.path.dirname(__file__), "data")
    assert os.path.isdir(data_dir)


@pytest.mark.parametrize(
    "sample_xlsx", glob.glob(os.path.join(os.path.dirname(__file__), "data", "*.xlsx"))
)
def test_audit_excel_table_output(sample_xlsx):
    runner = CliRunner()
    result = runner.invoke(main, ["audit-excel", "-e", sample_xlsx])
    assert result.exit_code == 0, result.output

    # first line should be our header
    first = result.output.splitlines()[0]
    assert first.startswith("SHEET"), "Expected table header"
    # each subsequent line should have 4 columns
    for line in result.output.splitlines()[1:]:
        parts = re.split(r"\s{2,}", line.strip())
        assert len(parts) >= 4, f"Bad line in audit table: {line}"


def test_audit_excel_json_output(tmp_path):
    # pick any test workbook
    sample = glob.glob(os.path.join(os.path.dirname(__file__), "data", "*.xlsx"))[0]
    runner = CliRunner()
    result = runner.invoke(main, ["audit-excel", "-e", sample, "-r"])
    assert result.exit_code == 0, result.output

    # JSON must parse to a list of dicts
    import json

    payload = json.loads(result.output)
    assert isinstance(payload, list)
    assert all(isinstance(obj, dict) for obj in payload)
    # check expected keys
    for obj in payload:
        assert {"step", "sheet", "level", "message"}.issubset(obj.keys())
