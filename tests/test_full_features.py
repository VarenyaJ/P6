import os
import re
import glob
import json
from click.testing import CliRunner
from P6.__main__ import main

# Path to the “full” example workbook (e.g. created under tests/data/full_example.xlsx)
FULL_XLSX = os.path.join(
    os.path.dirname(__file__), "data", "Sydney_Python_transformation.xlsx"
)


def test_full_features_parse_creates_all_blocks(decompressed_hpo):
    runner = CliRunner()
    result = runner.invoke(
        main, ["parse-excel", "-e", FULL_XLSX, "--custom-hpo", decompressed_hpo]
    )
    assert result.exit_code == 0, result.output

    # Extract the output directory from the summary line
    m = re.search(r"Wrote \d+ phenopacket files to (.+)", result.output)
    assert m, "Did not find output directory in CLI summary"
    outdir = m.group(1).strip()

    # There should be at least one JSON file
    json_files = glob.glob(os.path.join(outdir, "*.json"))
    assert json_files, f"No JSON files written to {outdir}"

    # Load the first phenopacket and check new keys
    pkt = json.load(open(json_files[0], encoding="utf-8"))
    assert isinstance(pkt.get("diseases", []), list), "Missing 'diseases' block"
    assert isinstance(pkt.get("measurements", []), list), "Missing 'measurements' block"
    assert isinstance(pkt.get("biosamples", []), list), "Missing 'biosamples' block"
