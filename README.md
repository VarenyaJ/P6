# P6
**Peter's Parse and Processing of Prenatal Particulars via Pandas**

A simple, extensible CLI for downloading the Human Phenotype Ontology, parsing genotype/phenotype Excel workbooks, and producing [GA4GH Phenopackets](https://phenopacket-schema.readthedocs.io/en/latest/schema.html#version-2-0) as specified [here](https://phenopacket-schema.readthedocs.io/_/downloads/en/stable/pdf/).

## Table of Contents

1. [Features](#features)  
2. [Prerequisites](#prerequisites)  
3. [Installation](#installation)  
4. [Quickstart](#quickstart)  
   - [Download HPO JSON](#download-hpo-json)  
   - [Parse Excel to Phenopackets](#parse-excel-to-phenopackets)  
5. [CLI Reference](#cli-reference)  
   - [`p6 download`](#p6-download)  
   - [`p6 parse-excel`](#p6-parse-excel)  
6. [Development & Testing](#development--testing)  
7. [Contributing](#contributing)  
8. [License](#license)  
9. [Contact](#contact)

## Features

- **Download**: fetch the latest or a specific `hp.json` release from GitHub  
- **Parse**: autodetect genotype vs phenotype sheets in any Excel workbook  
- **Normalize**: clean up column names, HPO IDs, timestamps, and data types  
- **Generate**: emit individual Phenopacket files, one per record (will change the file extension later)

## Installation

1.  **Clone** the repo:  
    ```bash
    git clone https://github.com/VarenyaJ/P6.git
    cd P6
    ```

2.  (Recommended) Create a virtual environment (venv or Conda):
    # === Simple Venv setup ===
    ```bash
    python3 -m venv .venv
    source .venv/bin/activate
    ```
    # === or with Conda ===
    ```bash
    conda env create -f requirements/environment.yml -y
    conda activate P6
    ```

3.  Install via pip:
    ```bash
    python3 -m pip install -r requirements/requirements.txt .
    ```

4.  Verify the installation:
    ```bash
    p6 --help
    ```
    
    You should see something like:
    ```bash
    Usage: p6 [OPTIONS] COMMAND [ARGS]...
    
      P6: Peter's Parse and Processing of Prenatal Particulars via Pandas.
    
    Options:
      --help  Show this message and exit.
    
    Commands:
      download    Download a specific or the latest HPO JSON release into...
      parse-excel Read each sheet, check column order, then: - Identify as a...
    ```

## Quickstart

### Download HPO JSON

Fetch the latest release into `tests/data/` (the default directory):
```bash
p6 download
```

After running, you’ll have `tests/data/hp.json`.


### Parse Excel to Phenopackets

With your HPO JSON in place at `tests/data/hp.json`, run:
```bash
p6 parse-excel -e tests/data/Sydney_Python_transformation.xlsx
```

Resulting phenopacket files will be under:
```plaintext
phenopacket-from-excel/$(date "+%Y-%m-%d_%H-%M-%S")/phenopackets/
```

## CLI Reference

### p6 download

Usage:
```markdown
p6 download [OPTIONS]

Options:
    -d, --data-path PATH    where to save HPO JSON (default: tests/data)
    -v, --hpo-version TEXT  exact HPO release tag (e.g. 2025-03-03 or v2025-03-03)
    --help                  Show this help message and exit.
```

Examples:

Fetch a specific release tag (e.g. v2025-03-03 or 2025-03-03) into `tests/data/` (the default directory):
```bash
p6 download -v 2025-03-03
p6 download --hpo-version 2025-03-03
```

Fetch a specific release tag (e.g. v2025-03-03 or 2025-03-03) into a custom directory:
```bash
p6 download -d src/P6 -v 2025-03-03
p6 download --data-path src/P6 --hpo-version 2025-03-03
```

### p6 parse-excel

Read an Excel workbook, classify sheets, normalize fields, and emit Phenopacket protobuffers.

Usage: `p6 parse-excel [OPTIONS] EXCEL_FILE`

Options:
```markdown
    -e, --excel-path FILE    path to the Excel workbook  [required]
    -hpo, --custom-hpo FILE  path to a custom HPO JSON file (defaults to `tests/data/hp.json`)
    --help                  Show this message and exit.
```

Example:

Explicitly point at a custom HPO file:
```bash
p6 parse-excel -e tests/data/Sydney_Python_transformation.xlsx -hpo src/P6/hp.json
```

## Development & Testing

Install dev requirements:
```bash
python3 -m pip install -r requirements/requirements.txt -r requirements/requirements_test.txt .
```
This will install `P6` along with the dependencies needed for the development.

Run the full test suite:
```bash
pytest -q
```

Lint & type-check (via ruff and built-in assertions):
```bash
ruff check .
ruff format .
```

## Contributing

1. Fork the repo & create a feature branch
2. Make your changes & add tests
3. Ensure all tests pass & lint is clean
4. Submit a pull request against main
5. Please follow the AGPL-3.0 code of conduct.

## License

This project is licensed under the AGPL-3.0. See LICENSE for details.

## Contact

Varenya Jain
varenyajj@gmail.com
GitHub: @VarenyaJ
