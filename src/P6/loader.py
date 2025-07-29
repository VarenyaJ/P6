import pandas as pd

# Columns that need renaming → target dataclass fields
RENAME_MAP = {
    # genotype columns
    "ref": "reference",
    "alt": "alternate",
    "gene": "gene_symbol",
    "start": "start_position",
    "end": "end_position",
    "chrom": "chromosome",
    # phenotype columns
    "hpo": "hpo_id",
    "timestamp": "date_of_observation",
}

def load_sheets_as_tables(workbook_path: str) -> dict[str, pd.DataFrame]:
    """
    Read each worksheet into a DataFrame:
      - first row = header
      - first column = index
      - normalize all headers to snake_case lowercase
      - apply renames from RENAME_MAP
    """
    
    excel = pd.ExcelFile(workbook_path, engine="openpyxl")
    tables: dict[str, pd.DataFrame] = {}

    for sheet_name in excel.sheet_names:
        df = pd.read_excel(
            excel, sheet_name=sheet_name, header=0, index_col=0, engine="openpyxl"
        )

        # CLEAN & NORMALIZE headers:
        df.columns = (
            df.columns.str.strip()
            .str.replace(r"\s*\(.*?\)", "", regex=True)  # drop any "(…)"
            .str.replace(r"\s+", "_", regex=True)  # spaces → underscore
            .str.replace(":", "", regex=False)  # drop colons
            .str.lower()
        )

        # apply specific renames (e.g. "ref" → "reference")
        df = df.rename(
            columns={
                orig: target
                for orig, target in RENAME_MAP.items()
                if orig in df.columns
            }
        )

        tables[sheet_name] = df

    return tables