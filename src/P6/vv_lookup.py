# src/P6/vv_lookup.py
"""
VariantValidator-powered gene cross-reference lookup helpers.

This module provides a small, defensive wrapper around VariantValidator's
`gene2transcripts_v2` endpoint to retrieve:

- HGNC ID           (e.g., "HGNC:36")
- Ensembl Gene ID   (e.g., "ENSG00000101331")
- NCBI Gene ID      (e.g., "368")
- RefSeq transcripts (e.g., ["NM_001171.6", ...])
- Ensembl transcripts (e.g., ["ENST00000374403.4", ...])
- Any MANE 'select' IDs present

Typical usage
-------------
>>> from P6.vv_lookup import get_gene_xrefs_vv
>>> get_gene_xrefs_vv("ABCC6")
{
  "query": "ABCC6",
  "hgnc_id": "HGNC:36",
  "ensembl_gene": "ENSG00000101331",
  "ncbi_gene_id": "368",
  "refseq_transcripts": ["NM_001171.6", ...],
  "ensembl_transcripts": ["ENST00000374403.4", ...],
  "mane_select": ["NM_001171.6", "ENST00000374403.4"],
  "source": "https://rest.variantvalidator.org/VariantValidator/tools/..."
}

Notes
-----
- **Rate limiting**: VariantValidator limits this endpoint to ~1 request/second.
  We add a small sleep and retry on 429/5xx.
- **Resilience**: The response schema can vary slightly. We parse defensively and
  return a consistent dict. If parsing fails, a VVLookupError is raised.
- **Configuration** (environment variables):
    VV_BASE         : Base URL (default: https://rest.variantvalidator.org)
    VV_DELAY_S      : Base delay seconds between requests/retries (default: 1.0)
    VV_TIMEOUT_S    : HTTP timeout seconds (default: 15)
    VV_MAX_RETRIES  : Retries on transient errors (default: 3)

"""

from __future__ import annotations

import os
import time
import typing as t
import requests
from functools import lru_cache

# -----------------------
# Configuration & limits
# -----------------------

# VariantValidator public base; allow override for tests/self-hosting.
VV_BASE = os.getenv("VV_BASE", "https://rest.variantvalidator.org")

# Endpoint doc states ~1 req/sec; keep polite.
REQUEST_DELAY_S = float(os.getenv("VV_DELAY_S", "1.0"))

# Networking guardrails.
TIMEOUT_S = int(os.getenv("VV_TIMEOUT_S", "15"))
MAX_RETRIES = int(os.getenv("VV_MAX_RETRIES", "3"))


class VVLookupError(RuntimeError):
    """Raised when VariantValidator lookups fail or return unexpected shapes."""

    pass


def _sleep_backoff(try_idx: int) -> None:
    """
    Light linear backoff between attempts.

    Parameters
    ----------
    try_idx : int
        Zero-based attempt index (0..MAX_RETRIES-1).
    """
    time.sleep(REQUEST_DELAY_S * (try_idx + 1))


def _get_json(url: str, params: dict | None = None) -> dict:
    """
    GET JSON with retries/backoff.

    Retries on 429, 500, 502, 503, 504 and on generic exceptions. Raises
    VVLookupError after exhausting retries.

    Parameters
    ----------
    url : str
        Fully-formed URL (without `content-type`, added automatically).
    params : dict | None
        Extra query params, merged with required `content-type`.

    Returns
    -------
    dict
        Parsed JSON object from response.

    Raises
    ------
    VVLookupError
        If request or JSON parsing repeatedly fails.
    """
    last_exc: Exception | None = None
    for i in range(MAX_RETRIES):
        try:
            # VV likes ?content-type=application/json in query
            qp = {"content-type": "application/json"}
            if params:
                qp.update(params)
            resp = requests.get(url, params=qp, timeout=TIMEOUT_S)
            # transient problems: retry
            if resp.status_code in (429, 500, 502, 503, 504):
                _sleep_backoff(i)
                continue
            resp.raise_for_status()
            return resp.json()
        except Exception as e:  # network, decode, etc.
            last_exc = e
            _sleep_backoff(i)
    raise VVLookupError(f"Failed to GET {url}: {last_exc!r}")


def _normalize_transcript_id(tx: str) -> str:
    """
    Normalize a transcript ID string (currently: trim whitespace, keep version).

    Parameters
    ----------
    tx : str
        Transcript identifier, possibly including version (e.g., "NM_000088.4").

    Returns
    -------
    str
        Trimmed transcript ID.
    """
    return tx.strip()


def _pick_list(obj: t.Any) -> list[str]:
    """
    Coerce common container shapes into a list of strings.

    VariantValidator sometimes returns:
    - a simple string,
    - a list of strings,
    - or a pipe/comma-separated string.

    Parameters
    ----------
    obj : Any

    Returns
    -------
    list[str]
        List of non-empty string tokens.
    """
    if obj is None:
        return []
    if isinstance(obj, list):
        return [str(x) for x in obj]
    if isinstance(obj, str):
        if "|" in obj:
            return [s for s in obj.split("|") if s]
        if "," in obj:
            return [s for s in obj.split(",") if s]
        return [obj]
    return []


@lru_cache(maxsize=2048)
def get_gene_xrefs_vv(
    gene_query: str,
    *,
    genome_build: str = "GRCh38",
    transcript_set: str = "all",  # "all" | "refseq" | "ensembl"
    limit_transcripts: str = "mane",  # "all" | "select" | "mane" | "mane_select"
    show_exon_info: bool = False,
) -> dict:
    """
    Fetch cross-references and transcript lists for a gene via VariantValidator.

    Parameters
    ----------
    gene_query : str
        HGNC gene symbol ("ABCC6"), an HGNC:ID ("HGNC:36"), or a transcript ID.
    genome_build : str, optional
        "GRCh37" or "GRCh38" (default "GRCh38").
    transcript_set : str, optional
        "all", "refseq", or "ensembl" (default "all").
    limit_transcripts : str, optional
        "all" | "select" | "mane" | "mane_select" (default "mane").
        - "mane"      → MANE Select + Plus Clinical (per VV docs)
        - "mane_select" → only MANE Select
    show_exon_info : bool, optional
        If true, asks VV to include exon structures/alignment data (slower).

    Returns
    -------
    dict
        Normalized dictionary with (when available):
        {
          "query": str,
          "hgnc_id": str | None,
          "ensembl_gene": str | None,      # e.g. "ENSG..."
          "ncbi_gene_id": str | None,      # Entrez Gene ID as string
          "refseq_transcripts": list[str], # NM_/NR_/XM_/XR_
          "ensembl_transcripts": list[str],# ENST...
          "mane_select": list[str],        # MANE transcript IDs
          "source": str,                   # URL used
        }

    Raises
    ------
    VVLookupError
        If response retrieval or parsing fails.
    """
    url = (
        f"{VV_BASE}/VariantValidator/tools/"
        f"gene2transcripts_v2/{gene_query}/{limit_transcripts}/{transcript_set}/{genome_build}"
    )
    payload = _get_json(url, params={"show_exon_info": str(show_exon_info).lower()})

    # The response schema may be:
    # - a dict with keys like 'data', 'gene', 'genes', 'result'
    # - a list of genes
    node: t.Any = payload
    for key in ("data", "gene", "genes", "result", "response"):
        if isinstance(node, dict) and key in node:
            node = node[key]

    # If we got multiple genes, try exact symbol match, else pick first.
    if isinstance(node, list) and node:
        exact = None
        for g in node:
            sym = str(g.get("hgnc_symbol") or g.get("symbol") or "").upper()
            if sym == gene_query.upper():
                exact = g
                break
        node = exact or node[0]

    if not isinstance(node, dict):
        raise VVLookupError(
            f"Unexpected VV response shape for {gene_query}: {type(node)}"
        )

    # Helper to grab the first populated scalar among candidate keys.
    def first(*candidates: str) -> str | None:
        for c in candidates:
            v = node.get(c)
            if isinstance(v, (str, int)):
                s = str(v).strip()
                if s:
                    return s
        return None

    # IDs (favor explicit fields but remain flexible to case changes)
    hgnc_id = first("hgnc_id", "HGNC_ID", "hgnc", "HGNC")
    ensg = first("ensembl_gene_id", "ensembl_id", "ensembl_gene", "ENSG")
    entrez = first("entrez_id", "ncbi_gene_id", "NCBI_Gene_ID", "entrez")

    # Transcript collections—VV uses different key names across tools.
    refseq_txs: set[str] = set()
    ensembl_txs: set[str] = set()
    mane_txs: set[str] = set()

    # Common direct lists
    for k in ("refseq", "refseq_transcripts", "refseq_tx", "refseq_ids"):
        for tx in _pick_list(node.get(k)):
            if tx.startswith(("NM_", "NR_", "XM_", "XR_")):
                refseq_txs.add(_normalize_transcript_id(tx))

    for k in ("ensembl", "ensembl_transcripts", "ensembl_tx", "ensembl_ids"):
        for tx in _pick_list(node.get(k)):
            if tx.startswith("ENST"):
                ensembl_txs.add(_normalize_transcript_id(tx))

    for k in ("mane", "mane_select", "mane_plus_clinical"):
        for tx in _pick_list(node.get(k)):
            mane_txs.add(_normalize_transcript_id(tx))

    # Some responses include transcript objects with an 'id'/'name' field.
    for k, v in node.items():
        if "transcript" in k.lower() and isinstance(v, list):
            for item in v:
                if isinstance(item, dict):
                    tid = item.get("transcript") or item.get("id") or item.get("name")
                    if isinstance(tid, str):
                        tid = _normalize_transcript_id(tid)
                        if tid.startswith(("NM_", "NR_", "XM_", "XR_")):
                            refseq_txs.add(tid)
                        elif tid.startswith("ENST"):
                            ensembl_txs.add(tid)

    return {
        "query": gene_query,
        "hgnc_id": hgnc_id,
        "ensembl_gene": ensg,
        "ncbi_gene_id": entrez,
        "refseq_transcripts": sorted(refseq_txs),
        "ensembl_transcripts": sorted(ensembl_txs),
        "mane_select": sorted(mane_txs) if mane_txs else [],
        "source": url,
    }
