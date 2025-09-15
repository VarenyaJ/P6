"""
VariantValidator cross-reference helpers.

High level (P6 perspective)
---------------------------
This module is an optional enrichment layer used by P6 to attach gene
cross-references (HGNC ID, Ensembl gene ID, and canonical transcript
accessions) to the VariationDescriptor's gene_context, *after* the core
variant normalization path has succeeded.

Key behaviors
-------------
- Calls VariantValidator's `gene2transcripts` endpoints (v2 preferred, v1 fallback).
- Normalizes the sometimes-variable VV payloads into a small, stable dict:
  {
      "hgnc_id": "HGNC:####",
      "ensembl_gene_id": "ENSG###########",
      "refseq_transcripts": [...],
      "ensembl_transcripts": [...]
  }
- Uses small retry/backoff for resilience.
- Is deliberately decoupled from the main parsing flow: any failure raises
  `VVLookupError`; the caller should catch and ignore if enrichment isn't critical.

Environment
-----------
VV_BASE_URL   : Optional base URL override (default "https://rest.variantvalidator.org")
"""

from __future__ import annotations

from functools import lru_cache
from typing import Any, Dict, List
import json
import os
import time
from urllib.parse import quote as _urlencode

import requests


class VVLookupError(RuntimeError):
    """Raised when VariantValidator enrichment lookups fail."""


# ------------------------------------------------------------------------------
# Module configuration
# ------------------------------------------------------------------------------

_VV_BASE = os.getenv("VV_BASE_URL", "https://rest.variantvalidator.org").rstrip("/")


# ------------------------------------------------------------------------------
# Small utilities
# ------------------------------------------------------------------------------


def _sleep_backoff(i: int) -> None:
    """
    Sleep using a small exponential backoff (polite to the VV API).
    Sequence ~ 0.25s, 0.5s, 1s, 2s.
    """
    time.sleep(0.25 * (2**i))


def _request_json(url: str, *, timeout: float = 10.0) -> dict:
    """
    GET JSON with simple retry/backoff for VV endpoints.

    Retries a few times on network/HTTP/JSON decode problems and raises
    VVLookupError if all attempts fail.
    """
    last_exc: Exception | None = None
    for i in range(4):  # attempts: 0,1,2,3
        try:
            resp = requests.get(url, timeout=timeout)
            resp.raise_for_status()
            return resp.json()
        except (requests.RequestException, json.JSONDecodeError, ValueError) as e:
            last_exc = e
            _sleep_backoff(i)
    assert last_exc is not None
    raise VVLookupError(f"Failed GET {url}: {last_exc}") from last_exc


# ------------------------------------------------------------------------------
# Payload normalizers (keep public API stable even if VV changes shape)
# ------------------------------------------------------------------------------


def _parse_v2_payload(payload: Dict[str, Any]) -> Dict[str, Any]:
    """
    Normalize the gene2transcripts_v2 response into a compact dict.

    Expected output keys:
      - hgnc_id (str)
      - ensembl_gene_id (str)
      - refseq_transcripts (List[str])
      - ensembl_transcripts (List[str])
    """
    out: Dict[str, Any] = {
        "hgnc_id": "",
        "ensembl_gene_id": "",
        "refseq_transcripts": [],
        "ensembl_transcripts": [],
    }
    if not isinstance(payload, dict):
        return out

    hgnc = payload.get("hgnc", {})
    if isinstance(hgnc, dict):
        out["hgnc_id"] = hgnc.get("hgnc_id", "") or hgnc.get("HGNC_ID", "")
        out["ensembl_gene_id"] = hgnc.get("ensembl_gene_id", "") or hgnc.get(
            "ensembl", ""
        )

    def _collect(lst: Any) -> List[str]:
        accs: List[str] = []
        if isinstance(lst, list):
            for item in lst:
                if isinstance(item, dict) and item.get("accession"):
                    accs.append(str(item["accession"]))
        return accs

    out["refseq_transcripts"] = _collect(payload.get("refseq"))
    out["ensembl_transcripts"] = _collect(payload.get("ensembl"))
    return out


def _parse_v1_payload(payload: Dict[str, Any]) -> Dict[str, Any]:
    """
    Normalize the gene2transcripts (v1) response into a compact dict.

    v1 is simpler; transcripts are often plain string lists.
    """
    out: Dict[str, Any] = {
        "hgnc_id": "",
        "ensembl_gene_id": "",
        "refseq_transcripts": [],
        "ensembl_transcripts": [],
    }
    if not isinstance(payload, dict):
        return out

    out["hgnc_id"] = payload.get("hgnc_id", "") or payload.get("HGNC_ID", "")
    out["ensembl_gene_id"] = payload.get("ENSEMBL", "") or payload.get("ensembl", "")

    rs = payload.get("refseq") or payload.get("RefSeq") or []
    if isinstance(rs, list):
        out["refseq_transcripts"] = [str(r) for r in rs if isinstance(r, str)]

    es = payload.get("ensembl_transcripts") or payload.get("ensembl") or []
    if isinstance(es, list):
        out["ensembl_transcripts"] = [str(e) for e in es if isinstance(e, str)]

    return out


# ------------------------------------------------------------------------------
# Public API
# ------------------------------------------------------------------------------


@lru_cache(maxsize=2048)
def get_gene_xrefs_vv(
    gene_query: str,
    *,
    genome_build: str = "GRCh38",
    transcript_set: str = "refseq",  # {refseq, ensembl, all}
    limit_transcripts: str = "mane",  # {mane, mane_select, select, raw, all}
) -> Dict[str, Any]:
    """
    Fetch HGNC, Ensembl, and transcript xrefs for an HGNC symbol/ID or transcript.

    Parameters
    ----------
    gene_query : str
        HGNC symbol/ID (e.g., 'ABCC6' or 'HGNC:36') or a transcript ID (NM_/ENST_).
    genome_build : str, optional
        'GRCh37' or 'GRCh38' (default 'GRCh38').
    transcript_set : str, optional
        'refseq', 'ensembl', or 'all' (default 'refseq').
    limit_transcripts : str, optional
        VV's limiting switch (default 'mane'). Useful values: 'mane', 'mane_select', 'select'.

    Returns
    -------
    dict
        Compact dict with keys: 'hgnc_id', 'ensembl_gene_id',
        'refseq_transcripts', 'ensembl_transcripts'.

    Raises
    ------
    VVLookupError
        If VV is unreachable or returns an unparseable/empty payload.
    """
    if not gene_query or not isinstance(gene_query, str):
        raise VVLookupError("gene_query must be a non-empty string")
    gene_query = gene_query.strip()

    # Preferred v2 endpoint
    v2_url = (
        f"{_VV_BASE}/VariantValidator/tools/gene2transcripts_v2/"
        f"{_urlencode(gene_query)}/{_urlencode(limit_transcripts)}/"
        f"{_urlencode(transcript_set)}/{_urlencode(genome_build)}"
        "?content-type=application%2Fjson"
    )
    v2_data = _request_json(v2_url)
    v2_norm = _parse_v2_payload(v2_data)
    if any(v2_norm.values()):
        return v2_norm

    # Fallback to the simpler v1 endpoint
    v1_url = (
        f"{_VV_BASE}/VariantValidator/tools/gene2transcripts/"
        f"{_urlencode(gene_query)}?content-type=application%2Fjson"
    )
    v1_data = _request_json(v1_url)
    v1_norm = _parse_v1_payload(v1_data)
    if any(v1_norm.values()):
        return v1_norm

    raise VVLookupError(f"No xrefs found for {gene_query!r}")
