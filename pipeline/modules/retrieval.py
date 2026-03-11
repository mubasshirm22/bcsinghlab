"""
Sequence retrieval module.

Supports three input modes:
  1. raw_fasta  — user pasted a FASTA or raw amino-acid string
  2. uniprot    — fetch by UniProt accession via official REST API
  3. ncbi       — fetch by NCBI protein accession via Entrez efetch

APIs used (all official, confirmed):
  UniProt: GET https://rest.uniprot.org/uniprotkb/{accession}.fasta  (no key needed)
  NCBI:    Biopython Bio.Entrez.efetch  (official E-utilities wrapper)  (no key needed)
"""

import re
import requests

try:
    from Bio import Entrez, SeqIO
    _BIOPYTHON = True
except ImportError:
    _BIOPYTHON = False

from pipeline.utils import fasta as fasta_utils

# NCBI requires an email for Entrez calls.  The value is informational — any
# valid-looking address works.  Set to the lab address once known.
ENTREZ_EMAIL = "singhlab.notify@gmail.com"


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run(input_data: dict, job_dir: str) -> dict:
    """
    Args:
        input_data: {
            "input_type": "raw_fasta" | "uniprot" | "ncbi",
            "sequence_input": str   # the raw text from the form
        }
        job_dir: absolute path to this job's output folder

    Returns:
        {
            "status": "ok" | "error",
            "sequence": str,        # clean uppercase sequence
            "header": str,          # description / accession
            "source": str,          # "raw" | "uniprot" | "ncbi"
            "organism": str,        # empty if unavailable
            "error": str
        }
    """
    input_type = input_data.get("input_type", "raw_fasta")
    raw = (input_data.get("sequence_input") or "").strip()

    if not raw:
        return _err("No input provided.")

    if input_type == "uniprot":
        return _fetch_uniprot(raw)
    if input_type == "ncbi":
        return _fetch_ncbi(raw)
    # Default: treat as raw FASTA / sequence
    return _parse_raw(raw)


# ---------------------------------------------------------------------------
# Raw FASTA
# ---------------------------------------------------------------------------

def _parse_raw(raw: str) -> dict:
    result = fasta_utils.parse_input(raw)
    if not result["ok"]:
        return _err(result["error"])
    return {
        "status": "ok",
        "sequence": result["sequence"],
        "header": result["header"] or "User-submitted sequence",
        "source": "raw",
        "organism": "",
        "error": "",
    }


# ---------------------------------------------------------------------------
# UniProt REST API
# Endpoint: https://rest.uniprot.org/uniprotkb/{accession}.fasta
# Method: GET   Auth: none   Confirmed: HIGH confidence
# ---------------------------------------------------------------------------

def _fetch_uniprot(accession: str) -> dict:
    accession = accession.strip().upper()
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    try:
        r = requests.get(url, timeout=20)
        if r.status_code == 404:
            return _err(f"UniProt accession not found: {accession}")
        r.raise_for_status()
    except requests.RequestException as e:
        return _err(f"UniProt fetch failed: {e}")

    text = r.text.strip()
    if not text or not text.startswith(">"):
        return _err(f"UniProt returned unexpected content for {accession}")

    parsed = fasta_utils.parse_input(text)
    if not parsed["ok"]:
        return _err(f"Could not parse UniProt FASTA: {parsed['error']}")

    # Extract organism from header  e.g. "sp|P12345|... OS=Homo sapiens OX=9606 ..."
    organism = ""
    m = re.search(r"OS=([^=]+?)(?:\s+OX=|\s+GN=|\s*$)", parsed["header"])
    if m:
        organism = m.group(1).strip()

    return {
        "status": "ok",
        "sequence": parsed["sequence"],
        "header": parsed["header"],
        "source": "uniprot",
        "organism": organism,
        "error": "",
    }


# ---------------------------------------------------------------------------
# NCBI Entrez efetch
# Uses Biopython wrapper around the official E-utilities API.
# Endpoint: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
# Method: GET/POST   Auth: none (email recommended)   Confirmed: HIGH confidence
# ---------------------------------------------------------------------------

def _fetch_ncbi(accession: str) -> dict:
    if not _BIOPYTHON:
        return _err("Biopython is not installed. Run: pip install biopython")

    Entrez.email = ENTREZ_EMAIL
    accession = accession.strip()

    try:
        handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
    except Exception as e:
        return _err(f"NCBI Entrez fetch failed for '{accession}': {e}")

    sequence = str(record.seq).upper()
    parsed = fasta_utils.parse_input(sequence)
    if not parsed["ok"]:
        return _err(f"Sequence validation failed: {parsed['error']}")

    return {
        "status": "ok",
        "sequence": parsed["sequence"],
        "header": record.description,
        "source": "ncbi",
        "organism": "",   # Entrez FASTA header doesn't always include organism
        "error": "",
    }


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _err(msg: str) -> dict:
    return {"status": "error", "sequence": "", "header": "", "source": "", "organism": "", "error": msg}
