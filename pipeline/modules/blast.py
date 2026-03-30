"""
Homology search via NCBI remote BLAST.

Uses Biopython's NCBIWWW wrapper around the official NCBI BLAST URL API.
  Submission endpoint: https://blast.ncbi.nlm.nih.gov/Blast.cgi  CMD=Put
  Polling endpoint:    https://blast.ncbi.nlm.nih.gov/Blast.cgi  CMD=Get
  Auth: none (email+tool params recommended for identification)
  Confirmed: HIGH confidence

Database: swissprot (curated, fast, ~200 MB on NCBI side — no local storage needed)
Alternative: "nr" for broader search but much slower (10–30+ min)

⚠️  This call is SLOW (typically 2–10 minutes for a 150-residue sequence against
    swissprot). It runs in a background thread inside runner.py, so the web UI
    is not blocked.
"""

import re

try:
    from Bio.Blast import NCBIWWW, NCBIXML
    from Bio import Entrez as _Entrez
    _BIOPYTHON = True
except ImportError:
    _BIOPYTHON = False

_MAX_HITS = 10
_DATABASE = "swissprot"
_PROGRAM  = "blastp"
_EMAIL    = "singhlab.notify@gmail.com"
_TOOL     = "singhlab_protpipe"


def run(sequence: str, job_dir: str, database: str = _DATABASE, max_hits: int = _MAX_HITS) -> dict:
    """
    Run blastp against swissprot (or specified database) and return top hits.

    Returns:
        {
          "status": "ok" | "error",
          "data": {
              "database": str,
              "hits": [
                  {
                    "accession": str,
                    "title": str,
                    "length": int,
                    "e_value": float,
                    "score": float,
                    "identity_pct": float,
                    "coverage_pct": float,
                    "query_start": int,
                    "query_end": int,
                  }, ...
              ]
          },
          "error": str
        }
    """
    if not _BIOPYTHON:
        return {"status": "error", "data": {}, "error": "Biopython not installed (pip install biopython)"}

    # Set identification params on Entrez module — required by newer Biopython (1.79+)
    # qblast() removed its own email/tool kwargs in favour of this global setting
    _Entrez.email = _EMAIL
    _Entrez.tool  = _TOOL

    print(f"[BLAST] submitting blastp to {database} ({len(sequence)} residues)…")

    try:
        result_handle = NCBIWWW.qblast(
            program=_PROGRAM,
            database=database,
            sequence=sequence,
            hitlist_size=max_hits,
            format_type="XML",
        )
    except Exception as e:
        return {"status": "error", "data": {}, "error": f"BLAST submission failed: {e}"}

    try:
        blast_record = NCBIXML.read(result_handle)
        result_handle.close()
    except Exception as e:
        return {"status": "error", "data": {}, "error": f"BLAST XML parse failed: {e}"}

    hits = []
    seq_len = len(sequence)

    for alignment in blast_record.alignments:
        if not alignment.hsps:
            continue
        hsp = alignment.hsps[0]  # best HSP for this alignment

        align_len = hsp.align_length if hsp.align_length else 1
        identity_pct = round(hsp.identities / align_len * 100, 1)
        coverage_pct  = round((hsp.query_end - hsp.query_start + 1) / seq_len * 100, 1)

        protein_name, organism = _parse_title(alignment.title)
        hits.append({
            "accession":    alignment.accession,
            "title":        alignment.title[:120],   # full raw title for reference
            "protein_name": protein_name,
            "organism":     organism,
            "length":       alignment.length,
            "e_value":      hsp.expect,
            "score":        hsp.score,
            "identity_pct": identity_pct,
            "coverage_pct": coverage_pct,
            "query_start":  hsp.query_start,
            "query_end":    hsp.query_end,
        })

    print(f"[BLAST] complete — {len(hits)} hits returned")
    return {
        "status": "ok",
        "data": {
            "database": database,
            "hits": hits,
        },
        "error": "",
    }


def _parse_title(title: str) -> tuple:
    """
    Extract (protein_name, organism) from a BLAST hit title.

    Handles formats:
      sp|P08473|NEP_HUMAN RecName: Full=Neprilysin; Flags: Precursor [Homo sapiens]
      sp|W4VS99.1| RecName: Full=Neprilysin-1; Flags: Precursor [Trittame loki]
      gi|12345|ref|NP_001234.1| Hypothetical protein [Mus musculus]
      Simple description without prefix [Organism]

    Returns:
        (protein_name: str, organism: str)
        Both stripped; empty string if not found.
    """
    # Strip leading db|acc| prefixes (sp|ACC|GENE  or  gb|ACC|  etc.)
    # These look like:  sp|P08473|NEP_HUMAN  or  sp|W4VS99.1|
    stripped = re.sub(r'^(?:[a-z]{2,4}\|[^\|]+\|[^\s]*\s*)+', '', title, flags=re.IGNORECASE).strip()
    if not stripped:
        stripped = title

    # Extract organism from [...] at the end
    organism = ""
    org_m = re.search(r'\[([^\[\]]+)\]\s*$', stripped)
    if org_m:
        organism = org_m.group(1).strip()
        stripped = stripped[:org_m.start()].strip()

    # Extract protein name from RecName: Full=...
    recname_m = re.search(r'RecName:\s*Full=([^;{[]+)', stripped, re.IGNORECASE)
    if recname_m:
        protein_name = recname_m.group(1).strip().rstrip(';').strip()
    else:
        # Remove flag-like suffixes (Flags: ..., SubName: ..., AltName: ...)
        clean = re.sub(r'(?:Flags|SubName|AltName|Contains|Includes)\s*:.*', '', stripped, flags=re.IGNORECASE).strip()
        protein_name = clean.rstrip(';').strip() or stripped.strip()

    # Cap length for display
    if len(protein_name) > 80:
        protein_name = protein_name[:77] + "…"

    return protein_name, organism
