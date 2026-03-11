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

        hits.append({
            "accession":    alignment.accession,
            "title":        alignment.title[:120],   # truncate for display
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
