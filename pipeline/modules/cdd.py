"""
NCBI CDD (Conserved Domain Database) adapter via Batch Web CD-Search.

API: https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi
  Step 1 — POST FASTA, get back a CDSID (job token)
  Step 2 — Poll GET until status == 0 (complete)
  Step 3 — Retrieve tab-delimited hit table

Official docs: https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#BatchWebCDSearch

Result format (tab-separated columns):
  Session   RID   Query#  Hit#  PSSM_ID  From  To  E-Value  Bitscore
  Name  Incomplete  Superfamily  (17 cols total)

dmode=rep returns representative non-overlapping domains only —
the cleanest output for architectural display.
"""

import re
import time
import requests

_BASE_URL      = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"
_POLL_INTERVAL = 10    # seconds between poll attempts
_TIMEOUT       = 300   # 5 minutes max


def run(sequence: str, job_dir: str) -> dict:
    """
    Run CDD search against the sequence.

    Returns:
        {
          "status": "ok" | "error",
          "data": {
              "annotations": [
                  {
                    "source": "CDD",
                    "feature_type": str,
                    "start": int,
                    "end": int,
                    "label": str,
                    "accession": str,
                    "description": str,
                    "e_value": float,
                    "score": float,
                    "evidence": str,
                  }, ...
              ]
          },
          "error": str
        }
    """
    fasta = f">query\n{sequence}"
    print(f"[CDD] submitting Batch CD-Search ({len(sequence)} residues)…")

    # ------------------------------------------------------------------
    # Step 1 — Submit
    # ------------------------------------------------------------------
    try:
        r = requests.post(
            _BASE_URL,
            data={
                "queries":        fasta,
                "db":             "cdd",
                "evalue":         "0.01",
                "maxhit":         "500",
                "useid1":         "T",
                "compbasedadj":   "T",
                "filter":         "true",
                "dmode":          "rep",   # representative (non-overlapping) domains
                "tdata":          "hits",
            },
            timeout=30,
        )
    except requests.RequestException as e:
        return _err(f"CDD submission network error: {e}")

    if not r.ok:
        return _err(f"CDD submission HTTP {r.status_code}: {r.text[:200]}")

    cdsid = _extract_cdsid(r.text)
    if not cdsid:
        return _err(f"CDD: could not extract CDSID from response:\n{r.text[:400]}")

    print(f"[CDD] job submitted: {cdsid}")

    # ------------------------------------------------------------------
    # Step 2 — Poll
    # ------------------------------------------------------------------
    start = time.time()
    while time.time() - start < _TIMEOUT:
        time.sleep(_POLL_INTERVAL)
        try:
            pr = requests.get(
                _BASE_URL,
                params={"cdsid": cdsid, "tdata": "hits", "cddefl": "true"},
                timeout=30,
            )
            pr.raise_for_status()
        except requests.RequestException as e:
            print(f"[CDD] poll error: {e}")
            continue

        status_code = _extract_status(pr.text)
        if status_code == 3:   # still running
            print(f"[CDD] {cdsid} still running…")
            continue
        if status_code != 0:   # error
            return _err(f"CDD job failed with status code {status_code}")
        # status 0 = complete
        annotations = _parse_hits(pr.text)
        print(f"[CDD] complete — {len(annotations)} hits")
        return {"status": "ok", "data": {"annotations": annotations}, "error": ""}

    return _err(f"CDD: timed out after {_TIMEOUT // 60} minutes")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _extract_cdsid(text: str) -> str:
    """Extract CDSID from the submission response."""
    # Response contains a line like: #cdsid   QM3-qcdsearch-XXXXXXXXXXXXXXXX-XXXXXXXXXX
    m = re.search(r"#cdsid\s+(\S+)", text)
    if m:
        return m.group(1)
    # Fallback: any token that looks like a CDD session ID
    m = re.search(r"(QM3-qcdsearch-\S+)", text)
    if m:
        return m.group(1)
    return ""


def _extract_status(text: str) -> int:
    """
    Extract status code from CDD poll response.
    #status 0  = complete
    #status 3  = running
    #status 1,2,4,5,6 = error states
    """
    m = re.search(r"#status\s+(\d+)", text)
    if m:
        return int(m.group(1))
    # If the response contains actual data lines (not just headers/status), treat as complete
    for line in text.splitlines():
        if line and not line.startswith("#") and "\t" in line:
            return 0
    return 3  # assume still running if we can't parse


def _parse_hits(text: str) -> list:
    """
    Parse CDD tab-delimited hit table.

    Column indices (0-based) from NCBI docs:
      0  Session
      1  RID (run ID)
      2  Query#
      3  Hit#
      4  PSSM_ID
      5  From
      6  To
      7  E-Value
      8  Bitscore
      9  Accession
      10 Short_name
      11 Incomplete
      12 Superfamily

    Some responses have slight column variations; we locate From/To/E-Value
    by header inspection first.
    """
    annotations = []
    header_cols = None

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            continue

        parts = line.split("\t")

        # Detect the column header row (does not start with #)
        if "From" in parts and "To" in parts and "E-Value" in parts:
            header_cols = [c.strip() for c in parts]
            continue

        if not header_cols:
            continue

        # Use header to locate columns
        if len(parts) < len(header_cols):
            continue
        try:
            idx        = {c: i for i, c in enumerate(header_cols)}
            start_val  = parts[idx["From"]].strip()
            end_val    = parts[idx["To"]].strip()
            evalue_str = parts[idx["E-Value"]].strip()
            accession  = parts[idx["Accession"]].strip()
            sn_col     = idx.get("Short_name", idx.get("Short name"))
            short_name = parts[sn_col].strip() if sn_col is not None else ""
            incomp     = parts[idx["Incomplete"]].strip() if "Incomplete" in idx else ""
            superfam   = parts[idx["Superfamily"]].strip() if "Superfamily" in idx else ""
        except (IndexError, KeyError):
            continue

        try:
            start  = int(start_val)
            end    = int(end_val)
            evalue = float(evalue_str)
        except ValueError:
            continue

        # Filter: superfamily-only hits are less specific; include but mark lower priority
        is_superfamily = bool(superfam and superfam != accession)

        # Skip hits that are pure superfamily calls with poor e-value
        if is_superfamily and evalue > 1e-3:
            continue

        # Quality threshold for specific domain hits
        if not is_superfamily and evalue > 1e-2:
            continue

        feature_type = "domain"  # CDD is always domain-level
        label = short_name or accession
        evidence = f"E-value: {_fmt_evalue(evalue)}"
        if incomp:
            evidence += f"  ({incomp})"

        annotations.append({
            "source":       "CDD",
            "feature_type": feature_type,
            "start":        start,
            "end":          end,
            "label":        label,
            "accession":    accession,
            "description":  short_name,
            "e_value":      evalue,
            "score":        None,
            "evidence":     evidence,
        })

    return annotations


def _fmt_evalue(e: float) -> str:
    if e == 0.0:
        return "0.0"
    if e < 1e-4:
        exp = int(f"{e:.0e}".split("e")[1])
        mant = e / (10 ** exp)
        return f"{mant:.2f} × 10^{exp}"
    return f"{e:.2e}"


def _err(msg: str) -> dict:
    print(f"[CDD] ERROR: {msg}")
    return {"status": "error", "data": {"annotations": []}, "error": msg}
