"""
Signal peptide and transmembrane topology prediction via Phobius.
Phobius predicts BOTH signal peptides AND TM topology in a single run.

Integration: EBI Job Dispatcher REST API
  Submit:  POST https://www.ebi.ac.uk/Tools/services/rest/phobius/run
  Status:  GET  https://www.ebi.ac.uk/Tools/services/rest/phobius/status/{jobId}
  Result:  GET  https://www.ebi.ac.uk/Tools/services/rest/phobius/result/{jobId}/out

⚠️  ENDPOINT VERIFICATION REQUIRED
Run the following curl and paste the response:

  curl -s -X POST \\
    "https://www.ebi.ac.uk/Tools/services/rest/phobius/run" \\
    -d "email=test@test.com&sequence=MKTIIALSYIFCLVFA&stype=protein" \\
    -D -

If the submit URL returns 404, the tool name in the EBI dispatcher path may differ
(e.g. "phobius" vs "pfa_phobius"). Report the response and I'll fix it.

Expected output format (Phobius "short" output):
  ID  SEQ_LENGTH  TM  SP  TOPOLOGY
  seq1  150  0  Y  n25c                ← signal peptide ending at residue 25
  seq1  200  2  0  i5-27o45-67i        ← 2 TM helices, no signal peptide

Fields returned:
  - has_signal_peptide: bool
  - signal_peptide_end: int (1-based cleavage site, 0 if none)
  - tm_count: int
  - topology: str (raw topology string, e.g. "i5-27o45-67i")
  - tm_helices: list of {start: int, end: int}  (1-based, empty if none)
"""

import re
from pipeline.utils.ebi_poll import submit_and_wait

_TOOL  = "phobius"
_EMAIL = "singhlab.notify@gmail.com"


def run(sequence: str, job_dir: str) -> dict:
    """
    Returns:
        {
          "status": "ok" | "error" | "endpoint_unverified",
          "data": {
              "has_signal_peptide": bool,
              "signal_peptide_end": int,
              "tm_count": int,
              "topology": str,
              "tm_helices": [{"start": int, "end": int}, ...]
          },
          "error": str
        }
    """
    print(f"[Phobius] submitting ({len(sequence)} residues)…")

    params = {
        "email":    _EMAIL,
        "sequence": f">query\n{sequence}",
        "stype":    "protein",
        "format":   "short",
    }

    result = submit_and_wait(
        tool=_TOOL,
        params=params,
        result_type="out",
        poll_interval=15,
        timeout=300,
        label="Phobius",
    )

    if not result["ok"]:
        # Check if this is a 404 / endpoint issue vs a real failure
        if "submission failed" in result["error"].lower():
            return {
                "status": "endpoint_unverified",
                "data": _empty_data(),
                "error": (
                    "Phobius EBI endpoint could not be reached. "
                    "Please run the curl test from the implementation plan and report the HTTP status code."
                ),
            }
        return {"status": "error", "data": _empty_data(), "error": result["error"]}

    return _parse_output(result["text"])


# ---------------------------------------------------------------------------
# Parser for Phobius "short" output format
# ---------------------------------------------------------------------------

def _parse_output(text: str) -> dict:
    """
    Phobius short output looks like:
      ID             SEQLEN TM  SP  PREDICTION
      query             150  0   Y  n25c
    or:
      query             200  2   0  i5-27o45-67i
    """
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("ID") or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 5:
            continue

        # Parts: [id, seqlen, tm_count, sp, topology]
        try:
            tm_count = int(parts[2])
        except ValueError:
            tm_count = 0

        sp_field = parts[3]
        has_signal_peptide = (sp_field.strip() not in ("0", "N", "NO", ""))
        topology = parts[4] if len(parts) > 4 else ""

        # Extract SP cleavage site from topology, e.g. "n25c" → end=25
        sp_end = 0
        if has_signal_peptide:
            m = re.match(r"n(\d+)c", topology, re.IGNORECASE)
            if m:
                sp_end = int(m.group(1))

        # Parse TM helices from topology string e.g. "i5-27o45-67i"
        # Pattern: digits-digits optionally preceded by i/o/n/c
        tm_helices = []
        for match in re.finditer(r"(\d+)-(\d+)", topology):
            start, end = int(match.group(1)), int(match.group(2))
            # Skip the signal peptide region if present
            if has_signal_peptide and end <= sp_end:
                continue
            tm_helices.append({"start": start, "end": end})

        data = {
            "has_signal_peptide": has_signal_peptide,
            "signal_peptide_end": sp_end,
            "tm_count":           tm_count,
            "topology":           topology,
            "tm_helices":         tm_helices,
        }
        print(f"[Phobius] parsed: SP={has_signal_peptide}, TM={tm_count}")
        return {"status": "ok", "data": data, "error": ""}

    # If we couldn't parse any result line, return a clear error rather than silently
    # returning zeros — the caller can decide whether to surface this to the user.
    return {
        "status": "error",
        "data": _empty_data(),
        "error": (
            "Phobius returned output but no parseable result lines were found. "
            f"Raw output (first 300 chars): {text[:300]!r}"
        ),
    }


def _empty_data() -> dict:
    return {
        "has_signal_peptide": False,
        "signal_peptide_end": 0,
        "tm_count":           0,
        "topology":           "",
        "tm_helices":         [],
    }
