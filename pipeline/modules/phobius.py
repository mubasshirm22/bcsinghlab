"""
Signal peptide and transmembrane topology prediction via Phobius.
Phobius predicts BOTH signal peptides AND TM topology in a single run.

Integration: EBI Job Dispatcher REST API
  Submit:  POST https://www.ebi.ac.uk/Tools/services/rest/phobius/run
  Status:  GET  https://www.ebi.ac.uk/Tools/services/rest/phobius/status/{jobId}
  Result:  GET  https://www.ebi.ac.uk/Tools/services/rest/phobius/result/{jobId}/out

Endpoint verified live (2026-06-25): submitting with tool="phobius" returns a
job ID and the result is reachable. The EBI wrapper only exposes result types
"out" / "sequence" / "submission" / "zip" — there is NO "short" result type,
so the previous parser (written for the standalone CLI's "-short" tabular
output) never matched anything real and silently mis-parsed "out" lines.

Actual "out" format is an EMBL-style feature table, e.g.:
  ID   EMBOSS_001
  FT   SIGNAL        1     16
  FT   DOMAIN        1      3       N-REGION.
  FT   DOMAIN        4     12       H-REGION.
  FT   DOMAIN       13     16       C-REGION.
  FT   DOMAIN       17     40       NON CYTOPLASMIC.
  //
or, for a TM protein:
  FT   TRANSMEM     45     67
  FT   DOMAIN        1     44       CYTOPLASMIC.
  FT   DOMAIN       68     90       NON CYTOPLASMIC.

Fields returned:
  - has_signal_peptide: bool
  - signal_peptide_end: int (1-based cleavage site = end of the SIGNAL region, 0 if none)
  - tm_count: int
  - topology: str (human-readable summary built from the FT lines)
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
# Parser for Phobius "out" feature-table output
# ---------------------------------------------------------------------------

_FT_LINE = re.compile(r"^FT\s+(\S+)\s+(\d+)\s+(\d+)\s*(.*)$")


def _parse_output(text: str) -> dict:
    """
    Phobius "out" output is an EMBL-style feature table:
      ID   EMBOSS_001
      FT   SIGNAL        1     16
      FT   DOMAIN        1      3       N-REGION.
      FT   TRANSMEM     45     67
      FT   DOMAIN        1     44       CYTOPLASMIC.
      //
    Exact 1-based start/end coordinates come directly from each FT line.
    """
    if "FT" not in text and "ID" not in text:
        return {
            "status": "error",
            "data": _empty_data(),
            "error": (
                "Phobius returned output but no feature-table lines were found. "
                f"Raw output (first 300 chars): {text[:300]!r}"
            ),
        }

    sp_end = 0
    tm_helices = []
    topology_parts = []

    for line in text.splitlines():
        line = line.strip()
        m = _FT_LINE.match(line)
        if not m:
            continue
        feature, start_s, end_s, desc = m.groups()
        start, end = int(start_s), int(end_s)
        desc = desc.strip().rstrip(".")

        if feature == "SIGNAL":
            sp_end = end
            topology_parts.append(f"SIGNAL {start}-{end}")
        elif feature == "TRANSMEM":
            tm_helices.append({"start": start, "end": end})
            topology_parts.append(f"TRANSMEM {start}-{end}")
        elif feature == "DOMAIN":
            topology_parts.append(f"{start}-{end} {desc}" if desc else f"DOMAIN {start}-{end}")

    has_signal_peptide = sp_end > 0
    data = {
        "has_signal_peptide": has_signal_peptide,
        "signal_peptide_end": sp_end,
        "tm_count":           len(tm_helices),
        "topology":           "; ".join(topology_parts),
        "tm_helices":         tm_helices,
    }
    print(f"[Phobius] parsed: SP={has_signal_peptide} (end={sp_end}), TM={len(tm_helices)} {tm_helices}")
    return {"status": "ok", "data": data, "error": ""}


def _empty_data() -> dict:
    return {
        "has_signal_peptide": False,
        "signal_peptide_end": 0,
        "tm_count":           0,
        "topology":           "",
        "tm_helices":         [],
    }
