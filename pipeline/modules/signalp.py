"""
Signal peptide prediction via SignalP (EBI Job Dispatcher).

Integration: EBI Job Dispatcher REST API
  Submit:  POST https://www.ebi.ac.uk/Tools/services/rest/signalp/run
  Status:  GET  https://www.ebi.ac.uk/Tools/services/rest/signalp/status/{jobId}
  Result:  GET  https://www.ebi.ac.uk/Tools/services/rest/signalp/result/{jobId}/out

Parameters:
  email:    contact email
  sequence: FASTA-formatted protein sequence
  organism: euk  (eukaryote — default; also gram+, gram-)
  format:   short

Expected output format (SignalP-4.1 short):
  # SignalP-4.1 euk predictions
  # name         Cmax  pos  Ymax  pos  Smax  pos  Smean   D      ?  Dmaxcut  Networks-used
  query          0.763   19  0.748   19  0.871   10  0.804  0.776  Y  0.450  SignalP-noTM

Column key (space-split, skip comment lines):
  [0]  name
  [1]  Cmax    — max cleavage-site score
  [2]  cmax_pos — position of cleavage site (last residue of signal peptide, 1-based)
  [3]  Ymax
  [4]  ymax_pos
  [5]  Smax
  [6]  smax_pos
  [7]  Smean
  [8]  D       — discriminant score
  [9]  ?       — Y = signal peptide detected
  [10] Dmaxcut
"""

from pipeline.utils.ebi_poll import submit_and_wait

_TOOL  = "signalp"
_EMAIL = "singhlab.notify@gmail.com"


def run(sequence: str, job_dir: str) -> dict:
    """
    Returns:
        {
          "status":  "ok" | "error" | "endpoint_unverified",
          "data":    {
              "has_signal_peptide": bool,
              "signal_peptide_start": 1,
              "signal_peptide_end":   int,   # last residue of SP (cmax_pos)
              "d_score":              float,
              "cmax_score":           float,
          },
          "error":   str
        }
    """
    print(f"[SignalP] submitting ({len(sequence)} residues)…")

    params = {
        "email":    _EMAIL,
        "sequence": f">query\n{sequence}",
        "organism": "euk",
        "format":   "short",
    }

    result = submit_and_wait(
        tool=_TOOL,
        params=params,
        result_type="out",
        poll_interval=15,
        timeout=300,
        label="SignalP",
    )

    if not result["ok"]:
        if "submission failed" in result["error"].lower():
            return {
                "status": "endpoint_unverified",
                "data": _empty_data(),
                "error": (
                    "SignalP EBI endpoint could not be reached. "
                    "The EBI Job Dispatcher may not currently host SignalP."
                ),
            }
        return {"status": "error", "data": _empty_data(), "error": result["error"]}

    return _parse_output(result["text"])


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def _parse_output(text: str) -> dict:
    """
    Parse SignalP-4.1 short output.
    Skip lines starting with '#'. First non-comment line is the result.
    """
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        fields = line.split()
        if len(fields) < 10:
            continue

        try:
            cmax_pos  = int(fields[2])
            d_score   = float(fields[8])
            cmax_score = float(fields[1])
            has_sp    = (fields[9].upper() == "Y")
        except (ValueError, IndexError):
            continue

        sp_end = cmax_pos if has_sp else 0
        data = {
            "has_signal_peptide":   has_sp,
            "signal_peptide_start": 1 if has_sp else 0,
            "signal_peptide_end":   sp_end,
            "d_score":              d_score,
            "cmax_score":           cmax_score,
        }
        print(f"[SignalP] parsed: SP={has_sp}, end={sp_end}, D={d_score:.3f}")
        return {"status": "ok", "data": data, "error": ""}

    return {
        "status": "error",
        "data":   _empty_data(),
        "error":  (
            "SignalP returned output but no parseable result lines were found. "
            f"Raw (first 300 chars): {text[:300]!r}"
        ),
    }


def _empty_data() -> dict:
    return {
        "has_signal_peptide":   False,
        "signal_peptide_start": 0,
        "signal_peptide_end":   0,
        "d_score":              0.0,
        "cmax_score":           0.0,
    }
