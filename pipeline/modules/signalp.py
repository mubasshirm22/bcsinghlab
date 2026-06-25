"""
Signal peptide prediction via SignalP-5.0 (DTU Health Tech).

The EBI Job Dispatcher does NOT host a "signalp" tool (verified live —
submitting returns "Tool 'signalp' was not found"), so this module talks
directly to DTU's own submission gateway instead, which is confirmed live:
https://services.healthtech.dtu.dk/services/SignalP-5.0/

Integration: DTU webface2.cgi job gateway (same system used by NetSurfP-2.0,
see services/netsurf.py for the sibling implementation).
  Submit:  POST https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi
           multipart: uploadfile=<fasta>, configfile=<SignalP-5.0 webface.cf>
           → redirects to a URL containing ?jobid=<JOBID>
  Poll:    GET  https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi?ajax=1&jobid=<JOBID>&wait=20
           → JSON job-queue status dict; wait for "status": "finished"
  Result:  GET  https://services.healthtech.dtu.dk/services/SignalP-5.0/tmp/<JOBID>/output.json

Verified live (2026-06-25) — example output.json for a sequence with a signal peptide:
  {
    "SEQUENCES": {
      "query": {
        "CS_pos": "Cleavage site between pos. 16 and 17: VFA-DY. Probability: 0.6551",
        "Likelihood": [0.791, 0.209],
        "Prediction": "Signal peptide (Sec/SPI)",
        "Protein_types": ["Signal Peptide (Sec/SPI)", "Other"]
      }
    }
  }

Fields returned:
  - has_signal_peptide: bool
  - signal_peptide_start: int (always 1 if present)
  - signal_peptide_end: int (1-based, last residue of the signal peptide / cleavage site)
  - probability: float (cleavage-site probability, 0-1)
  - prediction_label: str (e.g. "Signal peptide (Sec/SPI)" or "Other")
"""

import io
import re
import time
import json
import requests

_WEBFACE_URL = "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi"
_TMP_BASE_URL = "https://services.healthtech.dtu.dk/services/SignalP-5.0/tmp/"
_CONFIG_FILE = "/var/www/services/services/SignalP-5.0/webface.cf"
_POLL_SLEEP = 15
_CANCEL_AFTER = 600  # 10 minutes — SignalP is fast, usually finishes in <30s

_CS_POS_RE = re.compile(r"between pos\.\s*(\d+)\s*and\s*(\d+).*?Probability:\s*([\d.]+)", re.IGNORECASE)


def run(sequence: str, job_dir: str) -> dict:
    print(f"[SignalP] submitting ({len(sequence)} residues) to DTU SignalP-5.0…")

    fasta_bytes = f">query\n{sequence}\n".encode()
    files = {"uploadfile": ("query.fasta", io.BytesIO(fasta_bytes), "text/plain")}
    data = {"configfile": _CONFIG_FILE}

    try:
        r = requests.post(_WEBFACE_URL, data=data, files=files, allow_redirects=True, timeout=60)
    except requests.RequestException as e:
        return {"status": "error", "data": _empty_data(), "error": f"Network error during submission: {e}"}

    jobid = _extract_jobid(r)
    if not jobid:
        return {
            "status": "error",
            "data": _empty_data(),
            "error": f"Could not obtain job ID from SignalP server. Final URL: {r.url}",
        }

    print(f"[SignalP] jobid: {jobid}")
    output, error = _poll_and_fetch(jobid)
    if output is None:
        return {"status": "error", "data": _empty_data(), "error": error or "SignalP job timed out."}

    return _parse_output(output)


# ---------------------------------------------------------------------------
# Submission / polling
# ---------------------------------------------------------------------------

def _extract_jobid(response) -> str:
    m = re.search(r"[?&]jobid=([A-Za-z0-9_\-]+)", response.url)
    if m:
        return m.group(1)
    m = re.search(r"[?&]jobid=([A-Za-z0-9_\-]+)", response.text)
    return m.group(1) if m else ""


def _poll_and_fetch(jobid: str):
    deadline = time.time() + _CANCEL_AFTER
    while time.time() < deadline:
        time.sleep(_POLL_SLEEP)
        try:
            r = requests.get(
                _WEBFACE_URL,
                params={"ajax": "1", "jobid": jobid, "wait": "20"},
                timeout=60,
                allow_redirects=True,
            )
        except requests.RequestException as e:
            print("[SignalP] poll error:", e)
            continue

        text = r.text.strip()
        if not (text.startswith("{") or text.startswith("[")):
            continue
        try:
            job = json.loads(text)
        except json.JSONDecodeError:
            continue

        status = str(job.get("status", "")).lower()
        if status == "finished":
            return _fetch_output_json(jobid)
        if status in ("error", "failed"):
            return None, f"SignalP job reported status: {status}"

    return None, "SignalP job timed out or returned no result"


def _fetch_output_json(jobid: str):
    url = _TMP_BASE_URL + jobid + "/output.json"
    try:
        r = requests.get(url, timeout=60, allow_redirects=True)
        r.raise_for_status()
        return r.json(), ""
    except (requests.RequestException, json.JSONDecodeError) as e:
        return None, f"Could not fetch SignalP output.json: {e}"


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def _parse_output(output: dict) -> dict:
    sequences = (output or {}).get("SEQUENCES") or {}
    if not sequences:
        return {
            "status": "error",
            "data": _empty_data(),
            "error": f"SignalP output.json had no SEQUENCES entry: {str(output)[:300]!r}",
        }

    seq_result = next(iter(sequences.values()))
    prediction = str(seq_result.get("Prediction", ""))
    has_sp = "signal peptide" in prediction.lower()

    sp_end = 0
    probability = 0.0
    cs_pos = str(seq_result.get("CS_pos", ""))
    m = _CS_POS_RE.search(cs_pos)
    if m:
        sp_end = int(m.group(1))
        probability = float(m.group(3))
    elif has_sp:
        likelihood = seq_result.get("Likelihood") or []
        probability = float(likelihood[0]) if likelihood else 0.0

    data = {
        "has_signal_peptide":   has_sp,
        "signal_peptide_start": 1 if has_sp else 0,
        "signal_peptide_end":   sp_end,
        "probability":          probability,
        "prediction_label":     prediction,
    }
    print(f"[SignalP] parsed: SP={has_sp}, end={sp_end}, prob={probability:.3f}")
    return {"status": "ok", "data": data, "error": ""}


def _empty_data() -> dict:
    return {
        "has_signal_peptide":   False,
        "signal_peptide_start": 0,
        "signal_peptide_end":   0,
        "probability":          0.0,
        "prediction_label":     "",
    }
