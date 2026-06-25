"""
Coiled-coil prediction via Waggawagga (MPI-BPC, Goettingen).

LUPAS (NPSA-PRABI, Lyon) was replaced because its endpoint is unreliable
(observed returning HTTP 503). Waggawagga (https://waggawagga.motorprotein.de/)
aggregates several independent coiled-coil predictors (Marcoil, Multicoil2,
Ncoils, Paircoil2 by default) and was verified live (2026-06-25).

There is no documented REST/JSON API — this is a session-based Rails app.
Protocol (reverse-engineered and verified against a live submission):

  1. GET  /                                  → grab Rails session cookie +
                                                 <meta name="csrf-token"> value
  2. POST /request/compute                    (XHR: X-Requested-With header,
       multipart form: sequence=<FASTA>,       Accept: text/javascript)
       cctool_0=1 (Marcoil), cctool_2=1 (Multicoil2),
       cctool_3=1 (Ncoils), cctool_5=1 (Paircoil2),
       authenticity_token=<csrf>
       → response body is a bare job ID string, e.g. "178236446516"
  3. GET  /request/state?job=<id>             (same session + XHR headers)
       → "0".."100" while running, "done" or "failed" when finished
  4. GET  /request/result?job=<id>            (same session + XHR headers)
       → a Rails UJS "text/javascript" payload that sets innerHTML on the
         results page with a large escaped HTML blob (chart/SVG markup for
         the web UI, not structured data).

Within that HTML blob, each tool that predicted a coiled-coil region emits a
div with an id of the form "<tool>_<tool>_<start>-<end>" (e.g.
"marcoil_marcoil_5-86", "paircoil2_paircoil2_1-88") — confirmed via live
test submissions. Tools that found nothing simply have no such div. This is
the most structured signal available without a real API, so that's what
this module scrapes. Per-residue probabilities (as LUPAS provided) are not
recoverable this way, so `probabilities` is always returned empty; region
boundaries are exact per the underlying tool's own call.
"""

import os
import re
import time
import requests

_BASE_URL = "https://waggawagga.motorprotein.de"
_TIMEOUT = 30
_POLL_STEP = 3
_POLL_WAIT = 180  # multi-tool computation can take a while

# Default sub-tools — same ones checked by default on the Waggawagga web form.
_SUBMIT_FIELDS = {
    "cctool_0": "Marcoil",
    "cctool_2": "Multicoil2",
    "cctool_3": "Ncoils",
    "cctool_5": "Paircoil2",
}

# Display-name lookup for the lowercase tool slugs scraped from result div ids.
_TOOL_DISPLAY_NAMES = {
    "marcoil": "Marcoil",
    "multicoil": "Multicoil",
    "multicoil2": "Multicoil2",
    "ncoils": "Ncoils",
    "paircoil": "Paircoil",
    "paircoil2": "Paircoil2",
}

_REGION_RE = re.compile(r'(\w+)_\1_(\d+)-(\d+)')


def run(sequence: str, job_dir: str) -> dict:
    """
    Submit to Waggawagga, scrape predicted coiled-coil regions, return
    unified annotations.

    Returns:
        {
            "status":  "ok" | "error",
            "data":    {
                "annotations":   [ {unified annotation dict}, ... ],
                "probabilities": []   # not recoverable from this service
            },
            "error":   str
        }
    """
    clean_seq = sequence.strip().upper()
    if not clean_seq:
        return _err("Empty sequence passed to coils module.")

    session = requests.Session()
    try:
        home = session.get(_BASE_URL + "/", timeout=_TIMEOUT)
        home.raise_for_status()
    except requests.RequestException as e:
        return _err(f"Waggawagga homepage request failed: {e}")

    m = re.search(r'name="csrf-token" content="([^"]+)"', home.text)
    if not m:
        return _err("Could not find CSRF token on Waggawagga homepage — site markup may have changed.")
    csrf_token = m.group(1)

    headers = {
        "X-Requested-With": "XMLHttpRequest",
        "Accept": "text/javascript, */*; q=0.01",
        "Referer": _BASE_URL + "/",
    }
    data = {
        "utf8": "✓",
        "authenticity_token": csrf_token,
        "sequence": f">query\n{clean_seq}",
        "oligomerization_scorer2": "scorer2",
        "ccwindowsize": "21",
        "commit": "Compute",
    }
    for field in _SUBMIT_FIELDS:
        data[field] = "1"

    try:
        r = session.post(_BASE_URL + "/request/compute", data=data, headers=headers, timeout=_TIMEOUT)
        r.raise_for_status()
    except requests.RequestException as e:
        return _err(f"Waggawagga submission failed: {e}")

    jobid = r.text.strip()
    if not jobid or not jobid.isdigit():
        return _err(f"Waggawagga did not return a job ID. Response (first 200 chars): {r.text[:200]!r}")

    print(f"[coils] Waggawagga jobid: {jobid}")
    deadline = time.time() + _POLL_WAIT
    state = ""
    while time.time() < deadline:
        time.sleep(_POLL_STEP)
        try:
            rs = session.get(f"{_BASE_URL}/request/state", params={"job": jobid}, headers=headers, timeout=_TIMEOUT)
            state = rs.text.strip().strip('"')
        except requests.RequestException:
            continue
        if state in ("done", "failed"):
            break

    if state not in ("done", "failed"):
        return _err(f"Waggawagga job did not finish in time (last state: {state!r}).")

    try:
        rr = session.get(f"{_BASE_URL}/request/result", params={"job": jobid}, headers=headers, timeout=_TIMEOUT)
        rr.raise_for_status()
    except requests.RequestException as e:
        return _err(f"Could not fetch Waggawagga result: {e}")

    if state == "failed":
        # Waggawagga also reports "no coiled-coil region found" as job status
        # "failed" — that's a normal negative result, not a real error. Any
        # other notice() error text is a genuine failure.
        if "no coiled coil region" in rr.text.lower():
            print("[coils] Waggawagga: no coiled-coil region found in any sub-tool")
            return {"status": "ok", "data": {"annotations": [], "probabilities": []}, "error": ""}
        msg_match = re.search(r"notice\('(.+?)'\)", rr.text)
        return _err(f"Waggawagga reported failure: {msg_match.group(1) if msg_match else rr.text[:200]}")

    try:
        with open(os.path.join(job_dir, "coils_debug.html"), "w", encoding="utf-8", errors="replace") as fh:
            fh.write(rr.text)
    except OSError:
        pass

    return _parse_result(rr.text)


# ---------------------------------------------------------------------------
# Parse Waggawagga result payload
# ---------------------------------------------------------------------------

def _parse_result(text: str) -> dict:
    annotations = []
    seen = set()
    for tool_slug, start_s, end_s in _REGION_RE.findall(text or ""):
        tool_name = _TOOL_DISPLAY_NAMES.get(tool_slug.lower())
        if not tool_name:
            continue
        start, end = int(start_s), int(end_s)
        key = (tool_name, start, end)
        if key in seen:
            continue
        seen.add(key)
        annotations.append({
            "source":           tool_name,
            "accession":        "COIL",
            "label":            "Coiled coil",
            "feature_type":     "coiled_coil",
            "start":            start,
            "end":              end,
            "score":            None,
            "e_value":          None,
            "description":      f"Predicted coiled-coil region ({tool_name}, via Waggawagga)",
            "evidence":         f"{tool_name}  {start}–{end}",
            "source_support":   [tool_name],
            "display_priority": 45,
        })

    print(f"[coils] Waggawagga: {len(annotations)} coiled-coil region(s) across sub-tools")
    return {
        "status": "ok",
        "data": {"annotations": annotations, "probabilities": []},
        "error": "",
    }


def _err(msg: str) -> dict:
    print("[coils] Waggawagga failed:", msg)
    return {"status": "error", "data": {"annotations": [], "probabilities": []}, "error": msg}
