"""
SMART (Simple Modular Architecture Research Tool) adapter.

SMART has no official public REST API. This module submits a sequence
to the SMART web form, follows the redirect to results.cgi, and extracts
domain data from embedded JavaScript JSON variables.

Endpoint (from DevTools observation):
  POST https://smart.embl.de/smart/show_motifs.pl
  Fields: ID= (empty), SEQUENCE=<amino acid sequence>
  → 302 redirect to https://smart.embl.de/results.cgi?id=<jobid>

Data extraction:
  Results are embedded as JavaScript variables in the results page:
    var pData   = {"length":N,"paths":[[{"t":type,"s":start,"e":end,...}]]}
    var visData = [{"nm":name,"st":start,"en":end,"ev":evalue,"id":id}]
    var hidData = [{"nm":name,"st":start,"en":end,"ev":evalue,"id":id,"re":reason}]
    var domNfo  = {"id":{"n":name,"e":evalue,...},...}

  visData = visible (best) hits
  hidData = hidden hits (overlap, redundancy, low quality)
"""

import re
import json
import time
import os
import requests

_SUBMIT_URL = "https://smart.embl.de/smart/show_motifs.pl"
_RESULT_URL = "https://smart.embl.de/results.cgi"
_TIMEOUT    = 120   # seconds
_POLL_MAX   = 8
_POLL_SLEEP = 10    # seconds between polls

# Map SMART type codes → unified feature_type
_SMART_TYPE_MAP = {
    "INTRINSIC":   "transmembrane",
    "SIGNAL":      "signal_peptide",
    "TM":          "transmembrane",
    "COILED":      "coiled_coil",
    "COIL":        "coiled_coil",
    "LOW_COMPLEX": "low_complexity",
    "LOW":         "low_complexity",
    "PFAM":        "domain",
    "SMART":       "domain",
}

# Keywords in domain names → feature_type
_NAME_KEYWORDS = {
    "transmembrane":  "transmembrane",
    "signal peptide": "signal_peptide",
    "coiled coil":    "coiled_coil",
    "coiled-coil":    "coiled_coil",
    "low complexity": "low_complexity",
}


def run(sequence: str, job_dir: str) -> dict:
    """
    Submit sequence to SMART and parse JavaScript-embedded result data.

    Returns:
        {
          "status": "ok" | "error" | "parse_warning",
          "data": { "annotations": [...] },
          "error": str
        }
    """
    print(f"[SMART] submitting sequence ({len(sequence)} residues)…")

    try:
        r = requests.post(
            _SUBMIT_URL,
            data={
                "ID":       "",
                "SEQUENCE": sequence,
            },
            headers={
                "User-Agent": "Mozilla/5.0 (compatible; ProtPipe/1.0; +https://singhlab.net)",
                "Referer":    "https://smart.embl.de/",
                "Accept":     "text/html,application/xhtml+xml",
            },
            timeout=_TIMEOUT,
            allow_redirects=True,
        )
    except requests.RequestException as e:
        return _err(f"SMART network error: {e}")

    if not r.ok:
        return _err(f"SMART returned HTTP {r.status_code}")

    html = r.text
    debug_path = os.path.join(job_dir, "smart_debug.html")
    _save_debug(html, debug_path)

    # Extract job ID from final URL or page content
    job_id = _extract_job_id(r.url, html)

    # Poll if results not yet ready (queue page)
    if job_id and _is_queued(html):
        print(f"[SMART] job {job_id} queued — polling…")
        for i in range(_POLL_MAX):
            time.sleep(_POLL_SLEEP)
            try:
                pr = requests.get(
                    f"{_RESULT_URL}?id={job_id}",
                    headers={"User-Agent": "Mozilla/5.0 (compatible; ProtPipe/1.0)"},
                    timeout=_TIMEOUT,
                )
                if pr.ok:
                    html = pr.text
                    _save_debug(html, debug_path)
                    if not _is_queued(html):
                        print(f"[SMART] job {job_id} ready after {(i+1)*_POLL_SLEEP}s")
                        break
            except requests.RequestException as e:
                print(f"[SMART] poll error: {e}")
        else:
            return _err(f"SMART: timed out after {_POLL_MAX * _POLL_SLEEP}s")

    annotations = _parse_smart_js(html)
    if annotations is None:
        if "no smart domains" in html.lower() or "no domains" in html.lower():
            annotations = []
        else:
            print("[SMART] WARNING: could not parse JS result data. Saving debug HTML.")
            return {
                "status": "parse_warning",
                "data":   {"annotations": []},
                "error":  "SMART parse failed — JS data format may have changed",
            }

    print(f"[SMART] complete — {len(annotations)} annotations")
    return {"status": "ok", "data": {"annotations": annotations}, "error": ""}


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _parse_smart_js(html: str) -> list | None:
    """
    Extract annotations from SMART's embedded JavaScript JSON variables.
    Returns None if the expected JS structure is not found.
    """
    vis_m   = re.search(r'var\s+visData\s*=\s*(\[.*?\]);',  html, re.DOTALL)
    hid_m   = re.search(r'var\s+hidData\s*=\s*(\[.*?\]);',  html, re.DOTALL)
    dom_m   = re.search(r'var\s+domNfo\s*=\s*(\{.*?\});',   html, re.DOTALL)
    pdata_m = re.search(r'var\s+pData\s*=\s*(\{.*?\});',    html, re.DOTALL)

    if not vis_m and not pdata_m:
        return None

    # Parse domNfo for enriched labels/accessions
    dom_info = {}
    if dom_m:
        try:
            dom_info = json.loads(dom_m.group(1))
        except Exception:
            pass

    annotations = []

    # visData — visible (best) annotations
    if vis_m:
        try:
            vis_entries = json.loads(vis_m.group(1))
        except Exception:
            vis_entries = []
        for entry in vis_entries:
            ann = _entry_to_annotation(entry, dom_info)
            if ann:
                annotations.append(ann)

    # hidData — hidden entries; include only if they add unique structural info
    if hid_m:
        try:
            hid_entries = json.loads(hid_m.group(1))
        except Exception:
            hid_entries = []
        existing_coords = {(a["start"], a["end"]) for a in annotations}
        for entry in hid_entries:
            reason = (entry.get("re") or "").lower()
            if "overlap" in reason:
                continue   # already represented by a better hit
            ann = _entry_to_annotation(entry, dom_info)
            if ann and (ann["start"], ann["end"]) not in existing_coords:
                annotations.append(ann)

    return annotations


def _entry_to_annotation(entry: dict, dom_info: dict) -> dict | None:
    """Convert a visData/hidData entry to a unified annotation dict."""
    try:
        start = int(entry.get("st", 0))
        end   = int(entry.get("en", 0))
    except (ValueError, TypeError):
        return None

    if start <= 0 or end <= 0 or end < start:
        return None

    name     = entry.get("nm", "")
    ev_raw   = entry.get("ev", "N/A")
    entry_id = str(entry.get("id", ""))

    details = dom_info.get(entry_id, {})

    evalue = None
    if ev_raw and ev_raw != "N/A":
        try:
            evalue = float(ev_raw)
        except (ValueError, TypeError):
            pass

    feature_type = _classify_entry(name, details)

    # Note: quality filtering is handled by annotation_merger, not here.
    # We return all real parsed hits so nothing is silently dropped.

    label     = (name or "")[:50]
    accession = details.get("ac", "") or _extract_accession(name)
    desc      = details.get("n", name) or name

    evidence = f"SMART  {name}" if name else "SMART"
    if evalue is not None:
        evidence += f"  E: {evalue:.2e}"

    return {
        "source":       "SMART",
        "feature_type": feature_type,
        "start":        start,
        "end":          end,
        "label":        label,
        "accession":    accession,
        "description":  desc,
        "e_value":      evalue,
        "score":        None,
        "evidence":     evidence,
    }


def _classify_entry(name: str, details: dict) -> str:
    name_lower = (name or "").lower()
    for kw, ftype in _NAME_KEYWORDS.items():
        if kw in name_lower:
            return ftype
    dtype = str(details.get("t", "") or details.get("type", "")).upper()
    for key, ftype in _SMART_TYPE_MAP.items():
        if key in dtype:
            return ftype
    return "domain"


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def _extract_job_id(url: str, html: str) -> str:
    # From redirect URL: results.cgi?id=<jobid>
    m = re.search(r"[?&]id=([A-Za-z0-9_\-]+)", url)
    if m:
        return m.group(1)
    # From page content (hidden field or link)
    m = re.search(r'[?&]id=([A-Za-z0-9_\-]{6,})', html)
    if m:
        return m.group(1)
    return ""


def _is_queued(html: str) -> bool:
    lower = html.lower()
    return "queued" in lower or "please wait" in lower or "processing" in lower


def _save_debug(html: str, path: str) -> None:
    try:
        with open(path, "w", encoding="utf-8", errors="replace") as f:
            f.write(html)
    except Exception:
        pass


def _extract_accession(name: str) -> str:
    m = re.search(r"(SM\d+|PF\d+)", name or "")
    if m:
        return m.group(1)
    return (name or "").split("(")[0].strip()


def _err(msg: str) -> dict:
    print(f"[SMART] ERROR: {msg}")
    return {"status": "error", "data": {"annotations": []}, "error": msg}
