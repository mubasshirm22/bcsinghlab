"""
EBI InterProScan REST adapter.

Submits a sequence to InterProScan 5 via the EBI Job Dispatcher REST API,
polls until complete, and parses the JSON result.

Confirmed valid application identifiers (from DevTools capture):
  CDD, HAMAP, NCBIfam, Panther, PrositeProfiles, PrositePatterns

API docs: https://www.ebi.ac.uk/Tools/services/rest/iprscan5
"""

import time
import requests

_BASE_URL      = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
_EMAIL         = "singhlab.notify@gmail.com"
_POLL_INTERVAL = 20    # seconds
_TIMEOUT       = 600   # 10 minutes (InterProScan is the slowest tool)

# Valid appl values confirmed from browser DevTools capture
_APPLICATIONS = [
    "CDD",
    "HAMAP",
    "NCBIfam",
    "Panther",
    "PrositeProfiles",
    "PrositePatterns",
]

# Map InterProScan entry type → our unified feature_type
_ENTRY_TYPE_MAP = {
    "DOMAIN":                 "domain",
    "FAMILY":                 "family",
    "HOMOLOGOUS_SUPERFAMILY": "domain",
    "REPEAT":                 "repeat",
    "SITE":                   "site",
    "ACTIVE_SITE":            "active_site",
    "BINDING_SITE":           "binding_site",
    "CONSERVED_SITE":         "motif",
    "PTM":                    "motif",
    "COILED_COIL":            "coiled_coil",
    "SIGNAL_PEPTIDE":         "signal_peptide",
    "TRANSMEMBRANE":          "transmembrane",
}


def run(sequence: str, job_dir: str) -> dict:
    """
    Submit sequence to InterProScan, wait for results, return unified annotations.

    Returns:
        {
          "status": "ok" | "error",
          "data": { "annotations": [ {unified annotation dict}, ... ] },
          "error": str
        }
    """
    print(f"[InterProScan] submitting ({len(sequence)} residues)…")

    # ------------------------------------------------------------------
    # Step 1 — Submit (appl sent as repeated form fields, not comma-joined)
    # ------------------------------------------------------------------
    form_data = [
        ("email",    _EMAIL),
        ("title",    "ProtPipe-job"),
        ("sequence", f">Sequence1\n{sequence}"),
        ("goterms",  "false"),
        ("pathways", "false"),
    ] + [("appl", app) for app in _APPLICATIONS]

    try:
        r = requests.post(
            f"{_BASE_URL}/run",
            data=form_data,
            timeout=30,
        )
    except requests.RequestException as e:
        return _err(f"InterProScan submission network error: {e}")

    if not r.ok:
        return _err(f"InterProScan submission HTTP {r.status_code}: {r.text[:300]}")

    job_id = r.text.strip()
    if not job_id:
        return _err(f"InterProScan: empty job ID in response: {r.text[:100]}")

    print(f"[InterProScan] job submitted: {job_id}")

    # ------------------------------------------------------------------
    # Step 2 — Poll status
    # ------------------------------------------------------------------
    start = time.time()
    while time.time() - start < _TIMEOUT:
        time.sleep(_POLL_INTERVAL)
        try:
            sr = requests.get(f"{_BASE_URL}/status/{job_id}", timeout=20)
            sr.raise_for_status()
            status = sr.text.strip().upper()
        except requests.RequestException as e:
            print(f"[InterProScan] poll error: {e}")
            continue

        if status == "RUNNING":
            elapsed = int(time.time() - start)
            print(f"[InterProScan] {job_id} still running… ({elapsed}s)")
            continue
        if status == "FINISHED":
            break
        if status in ("FAILED", "ERROR", "NOT_FOUND"):
            return _err(f"InterProScan job {job_id} ended with status: {status}")
    else:
        return _err(f"InterProScan: timed out after {_TIMEOUT // 60} minutes")

    # ------------------------------------------------------------------
    # Step 3 — Retrieve JSON result
    # ------------------------------------------------------------------
    try:
        rr = requests.get(
            f"{_BASE_URL}/result/{job_id}/json",
            headers={"Accept": "application/json"},
            timeout=60,
        )
        rr.raise_for_status()
        data = rr.json()
    except requests.RequestException as e:
        return _err(f"InterProScan result fetch error: {e}")
    except ValueError as e:
        return _err(f"InterProScan JSON parse error: {e}")

    annotations = _parse_result(data)
    print(f"[InterProScan] complete — {len(annotations)} annotations")
    return {"status": "ok", "data": {"annotations": annotations}, "error": ""}


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def _parse_result(data: dict) -> list:
    annotations = []

    # InterProScan JSON: {"results": [{"matches": [...], ...}]}
    results = data.get("results") or []
    if not results:
        return []

    matches = results[0].get("matches") or []

    for match in matches:
        signature  = match.get("signature") or {}
        entry      = signature.get("entry") or {}

        entry_acc  = entry.get("accession", "")
        entry_desc = entry.get("description", "")
        entry_type = entry.get("type", "")

        sig_acc  = signature.get("accession", "")
        sig_name = signature.get("name", "")
        sig_db   = (signature.get("signatureLibraryRelease") or {}).get("library", "")

        label        = entry_desc or sig_name or sig_acc
        accession    = entry_acc or sig_acc
        feature_type = _ENTRY_TYPE_MAP.get(entry_type.upper(), "domain")

        for loc in (match.get("locations") or []):
            start = loc.get("start")
            end   = loc.get("end")
            if start is None or end is None:
                continue

            evalue = loc.get("evalue") or match.get("evalue")
            score  = loc.get("score")  or match.get("score")

            # Quality filter
            if evalue is not None and evalue > 1e-2:
                continue

            evidence = _build_evidence(sig_db, accession, evalue, score)

            annotations.append({
                "source":       f"InterProScan/{sig_db}" if sig_db else "InterProScan",
                "feature_type": feature_type,
                "start":        int(start),
                "end":          int(end),
                "label":        (label or accession)[:60],
                "accession":    accession,
                "description":  entry_desc or sig_name,
                "e_value":      float(evalue) if evalue is not None else None,
                "score":        float(score)  if score  is not None else None,
                "evidence":     evidence,
            })

    return annotations


def _build_evidence(db: str, accession: str, evalue, score) -> str:
    parts = []
    if db:
        parts.append(db)
    if accession:
        parts.append(accession)
    if evalue is not None:
        parts.append(f"E: {evalue:.2e}")
    elif score is not None:
        parts.append(f"Score: {score:.1f}")
    return "  ".join(parts) if parts else "—"


def _err(msg: str) -> dict:
    print(f"[InterProScan] ERROR: {msg}")
    return {"status": "error", "data": {"annotations": []}, "error": msg}
