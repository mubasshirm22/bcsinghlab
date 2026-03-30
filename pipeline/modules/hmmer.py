"""
Pfam domain annotation via EBI HMMER web API.

⚠️  ENDPOINT VERIFICATION REQUIRED
The EBI HMMER REST API endpoint below is based on documented sources
(confidence: MEDIUM-HIGH) but MUST be verified with a live curl test before
relying on it in production:

  curl -s -X POST \\
    "https://www.ebi.ac.uk/Tools/hmmer/api/v1/search/hmmscan" \\
    -H "Content-Type: application/x-www-form-urlencoded" \\
    -d "input=MKTIIALSYIFCLVFA&database=pfam" \\
    -D -

If this returns 404 or the response shape differs, paste the actual response
and I will update the parser accordingly.

Expected submission response (JSON):
  {"id": "<uuid>", ...}

Expected result endpoint (GET):
  https://www.ebi.ac.uk/Tools/hmmer/api/v1/result/<uuid>

Expected result JSON shape (based on HMMER web server documentation):
  {
    "results": {
      "hits": [
        {
          "name": "PF00001.27",
          "acc":  "PF00001",
          "desc": "7 transmembrane receptor (rhodopsin family)",
          "pvalue": 1.2e-10,
          "evalue": 5.4e-7,
          "domains": [
            {"ienv": 10, "jenv": 300, "iali": 15, "jali": 295, "ievalue": 1e-8}
          ]
        }
      ]
    }
  }
"""

import time
import requests

_SUBMIT_URL  = "https://www.ebi.ac.uk/Tools/hmmer/api/v1/search/hmmscan"
_RESULT_BASE = "https://www.ebi.ac.uk/Tools/hmmer/api/v1/result"
_DATABASE    = "pfam"
_POLL_INTERVAL = 10   # seconds
_TIMEOUT       = 300  # 5 minutes is usually plenty for HMMER


def run(sequence: str, job_dir: str) -> dict:
    """
    Search sequence against Pfam via EBI HMMER API.

    Returns:
        {
          "status": "ok" | "error" | "endpoint_unverified",
          "data": {
              "domains": [
                  {
                    "name": str,        # Pfam ID e.g. PF00001
                    "description": str,
                    "e_value": float,
                    "seq_start": int,   # 1-based, on the query sequence
                    "seq_end": int,
                  }, ...
              ]
          },
          "error": str
        }
    """
    print(f"[HMMER] submitting hmmscan to Pfam ({len(sequence)} residues)…")

    # Step 1 — submit
    # The EBI HMMER API expects a JSON body (not URL-encoded form data).
    # Sending form data produces HTTP 400 "Cannot parse request body".
    try:
        r = requests.post(
            _SUBMIT_URL,
            json={"input": sequence, "database": _DATABASE},
            headers={"Accept": "application/json"},
            timeout=30,
        )
    except requests.RequestException as e:
        return _err(f"HMMER submission network error: {e}")

    if r.status_code == 404:
        return {
            "status": "endpoint_unverified",
            "data": {"domains": []},
            "error": (
                "EBI HMMER API returned 404. The endpoint URL may have changed. "
                "Please run the curl test from the implementation plan and report the response."
            ),
        }

    if not r.ok:
        return _err(f"HMMER submission returned HTTP {r.status_code}: {r.text[:200]}")

    try:
        job_info = r.json()
        job_uuid = job_info.get("id") or job_info.get("uuid") or job_info.get("job_id")
        if not job_uuid:
            # Some versions return the UUID as the direct response body
            job_uuid = r.text.strip()
        if not job_uuid:
            return _err(f"HMMER: could not extract job ID from response: {r.text[:300]}")
    except Exception:
        # Try treating the raw text as the job ID (some APIs do this)
        job_uuid = r.text.strip()
        if not job_uuid:
            return _err(f"HMMER: unparseable submission response: {r.text[:300]}")

    print(f"[HMMER] job submitted: {job_uuid}")

    # Step 2 — poll
    start = time.time()
    result_data = None
    while time.time() - start < _TIMEOUT:
        time.sleep(_POLL_INTERVAL)
        try:
            pr = requests.get(f"{_RESULT_BASE}/{job_uuid}",
                              headers={"Accept": "application/json"},
                              timeout=20)
            if pr.status_code == 202:  # still running
                print(f"[HMMER] {job_uuid} still running…")
                continue
            pr.raise_for_status()
            result_data = pr.json()
            break
        except requests.RequestException as e:
            print(f"[HMMER] poll error: {e}")
            continue

    if result_data is None:
        return _err(f"HMMER: timed out after {_TIMEOUT // 60} minutes")

    # Step 3 — parse
    return _parse_result(result_data)


# ---------------------------------------------------------------------------
# Parser — update here if the actual JSON shape differs from expectation
# ---------------------------------------------------------------------------

def _parse_result(data: dict) -> dict:
    domains = []

    # Navigate into the result structure — adjust path if actual shape differs
    hits = (
        data.get("result", {}).get("hits")
        or data.get("results", {}).get("hits")
        or data.get("hits")
        or []
    )

    for hit in hits:
        if not hit.get("is_included", True):
            continue
        acc = hit.get("acc") or hit.get("name") or "unknown"
        # Strip version suffix: "PF00049.23" → "PF00049"
        name = acc.split(".")[0] if "." in acc else acc
        desc = hit.get("desc") or hit.get("description") or ""
        evalue = hit.get("evalue") or hit.get("pvalue") or 1.0

        for dom in (hit.get("domains") or []):
            if not dom.get("is_included", True):
                continue
            seq_start = dom.get("iali") or dom.get("ienv") or 0
            seq_end   = dom.get("jali") or dom.get("jenv") or 0
            dom_evalue = dom.get("ievalue") or evalue

            domains.append({
                "name":        name,
                "description": desc,
                "e_value":     dom_evalue,
                "seq_start":   seq_start,
                "seq_end":     seq_end,
            })

    print(f"[HMMER] complete — {len(domains)} domain hits")
    return {"status": "ok", "data": {"domains": domains}, "error": ""}


def _err(msg: str) -> dict:
    return {"status": "error", "data": {"domains": []}, "error": msg}
