"""
Shared polling helper for EBI Job Dispatcher async jobs.

The EBI Job Dispatcher pattern is:
  1. POST to /run  →  returns plain-text job ID
  2. GET  /status/{jobId}  →  returns one of: RUNNING, FINISHED, ERROR, FAILURE, NOT_FOUND, QUEUED
  3. GET  /result/{jobId}/{resultType}  →  returns result content

⚠️  ENDPOINT VERIFICATION REQUIRED
These base URLs are based on documented EBI Job Dispatcher REST API endpoints.
Before relying on them in production, run the curl tests described in the
implementation plan to confirm live responses.

Dispatcher base: https://www.ebi.ac.uk/Tools/services/rest/
"""

import time
import requests

_DEFAULT_POLL_INTERVAL = 15   # seconds between status checks
_DEFAULT_TIMEOUT       = 2700  # 45 minutes total


def submit_and_wait(
    tool: str,
    params: dict,
    result_type: str = "out",
    poll_interval: int = _DEFAULT_POLL_INTERVAL,
    timeout: int = _DEFAULT_TIMEOUT,
    label: str = "",
) -> dict:
    """
    Submit a job to EBI Job Dispatcher, poll until FINISHED, return result.

    Args:
        tool:         Tool subdirectory name, e.g. "phobius", "clustalo", "iprscan5"
        params:       Form-encoded POST dict (must include 'email' and 'sequence')
        result_type:  Which result type to fetch when complete (default "out")
        poll_interval: Seconds between status polls
        timeout:      Give up after this many seconds
        label:        Human-readable label for log messages

    Returns:
        {
          "ok": True/False,
          "text": str,    # raw result text (if ok)
          "error": str    # error message (if not ok)
        }
    """
    base = f"https://www.ebi.ac.uk/Tools/services/rest/{tool}"

    # Step 1: submit
    try:
        r = requests.post(f"{base}/run", data=params, timeout=30)
        r.raise_for_status()
        job_id = r.text.strip()
        if not job_id:
            return {"ok": False, "text": "", "error": f"{label}: empty job ID returned by EBI"}
        print(f"[EBI/{tool}] job submitted: {job_id}")
    except requests.RequestException as e:
        return {"ok": False, "text": "", "error": f"{label}: submission failed: {e}"}

    # Step 2: poll
    start = time.time()
    while time.time() - start < timeout:
        time.sleep(poll_interval)
        try:
            sr = requests.get(f"{base}/status/{job_id}", timeout=15)
            sr.raise_for_status()
            status = sr.text.strip()
            print(f"[EBI/{tool}] {job_id} status: {status}")
        except requests.RequestException as e:
            print(f"[EBI/{tool}] status check failed: {e}")
            continue

        if status == "FINISHED":
            break
        if status in ("ERROR", "FAILURE", "NOT_FOUND"):
            return {"ok": False, "text": "", "error": f"{label}: EBI job ended with status {status}"}
        # RUNNING / QUEUED — keep waiting
    else:
        return {"ok": False, "text": "", "error": f"{label}: timed out after {timeout // 60} minutes"}

    # Step 3: fetch result
    try:
        rr = requests.get(f"{base}/result/{job_id}/{result_type}", timeout=30)
        rr.raise_for_status()
        return {"ok": True, "text": rr.text, "error": ""}
    except requests.RequestException as e:
        return {"ok": False, "text": "", "error": f"{label}: result fetch failed: {e}"}
