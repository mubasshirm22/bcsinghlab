"""
Pipeline orchestrator.

submit_job()  — called by the Flask route; starts a background thread and
                returns a job_id immediately so the user can be redirected.

_run()        — runs in the background thread; calls each module in sequence,
                writes results to the job folder, updates status.json at every step.

Module execution order:
  1. retrieval   — fetch/validate sequence (MUST succeed; others are skipped if it fails)
  2. properties  — instant local calculation
  3. blast       — SLOW (~2-10 min), runs concurrently with hmmer + phobius
  4. hmmer       — EBI HMMER API (flagged: endpoint needs curl verification)
  5. phobius     — EBI Phobius API (flagged: endpoint needs curl verification)
  6. figures     — generates SVG from hmmer + phobius results

BLAST, HMMER, and Phobius run in parallel sub-threads to minimise total wait time.
"""

import os
import json
import threading
from datetime import datetime, timezone

from pipeline.utils import jobs as job_store
from pipeline.modules import retrieval, properties, blast, hmmer, phobius, figures


def submit_job(input_data: dict) -> str:
    """
    Create a job directory, spawn a background thread, and return the job_id.
    Never raises — errors are written to status.json.
    """
    job_id = job_store.new_job_id()
    job_dir = job_store.create_job_dir(job_id)

    t = threading.Thread(target=_run, args=(input_data, job_id, job_dir), daemon=True)
    t.name = f"protpipe-{job_id}"
    t.start()
    return job_id


def get_status(job_id: str) -> dict:
    """Return status dict for the given job."""
    return job_store.get_status(job_id)


def get_summary(job_id: str) -> dict:
    """Return the combined summary dict once the job is complete."""
    return job_store.get_summary(job_id)


# ---------------------------------------------------------------------------
# Background worker
# ---------------------------------------------------------------------------

def _run(input_data: dict, job_id: str, job_dir: str):
    try:
        job_store.set_status(job_id, "running")

        # ----------------------------------------------------------------
        # Step 1 — sequence retrieval (serial, must succeed first)
        # ----------------------------------------------------------------
        job_store.set_module_status(job_id, "retrieval", "running")
        ret = retrieval.run(input_data, job_dir)

        if ret["status"] != "ok":
            job_store.write_result(job_id, "retrieval.json", ret)
            job_store.set_module_status(job_id, "retrieval", "error")
            job_store.set_status(job_id, "error", error=ret["error"])
            return

        job_store.write_result(job_id, "retrieval.json", ret)
        job_store.set_module_status(job_id, "retrieval", "complete")
        seq = ret["sequence"]

        # ----------------------------------------------------------------
        # Step 2 — basic properties (instant, serial)
        # ----------------------------------------------------------------
        job_store.set_module_status(job_id, "properties", "running")
        prop_result = properties.run(seq, job_dir)
        job_store.write_result(job_id, "properties.json", prop_result)
        job_store.set_module_status(
            job_id, "properties",
            "complete" if prop_result["status"] == "ok" else "error"
        )

        # ----------------------------------------------------------------
        # Steps 3, 4, 5 — BLAST / HMMER / Phobius in parallel sub-threads
        # ----------------------------------------------------------------
        results = {"blast": None, "hmmer": None, "phobius": None}
        run_blast   = input_data.get("run_blast",   True)
        run_hmmer   = input_data.get("run_hmmer",   True)
        run_phobius = input_data.get("run_phobius", True)

        sub_threads = []

        if run_blast:
            job_store.set_module_status(job_id, "blast", "running")
            sub_threads.append(threading.Thread(
                target=_run_module,
                args=(blast.run, seq, job_dir, "blast", results, job_id),
                daemon=True,
            ))
        else:
            job_store.set_module_status(job_id, "blast", "skipped")

        if run_hmmer:
            job_store.set_module_status(job_id, "hmmer", "running")
            sub_threads.append(threading.Thread(
                target=_run_module,
                args=(hmmer.run, seq, job_dir, "hmmer", results, job_id),
                daemon=True,
            ))
        else:
            job_store.set_module_status(job_id, "hmmer", "skipped")

        if run_phobius:
            job_store.set_module_status(job_id, "phobius", "running")
            sub_threads.append(threading.Thread(
                target=_run_module,
                args=(phobius.run, seq, job_dir, "phobius", results, job_id),
                daemon=True,
            ))
        else:
            job_store.set_module_status(job_id, "phobius", "skipped")

        for t in sub_threads:
            t.start()
        for t in sub_threads:
            t.join()   # wait for all three before generating the figure

        # ----------------------------------------------------------------
        # Step 6 — SVG figure (needs hmmer + phobius results)
        # ----------------------------------------------------------------
        job_store.set_module_status(job_id, "figures", "running")
        domain_list  = (results.get("hmmer") or {}).get("data", {}).get("domains", [])
        phobius_data = (results.get("phobius") or {}).get("data", {})
        fig_result   = figures.run(seq, domain_list, phobius_data, job_dir)
        job_store.write_result(job_id, "figures.json", fig_result)
        job_store.set_module_status(
            job_id, "figures",
            "complete" if fig_result["status"] == "ok" else "error"
        )

        # ----------------------------------------------------------------
        # Write combined summary
        # ----------------------------------------------------------------
        summary = {
            "job_id":       job_id,
            "completed_at": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
            "retrieval":    ret,
            "properties":   prop_result.get("data", {}),
            "blast":        (results.get("blast")   or {}).get("data", {}),
            "hmmer":        (results.get("hmmer")   or {}).get("data", {}),
            "phobius":      (results.get("phobius") or {}).get("data", {}),
            "figure_ok":    fig_result.get("status") == "ok",
            "module_errors": _collect_errors(results, prop_result, fig_result),
        }
        job_store.write_result(job_id, "summary.json", summary)
        job_store.set_status(job_id, "complete")
        print(f"[runner] job {job_id} complete")

    except Exception as e:
        import traceback
        msg = f"Unexpected pipeline error: {e}\n{traceback.format_exc()}"
        print(f"[runner] ERROR in job {job_id}: {msg}")
        job_store.set_status(job_id, "error", error=msg)


def _run_module(module_fn, seq, job_dir, key, results_dict, job_id):
    """Thread target: run one module function, store result."""
    try:
        result = module_fn(seq, job_dir)
        results_dict[key] = result
        job_store.write_result(job_id, f"{key}.json", result)
        status = "complete" if result["status"] in ("ok",) else "error"
        # endpoint_unverified counts as partial — show warning, not hard error
        if result["status"] == "endpoint_unverified":
            status = "error"
        job_store.set_module_status(job_id, key, status)
    except Exception as e:
        import traceback
        results_dict[key] = {"status": "error", "data": {}, "error": str(e)}
        job_store.set_module_status(job_id, key, "error")
        print(f"[runner] module {key} raised: {e}\n{traceback.format_exc()}")


def _collect_errors(results: dict, prop_result: dict, fig_result: dict) -> dict:
    errors = {}
    for key, val in results.items():
        if val and val.get("status") not in ("ok", None):
            errors[key] = val.get("error", "unknown error")
    if prop_result.get("status") != "ok":
        errors["properties"] = prop_result.get("error", "")
    if fig_result.get("status") != "ok":
        errors["figures"] = fig_result.get("error", "")
    return errors
