"""
Pipeline orchestrator.

submit_job()  — called by the Flask route; starts a background thread and
                returns a job_id immediately so the user can be redirected.

_run()        — runs in the background thread; calls each module in sequence,
                writes results to the job folder, updates status.json at every step.

Module execution order:
  1. retrieval      — fetch/validate sequence (MUST succeed; others are skipped if it fails)
  2. properties     — instant local calculation
  3. blast          — SLOW (~2-10 min), runs concurrently with all annotation tools
  4. hmmer          — EBI HMMER/Pfam API
  5. phobius        — EBI Phobius API (signal peptide + TM topology)
  6. cdd            — NCBI CDD Batch Web CD-Search
  7. scanprosite    — ExPASy ScanProsite (patterns, profiles, sites)
  8. smart          — SMART web adapter (HTML, best-effort)
  9. interproscan   — EBI InterProScan REST (integrative, slowest)
  10. annotation_merger — merges all annotation sources into unified schema

NOTE: Figure generation is NOT automatic. The user selects annotations on the
results page and clicks "Generate Figure" to invoke the mydomains module
on-demand via a separate Flask route.

HMMER, Phobius, CDD, ScanProsite, SMART, and InterProScan run in parallel
sub-threads. BLAST runs concurrently but is joined AFTER the annotation merge
so a slow BLAST run does not delay the results page.
"""

import os
import json
import time
import threading
from datetime import datetime, timezone

from pipeline.utils import jobs as job_store
from pipeline.modules import (
    retrieval, properties, blast, hmmer, phobius, signalp,
    cdd, scanprosite, smart, interproscan, coils,
    annotation_merger,
)

try:
    from pipeline.modules.motif_search import run_motif_analysis as _run_motifs
    _HAS_MOTIF = True
except ImportError:
    _HAS_MOTIF = False

# Keep old figures module import for backward compatibility (not used in new pipeline)
try:
    from pipeline.modules import figures as _figures_legacy
    _HAS_LEGACY_FIGURES = True
except ImportError:
    _HAS_LEGACY_FIGURES = False


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
        # Step 2b — motif search (instant, local — runs only if queries given)
        # ----------------------------------------------------------------
        motif_queries = input_data.get("motif_queries") or []
        motif_results_list = []
        if motif_queries and _HAS_MOTIF:
            try:
                motif_out = _run_motifs(seq, motif_queries)
                motif_results_list = motif_out.get("data", {}).get("motif_results", [])
                job_store.write_result(job_id, "motifs.json", motif_out)
                total_hits = sum(m.get("hit_count", 0) for m in motif_results_list)
                print(f"[runner] motif search: {len(motif_queries)} pattern(s), {total_hits} total hits")
            except Exception as _me:
                print(f"[runner] motif search error: {_me}")

        # ----------------------------------------------------------------
        # Steps 3–9 — annotation tools in parallel; BLAST runs independently
        #
        # BLAST does NOT feed annotation_merger — it goes straight into the
        # summary dict. So we start it in a separate thread and only join it
        # AFTER the merge step, so slow BLAST runs don't block annotations.
        # ----------------------------------------------------------------
        results = {
            "blast": None, "hmmer": None, "phobius": None, "signalp": None,
            "cdd": None, "scanprosite": None, "smart": None, "interproscan": None,
            "coils": None,
        }

        run_blast        = input_data.get("run_blast",        True)
        run_hmmer        = input_data.get("run_hmmer",        True)
        run_phobius      = input_data.get("run_phobius",      True)
        run_signalp      = input_data.get("run_signalp",      False)  # EBI endpoint unavailable
        run_cdd          = input_data.get("run_cdd",          True)
        run_scanprosite  = input_data.get("run_scanprosite",  True)
        run_smart        = input_data.get("run_smart",        True)
        run_interproscan = input_data.get("run_interproscan", False)   # very slow, opt-in only
        run_coils        = input_data.get("run_coils",        True)
        blast_max_hits   = input_data.get("blast_max_hits",   10)
        blast_database   = input_data.get("blast_database",   "swissprot")

        # BLAST thread — started immediately but joined after merge
        blast_thread = None
        if run_blast:
            job_store.set_module_status(job_id, "blast", "running")
            blast_thread = threading.Thread(
                target=_run_module,
                args=(blast.run, seq, job_dir, "blast", results, job_id),
                kwargs={"max_hits": blast_max_hits, "database": blast_database},
                daemon=True,
            )
            blast_thread.start()
        else:
            job_store.set_module_status(job_id, "blast", "skipped")

        # Annotation threads — joined before merge
        ann_threads = []

        def _add(module_fn, key, run_flag):
            if run_flag:
                job_store.set_module_status(job_id, key, "running")
                ann_threads.append(threading.Thread(
                    target=_run_module,
                    args=(module_fn, seq, job_dir, key, results, job_id),
                    daemon=True,
                ))
            else:
                job_store.set_module_status(job_id, key, "skipped")

        _add(hmmer.run,        "hmmer",        run_hmmer)
        _add(phobius.run,      "phobius",      run_phobius)
        _add(signalp.run,      "signalp",      run_signalp)
        _add(cdd.run,          "cdd",          run_cdd)
        _add(scanprosite.run,  "scanprosite",  run_scanprosite)
        _add(smart.run,        "smart",        run_smart)
        _add(interproscan.run, "interproscan", run_interproscan)
        _add(coils.run,        "coils",        run_coils)

        for t in ann_threads:
            t.start()

        # All annotation threads share a single 10-minute deadline.
        # args passed to _run_module: (module_fn, seq, job_dir, key, results, job_id)
        #                                          idx:  0       1     2        3  4       5
        # args tuple passed to _run_module: (module_fn, seq, job_dir, key, results, job_id)
        # so key is at index 3.
        _MODULE_TIMEOUT = 600   # 10 minutes total shared deadline
        _deadline = time.time() + _MODULE_TIMEOUT
        for t in ann_threads:
            remaining = max(1.0, _deadline - time.time())
            t.join(timeout=remaining)
            if t.is_alive():
                key = t._args[3] if hasattr(t, '_args') and len(t._args) > 3 else "unknown"
                print(f"[runner] module '{key}' timed out — continuing without it")
                job_store.set_module_status(job_id, key, "error")
                if isinstance(key, str) and key in results:
                    results[key] = {"status": "error", "data": {}, "error": f"Module timed out after {_MODULE_TIMEOUT // 60} min"}

        # ----------------------------------------------------------------
        # Step 10 — merge all annotation sources
        # ----------------------------------------------------------------
        job_store.set_module_status(job_id, "annotations", "running")
        merged = annotation_merger.merge(
            hmmer_result        = results.get("hmmer"),
            phobius_result      = results.get("phobius"),
            signalp_result      = results.get("signalp"),
            cdd_result          = results.get("cdd"),
            scanprosite_result  = results.get("scanprosite"),
            smart_result        = results.get("smart"),
            interproscan_result = results.get("interproscan"),
            coils_result        = results.get("coils"),
        )
        job_store.write_result(job_id, "annotations.json", merged)
        job_store.set_module_status(job_id, "annotations", "complete")

        # ----------------------------------------------------------------
        # Wait for BLAST — give it up to 10 more minutes after annotations finish.
        # If it's still alive after that, mark error and continue writing summary.
        # ----------------------------------------------------------------
        if blast_thread is not None:
            blast_thread.join(timeout=600)
            if blast_thread.is_alive():
                print(f"[runner] BLAST timed out — writing summary without BLAST results")
                results["blast"] = {
                    "status": "error",
                    "data":   {},
                    "error":  "BLAST timed out after 10 minutes. Try a smaller sequence or a faster database (SwissProt).",
                }
                job_store.set_module_status(job_id, "blast", "error")

        # ----------------------------------------------------------------
        # Write combined summary
        # ----------------------------------------------------------------
        annotation_list  = merged.get("annotations", [])
        low_conf_list    = merged.get("low_confidence_annotations", [])
        raw_ann_list     = merged.get("raw_annotations", [])
        ann_debug        = merged.get("debug", {})
        domain_count = sum(1 for a in annotation_list
                           if a.get("feature_type") in ("domain", "family", "repeat"))
        site_count   = sum(1 for a in annotation_list
                           if a.get("feature_type") in
                           ("active_site", "binding_site", "metal_binding", "disulfide", "site"))
        tm_count     = sum(1 for a in annotation_list if a.get("feature_type") == "transmembrane")
        has_sp       = any(a.get("feature_type") == "signal_peptide" for a in annotation_list)

        blast_result_data = (results.get("blast") or {}).get("data", {})
        blast_result_data["max_hits_requested"] = blast_max_hits

        # Carry through resolved protein candidates list (from mRNA/gene resolution)
        if ret.get("resolved_proteins"):
            pass  # already in ret dict, passed through below

        summary = {
            "job_id":        job_id,
            "completed_at":  datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
            "retrieval":     ret,
            "properties":    prop_result.get("data", {}),
            "blast":         blast_result_data,
            "hmmer":         (results.get("hmmer")   or {}).get("data", {}),
            "phobius":       (results.get("phobius") or {}).get("data", {}),
            "signalp":       (results.get("signalp") or {}).get("data", {}),
            "motif_results":           motif_results_list,
            "annotations":            annotation_list,
            "low_confidence_annotations": low_conf_list,
            "raw_annotations":        raw_ann_list,
            "annotation_summary": {
                "domain_count":   domain_count,
                "site_count":     site_count,
                "tm_helix_count": tm_count,
                "has_signal_peptide": has_sp,
                "source_summary": merged.get("source_summary", {}),
            },
            "annotation_debug": ann_debug,
            "figure_ok":     False,
            "figure_renderer": "none",
            "module_errors": _collect_errors(results, prop_result, merged),
        }
        job_store.write_result(job_id, "summary.json", summary)
        job_store.set_status(job_id, "complete")
        print(f"[runner] job {job_id} complete — {len(annotation_list)} annotations")

    except Exception as e:
        import traceback
        msg = f"Unexpected pipeline error: {e}\n{traceback.format_exc()}"
        print(f"[runner] ERROR in job {job_id}: {msg}")
        job_store.set_status(job_id, "error", error=msg)


def _run_module(module_fn, seq, job_dir, key, results_dict, job_id, **kwargs):
    """Thread target: run one module function, store result."""
    try:
        result = module_fn(seq, job_dir, **kwargs)
        results_dict[key] = result
        job_store.write_result(job_id, f"{key}.json", result)
        # parse_warning is a soft failure — don't mark as error
        if result.get("status") in ("ok", "parse_warning"):
            status = "complete"
        elif result.get("status") == "endpoint_unverified":
            status = "error"
        else:
            status = "error"
        job_store.set_module_status(job_id, key, status)
    except Exception as e:
        import traceback
        results_dict[key] = {"status": "error", "data": {}, "error": str(e)}
        job_store.set_module_status(job_id, key, "error")
        print(f"[runner] module {key} raised: {e}\n{traceback.format_exc()}")


def _collect_errors(results: dict, prop_result: dict, merged: dict) -> dict:
    errors = {}
    for key, val in results.items():
        if val and val.get("status") not in ("ok", "parse_warning", None):
            errors[key] = val.get("error", "unknown error")
    if prop_result.get("status") != "ok":
        errors["properties"] = prop_result.get("error", "")
    return errors
