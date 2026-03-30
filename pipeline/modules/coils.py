"""
Coiled-coil prediction via LUPAS (NPSA-PRABI, Lyon).

API: POST https://npsa-prabi.ibcp.fr/cgi-bin/primanal_lupas.pl
Parameters:
    notice:    protein sequence (raw single-letter codes)
    matrix:    2   (MTIDK matrix — default, most widely used)
    weight:    N   (no weighting)
    format:    p   (per-residue probability output)
    ali_width: 70
Response: HTML; predictions are in a <pre> block.
          Lines: "   <pos>   <AA>   <probability>"

Annotations returned: contiguous regions with probability >= 0.5,
merged if gap <= 5 residues, minimum 14 residues (one heptad repeat × 2).
"""

import os
import re
import requests

_URL       = "https://npsa-prabi.ibcp.fr/cgi-bin/primanal_lupas.pl"
_TIMEOUT   = 60
_THRESH    = 0.5   # per-residue probability threshold
_MERGE_GAP = 5     # merge segments separated by this many residues or fewer
_MIN_LEN   = 14    # discard coiled-coil segments shorter than this


def run(sequence: str, job_dir: str) -> dict:
    """
    POST sequence to LUPAS, parse per-residue probabilities, return coiled-coil annotations.

    Returns:
        {
            "status":  "ok" | "error",
            "data":    {
                "annotations":   [ {unified annotation dict}, ... ],
                "probabilities": [ float, ... ]   # per-residue, 0-indexed
            },
            "error":   str
        }
    """
    clean_seq = sequence.strip().upper()
    if not clean_seq:
        return _err("Empty sequence passed to coils module.")

    try:
        r = requests.post(
            _URL,
            data={
                "notice":    clean_seq,
                "matrix":    "2",
                "weight":    "N",
                "format":    "p",
                "ali_width": "70",
            },
            timeout=_TIMEOUT,
        )
        r.raise_for_status()
    except requests.RequestException as e:
        return _err(f"LUPAS request failed: {e}")

    # Server returns ISO-8859-1
    html = r.content.decode("iso-8859-1", errors="replace")

    # Save raw HTML for debugging
    try:
        with open(os.path.join(job_dir, "coils_debug.html"), "w", encoding="utf-8", errors="replace") as fh:
            fh.write(html)
    except OSError:
        pass

    probs = _parse_probabilities(html)
    if probs is None:
        return _err(
            "LUPAS response could not be parsed — the server may have changed its output format. "
            "Check coils_debug.html in the job directory."
        )

    regions = _find_regions(probs)
    annotations = []
    for start, end in regions:
        annotations.append({
            "source":           "LUPAS",
            "accession":        "COIL",
            "label":            "Coiled coil",
            "feature_type":     "coiled_coil",
            "start":            start,
            "end":              end,
            "score":            None,
            "evalue":           None,
            "description":      f"Predicted coiled-coil region (LUPAS MTIDK matrix, probability ≥ {_THRESH})",
            "source_support":   ["LUPAS"],
            "display_priority": 45,
        })

    print(f"[coils] LUPAS: {len(clean_seq)} aa → {len(regions)} coiled-coil region(s)")
    return {
        "status": "ok",
        "data": {
            "annotations":   annotations,
            "probabilities": probs,
        },
        "error": "",
    }


# ---------------------------------------------------------------------------
# Parse LUPAS HTML output
# ---------------------------------------------------------------------------

def _parse_probabilities(html: str) -> list | None:
    """
    Extract per-residue coiled-coil probabilities from LUPAS HTML response.

    Looks for <pre> blocks containing lines of the form:
        <pos>  <AA>  <probability>
    e.g.:
        1   M   0.000
        2   E   0.025

    Returns list[float] length == sequence_length (0-indexed), or None on failure.
    """
    pre_blocks = re.findall(r'<pre[^>]*>(.*?)</pre>', html, re.IGNORECASE | re.DOTALL)

    line_pat = re.compile(r'^\s*(\d+)\s+([A-Z])\s+([\d.]+)\s*$', re.MULTILINE)

    for block in pre_blocks:
        # Strip any residual HTML tags (e.g. <b>, <font>)
        text = re.sub(r'<[^>]+>', '', block)
        matches = line_pat.findall(text)
        if len(matches) < 5:
            continue

        max_pos = max(int(m[0]) for m in matches)
        # Build 1-indexed array then slice to 0-indexed
        raw = [0.0] * (max_pos + 1)
        for pos_str, _aa, prob_str in matches:
            pos = int(pos_str)
            try:
                raw[pos] = float(prob_str)
            except (ValueError, IndexError):
                pass
        return raw[1:]  # drop index 0 → result is 0-indexed

    return None


# ---------------------------------------------------------------------------
# Region detection
# ---------------------------------------------------------------------------

def _find_regions(probs: list) -> list:
    """
    Find contiguous high-probability coiled-coil regions.

    1. Mark positions with prob >= _THRESH.
    2. Merge adjacent segments separated by <= _MERGE_GAP residues.
    3. Discard regions shorter than _MIN_LEN residues.

    Returns list of (start, end) 1-indexed tuples.
    """
    # First pass: raw segments
    segments = []
    in_seg   = False
    seg_start = 0
    for i, p in enumerate(probs):
        pos = i + 1  # 1-indexed
        if p >= _THRESH and not in_seg:
            in_seg    = True
            seg_start = pos
        elif p < _THRESH and in_seg:
            segments.append((seg_start, pos - 1))
            in_seg = False
    if in_seg:
        segments.append((seg_start, len(probs)))

    if not segments:
        return []

    # Merge close segments
    merged = [segments[0]]
    for start, end in segments[1:]:
        if start - merged[-1][1] <= _MERGE_GAP:
            merged[-1] = (merged[-1][0], end)
        else:
            merged.append((start, end))

    # Discard too-short segments
    return [(s, e) for s, e in merged if e - s + 1 >= _MIN_LEN]


def _err(msg: str) -> dict:
    return {"status": "error", "data": {"annotations": [], "probabilities": []}, "error": msg}
