"""
ScanProsite adapter — searches a sequence against PROSITE patterns and profiles.

Endpoint:
  POST https://prosite.expasy.org/cgi-bin/prosite/PSScan.cgi
  Fields:
    meta=opt1, meta1_protein=opt1, seq=<sequence>,
    skip=on, lowscore=1, output=nice, submit=START THE SCAN

Response: HTML page (ScanView.cgi) with profile match table + predicted features.

Two-pass HTML parse:
  Pass 1 (regex): Find PSImage.cgi URLs to extract profile/pattern match coordinates
                  and accession numbers. Also searches surrounding context for score.
  Pass 2 (table): Walk <tr> rows looking for feature-level annotations:
                    DISULFID  → disulfide
                    ACT_SITE  → active_site
                    METAL     → metal_binding
                    BINDING   → metal_binding (if Zn/Fe/Cu etc.) or binding_site
                    MOD_RES   → site
                  Skip rows where condition cell contains "not true" or "incomplete group"
                  (those are features predicted absent from this protein).

Feature types emitted:
  domain, motif, active_site, binding_site, metal_binding, disulfide, site
"""

import re
import os
import time
import requests

_SCAN_URL  = "https://prosite.expasy.org/cgi-bin/prosite/PSScan.cgi"
_VIEW_BASE = "https://prosite.expasy.org/cgi-bin/prosite/ScanView.cgi"
_TIMEOUT   = 120

_UA = "Mozilla/5.0 (compatible; ProtPipe/1.0; +https://singhlab.net)"

_METAL_LIGANDS = {
    "zn", "zinc", "fe", "iron", "cu", "copper", "co", "cobalt",
    "mn", "manganese", "ni", "nickel", "mo", "molybdenum",
    "mg", "magnesium", "ca", "calcium",
}

# Known accessions → explicit feature type
_SITE_ACCESSIONS: dict[str, str] = {
    "PS00107": "active_site",    # Serine proteases, trypsin active site
    "PS00213": "active_site",    # Cysteine proteases active site
    "PS00142": "active_site",    # Alcohol dehydrogenase zinc-binding
    "PS00062": "active_site",    # Alkaline phosphatase active site
    "PS00154": "active_site",    # Aspartyl proteases active site
    "PS00634": "active_site",    # Zinc metallopeptidase active site
    "PS00139": "active_site",    # Serine/threonine kinase active site
    "PS00108": "active_site",    # Protein kinase active site
    "PS00028": "metal_binding",  # Zinc metallopeptidase zinc ligands
    "PS00214": "metal_binding",  # Zinc finger, C2H2
    "PS50157": "metal_binding",  # Zinc finger, C2H2 profile
    "PS01237": "metal_binding",  # Zinc finger, RING-type
    "PS51858": "disulfide",
    "PS51859": "disulfide",
    "PS00103": "binding_site",   # ATP/GTP binding (P-loop)
    "PS00017": "binding_site",   # ATP/GTP-binding site (A)
}


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run(sequence: str, job_dir: str) -> dict:
    """
    Search sequence against PROSITE patterns + profiles, including predicted
    feature-level annotations (DISULFID, ACT_SITE, BINDING, METAL).

    Returns:
        {
          "status": "ok" | "parse_warning" | "error",
          "data": { "annotations": [...] },
          "error": str
        }
    """
    print(f"[ScanProsite] submitting sequence ({len(sequence)} residues)…")

    try:
        r = requests.post(
            _SCAN_URL,
            data={
                "meta":          "opt1",
                "meta1_protein": "opt1",
                "seq":           sequence,
                "skip":          "on",
                "lowscore":      "1",
                "output":        "nice",
                "submit":        "START THE SCAN",
            },
            headers={"User-Agent": _UA, "Accept": "text/html,text/plain"},
            timeout=_TIMEOUT,
            allow_redirects=True,
        )
        r.raise_for_status()
        html = r.text
    except requests.RequestException as e:
        return _err(f"ScanProsite network error: {e}")

    # Handle JS redirect to ScanView.cgi (asynchronous result case)
    scanfile_m = re.search(r"ScanView\.cgi\?scanfile=([^\s'\"&<>]+)", html)
    if scanfile_m:
        fetched = _fetch_scanview(scanfile_m.group(1))
        if fetched is None:
            return _err("ScanProsite: timed out waiting for ScanView results")
        html = fetched

    debug_path = os.path.join(job_dir, "scanprosite_debug.html")
    _save_debug(html, debug_path)

    annotations = _parse_nice_html(html)

    if annotations is None:
        lower = html.lower()
        if "no match" in lower or "no hit" in lower or not html.strip():
            return {"status": "ok", "data": {"annotations": []}, "error": ""}
        print("[ScanProsite] WARNING: could not parse HTML result. Saved debug.")
        return {
            "status": "parse_warning",
            "data":   {"annotations": []},
            "error":  "ScanProsite HTML parse failed — format may have changed",
        }

    print(f"[ScanProsite] complete — {len(annotations)} annotations (profiles + features)")
    return {"status": "ok", "data": {"annotations": annotations}, "error": ""}


# ---------------------------------------------------------------------------
# ScanView polling
# ---------------------------------------------------------------------------

def _fetch_scanview(scanfile: str, max_wait: int = 120, interval: int = 5) -> str | None:
    """Poll ScanView.cgi until results are ready. Returns HTML or None on timeout."""
    url = f"{_VIEW_BASE}?scanfile={scanfile}"
    print(f"[ScanProsite] following ScanView scanfile: {scanfile}")
    deadline = time.time() + max_wait
    while time.time() < deadline:
        try:
            rv = requests.get(
                url,
                headers={"User-Agent": _UA},
                timeout=30,
                allow_redirects=True,
            )
            rv.raise_for_status()
            html = rv.text
            # Still-processing pages typically say "please wait" / "refresh"
            if re.search(r"please.{0,5}wait|still.{0,10}processing|refresh", html, re.IGNORECASE):
                time.sleep(interval)
                continue
            return html
        except requests.RequestException as e:
            print(f"[ScanProsite] ScanView poll error: {e}")
            time.sleep(interval)
    return None


# ---------------------------------------------------------------------------
# HTML parsing — two passes
# ---------------------------------------------------------------------------

def _parse_nice_html(html: str) -> list | None:
    """
    Parse ScanView.cgi HTML.

    Pass 1 — regex scan for PSImage.cgi URLs → profile match coordinates + accession.
              Also extracts score from surrounding context (~500 chars around the URL).
    Pass 2 — table row scan for feature-level annotations (DISULFID, ACT_SITE, etc.)
              skipping absent features (condition "not true" / "incomplete group").

    Returns [] for genuine no-hits pages, None if the page structure can't be parsed.
    """
    # Remove HTML comments (can contain PSImage-like text that confuses the parser)
    html_clean = re.sub(r"<!--.*?-->", "", html, flags=re.DOTALL)

    # -----------------------------------------------------------------------
    # Pass 1: profile matches from PSImage.cgi URLs
    # -----------------------------------------------------------------------
    profile_anns = []
    seen_profiles: set = set()

    psimage_re = re.compile(
        r"PSImage\.cgi\?hit=(\d+),(\d+),(PS\d{5}),([^&\"'\s<>]+)"
    )
    for m in psimage_re.finditer(html_clean):
        start = int(m.group(1))
        end   = int(m.group(2))
        acc   = m.group(3)
        label = m.group(4).replace("+", " ").replace("%20", " ").strip()

        key = (acc, start, end)
        if key in seen_profiles:
            continue
        seen_profiles.add(key)

        # Look for score in surrounding context
        pos     = m.start()
        ctx     = html_clean[max(0, pos - 600): pos + 600]
        score   = None
        score_m = re.search(r"score\s*[=:]\s*([\d.]+)", ctx, re.IGNORECASE)
        if score_m:
            try:
                score = float(score_m.group(1))
            except ValueError:
                pass

        feature_type = _classify_accession(acc)
        evidence     = acc if score is None else f"{acc}  Score {score:.3f}"

        profile_anns.append({
            "source":       "ScanProsite",
            "feature_type": feature_type,
            "start":        start,
            "end":          end,
            "label":        label or acc,
            "accession":    acc,
            "description":  "",
            "e_value":      None,
            "score":        score,
            "evidence":     evidence,
        })

    # Decide if this is a parseable page
    if not profile_anns:
        if re.search(r"PS\d{5}", html_clean):
            # Has PROSITE accessions but we couldn't find PSImage URLs — format changed
            return None
        # No accessions at all → treat as "no hits" (caller will check no-hit text)
        return None

    # -----------------------------------------------------------------------
    # Pass 2: feature rows from HTML tables
    # -----------------------------------------------------------------------
    feature_anns = _parse_feature_table(html_clean)

    return profile_anns + feature_anns


def _parse_feature_table(html: str) -> list:
    """
    Walk all <tr> rows in the HTML looking for predicted feature annotations.

    Feature row first cell is one of: ACT_SITE, BINDING, DISULFID, DISULFIDE,
    METAL, MOD_RES, LIPID, CARBOHYD, CROSSLNK, SITE.

    Rows where condition cell contains "not true" or "incomplete group" are skipped
    — those are features predicted to be ABSENT in this protein.
    """
    _FEATURE_KEYWORDS = frozenset({
        "ACT_SITE", "BINDING", "DISULFID", "DISULFIDE",
        "METAL", "MOD_RES", "LIPID", "CARBOHYD", "CROSSLNK", "SITE",
    })

    features = []

    row_re  = re.compile(r"<tr[^>]*>(.*?)</tr>",    re.DOTALL | re.IGNORECASE)
    cell_re = re.compile(r"<t[dh][^>]*>(.*?)</t[dh]>", re.DOTALL | re.IGNORECASE)
    tag_re  = re.compile(r"<[^>]+>")

    for row_m in row_re.finditer(html):
        row_html = row_m.group(1)
        cells = [
            tag_re.sub("", c.group(1)).strip()
            for c in cell_re.finditer(row_html)
        ]
        if not cells:
            continue

        ftype_raw = cells[0].upper().strip()
        if ftype_raw not in _FEATURE_KEYWORDS:
            continue

        # Skip absent features
        all_text = " ".join(cells).lower()
        if "not true" in all_text or "incomplete group" in all_text:
            continue

        # Parse start position from cell[1]
        if len(cells) < 2:
            continue
        start_m = re.search(r"\d+", cells[1])
        if not start_m:
            continue
        start = int(start_m.group())

        # Parse end position from cell[2] (sites may have end = start)
        end = start
        if len(cells) >= 3 and cells[2]:
            end_m = re.search(r"\d+", cells[2])
            if end_m:
                end = int(end_m.group())
        if end < start:
            end = start

        # Collect note / ligand text from subsequent cells
        note = ""
        for c in cells[3:]:
            if c and (c.startswith("/") or "ligand" in c.lower() or "note" in c.lower()):
                note = c
                break
        if not note and len(cells) > 3:
            # Take whatever is in the 4th cell as the note
            note = cells[3]

        feature_type, description = _classify_feature(ftype_raw, note)
        evidence_parts = [ftype_raw]
        if note:
            evidence_parts.append(note)

        features.append({
            "source":       "ScanProsite",
            "feature_type": feature_type,
            "start":        start,
            "end":          end,
            "label":        ftype_raw,
            "accession":    "",
            "description":  description,
            "e_value":      None,
            "score":        None,
            "evidence":     "  ".join(evidence_parts),
        })

    return features


# ---------------------------------------------------------------------------
# Classification helpers
# ---------------------------------------------------------------------------

def _classify_accession(accession: str) -> str:
    """Map a PS accession to a feature_type string."""
    if accession in _SITE_ACCESSIONS:
        return _SITE_ACCESSIONS[accession]
    # PS5xxxx = profiles → domain;  PS0xxxx = patterns → motif
    if re.match(r"PS5\d{4}", accession):
        return "domain"
    if re.match(r"PS0\d{4}", accession):
        return "motif"
    return "motif"


def _classify_feature(ftype_raw: str, note: str) -> tuple[str, str]:
    """Map raw feature keyword + note text to (feature_type, description)."""
    note_lower = note.lower()

    if ftype_raw in ("DISULFID", "DISULFIDE"):
        return "disulfide", "Disulfide bond"
    if ftype_raw == "ACT_SITE":
        return "active_site", note or "Active site"
    if ftype_raw == "METAL":
        return "metal_binding", note or "Metal binding"
    if ftype_raw == "BINDING":
        for metal in _METAL_LIGANDS:
            if metal in note_lower:
                return "metal_binding", note or "Metal binding"
        return "binding_site", note or "Binding site"
    # MOD_RES, LIPID, CARBOHYD, CROSSLNK, SITE
    return "site", note or ftype_raw


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def _save_debug(text: str, path: str) -> None:
    try:
        with open(path, "w", encoding="utf-8", errors="replace") as f:
            f.write(text)
    except Exception:
        pass


def _err(msg: str) -> dict:
    print(f"[ScanProsite] ERROR: {msg}")
    return {"status": "error", "data": {"annotations": []}, "error": msg}
