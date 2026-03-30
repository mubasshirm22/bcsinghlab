"""
Motif search engine for protein sequences.

Public API
----------
parse_prosite(text)          → (segments, regex_str)
parse_advanced(text)         → (segments, regex_str)
segments_to_regex(segments)  → str
search(sequence, pattern)    → result dict
run_motif_analysis(sequence, motif_queries)  → pipeline result dict
validate_sequence(raw)       → (cleaned_str, error_str)

Segment schema
--------------
  {"type": "residue", "value": "H"}           single amino acid
  {"type": "residue", "value": "[ST]"}         residue group
  {"type": "gap", "min": 0, "max": 0}          adjacent (no intervening residues)
  {"type": "gap", "min": 2, "max": 2}          exactly 2 any residues
  {"type": "gap", "min": 20, "max": 80}        between 20 and 80 any residues
  {"type": "gap", "min": 0, "max": -1}         zero or more  (max -1 = unlimited)
  {"type": "gap", "min": 1, "max": -1}         one or more

All positions in search results use 1-based residue numbering.

Supported input syntax (parse_advanced)
----------------------------------------
  PROSITE-like (contains '-' or 'x('):
      H-E-x(2)-H
      H-E-x(2)-H-x(20,80)-E
      N-x-[ST]
      [KR]-x(2)-[DE]

  Regex-like (no '-'):
      HExxH           dots and x's are any residue
      HE..H
      HE.{2}H
      HE.{20,80}E
      [ST]x{2}[DE]
"""

import re
from typing import Optional

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: All standard amino acid one-letter codes
AMINO_ACIDS: frozenset = frozenset("ACDEFGHIKLMNPQRSTVWY")

#: Sentinel for "no upper bound" on a gap
_UNLIMITED = -1

#: Regex character class that matches any standard amino acid
_ANY_AA = "[ACDEFGHIKLMNPQRSTVWY]"

#: Preset motifs exposed to the UI (id, display name, pattern, description)
PRESETS: list[dict] = [
    {
        "id":      "hexxh",
        "name":    "HExxH",
        "pattern": "H-E-x(2)-H",
        "desc":    "Zinc-binding metalloprotease motif. The two histidines and glutamate coordinate a catalytic zinc ion.",
    },
    {
        "id":      "hexxh_e",
        "name":    "HExxH…E (neprilysin-type)",
        "pattern": "H-E-x(2)-H-x(20,80)-E",
        "desc":    "Extended metalloprotease motif: HExxH plus a downstream third zinc ligand (E) within 20–80 residues.",
    },
    {
        "id":      "cxxc",
        "name":    "CxxC",
        "pattern": "C-x(2)-C",
        "desc":    "Two cysteines separated by exactly two residues — common in zinc fingers, disulfide pairs, and redox motifs.",
    },
    {
        "id":      "nxst",
        "name":    "N-x-S/T (N-glycosylation)",
        "pattern": "N-x(1)-[ST]",
        "desc":    "Core N-linked glycosylation sequon: Asn followed by any residue then Ser or Thr (but not Pro at position 2 in biology; this search does not exclude Pro).",
    },
    {
        "id":      "rxsxe",
        "name":    "R…S…E",
        "pattern": "R-x(0,200)-S-x(0,200)-E",
        "desc":    "Flexible three-residue pattern: Arg, then Ser, then Glu at variable spacing. Useful for searching custom catalytic triads.",
    },
    {
        "id":      "kdel",
        "name":    "KDEL / HDEL (ER retention)",
        "pattern": "[KH]-D-E-L",
        "desc":    "C-terminal ER retention signal. Matches both KDEL and HDEL variants.",
    },
]

# ---------------------------------------------------------------------------
# Segment builder helpers (used by visual builder on frontend and by parsers)
# ---------------------------------------------------------------------------

def _residue(value: str) -> dict:
    return {"type": "residue", "value": value.upper()}

def _gap(mn: int, mx: int) -> dict:
    return {"type": "gap", "min": mn, "max": mx}


# ---------------------------------------------------------------------------
# segments_to_regex
# ---------------------------------------------------------------------------

def segments_to_regex(segments: list) -> str:
    """
    Convert a segments list to a Python regex pattern string.

    Gap with max=-1 (unlimited) uses open-ended quantifier {n,}.
    Adjacent residues (gap min=0, max=0) produce no gap characters.
    """
    parts = []
    for seg in segments:
        if seg["type"] == "residue":
            val = seg["value"].upper()
            if val in ("X", "."):
                parts.append(_ANY_AA)
            elif val.startswith("[") and val.endswith("]"):
                # Character class — keep as-is but uppercase contents
                parts.append(val)
            elif len(val) == 1:
                parts.append(re.escape(val))
            else:
                # Multi-char literal (e.g. typed directly) — escape each char
                parts.append("".join(re.escape(c) for c in val))

        elif seg["type"] == "gap":
            mn = int(seg.get("min", 0))
            mx = int(seg.get("max", 0))
            if mn == 0 and mx == 0:
                pass  # directly adjacent — no characters between them
            elif mx == _UNLIMITED:
                if mn == 0:
                    parts.append(f"{_ANY_AA}*")
                elif mn == 1:
                    parts.append(f"{_ANY_AA}+")
                else:
                    parts.append(f"{_ANY_AA}{{{mn},}}")
            elif mn == mx:
                parts.append(f"{_ANY_AA}{{{mn}}}")
            else:
                parts.append(f"{_ANY_AA}{{{mn},{mx}}}")

    return "".join(parts)


# ---------------------------------------------------------------------------
# PROSITE-like syntax parser
# ---------------------------------------------------------------------------

def parse_prosite(text: str) -> tuple[list, str]:
    """
    Parse PROSITE-like syntax into (segments, regex_str).

    Tokens are separated by '-'. Each token is one of:
      H         single amino acid
      x         single any-residue (1 any AA)
      x(n)      exactly n any-residues (gap)
      x(n,m)    n to m any-residues (gap)
      [ST]      residue group
      {P}       NOT P (converted to [^P] in regex — NOT implemented here,
                just passed through as a residue with exclusion notation)

    Examples:
      H-E-x(2)-H                 → HExxH
      H-E-x(2)-H-x(20,80)-E     → HExxH with downstream E
      N-x-[ST]                   → N-glycosylation sequon
      [KR]-x(2)-[DE]             → bipartite signal
    """
    text = text.strip()
    tokens = [t.strip() for t in text.split("-") if t.strip()]
    segments = []

    for tok in tokens:
        # x(n) or x(n,m) — gap with explicit count
        xm = re.match(r'^[xX]\((\d+)(?:,(\d+))?\)$', tok)
        if xm:
            mn = int(xm.group(1))
            mx = int(xm.group(2)) if xm.group(2) is not None else mn
            segments.append(_gap(mn, mx))
            continue

        # bare x or X — single any-residue
        if tok.lower() == "x":
            segments.append(_gap(1, 1))
            continue

        # [ST] style group
        if tok.startswith("[") and tok.endswith("]"):
            segments.append(_residue(tok.upper()))
            continue

        # {P} exclusion group — pass through as [^P]
        if tok.startswith("{") and tok.endswith("}"):
            inner = tok[1:-1].upper()
            segments.append({"type": "residue", "value": f"[^{inner}]"})
            continue

        # Single amino acid letter
        if len(tok) == 1 and tok.upper() in (AMINO_ACIDS | {"X"}):
            if tok.upper() == "X":
                segments.append(_gap(1, 1))
            else:
                segments.append(_residue(tok.upper()))
            continue

        # Multi-character token without dashes (e.g. "HExxH" as a single token)
        # Parse character by character
        for ch in tok:
            if ch.upper() in AMINO_ACIDS:
                segments.append(_residue(ch.upper()))
            elif ch in ("x", "X", "."):
                segments.append(_gap(1, 1))

    # Merge consecutive single-char gaps (x x → gap(2,2), etc.)
    segments = _merge_consecutive_gaps(segments)
    return segments, segments_to_regex(segments)


# ---------------------------------------------------------------------------
# Advanced / regex-like syntax parser
# ---------------------------------------------------------------------------

def parse_advanced(text: str) -> tuple[list, str]:
    """
    Parse simplified regex-like or PROSITE syntax.

    Routes to parse_prosite() if the text looks like PROSITE (contains '-'
    or 'x('), otherwise parses as regex-like:
      HExxH       → HE + 2 gaps + H
      HE..H       → same
      HE.{2}H     → same (explicit repeat)
      HE.{20,80}E → range repeat
      [ST]        → residue group

    Returns (segments, regex_str).
    """
    text = text.strip()
    if not text:
        return [], ""

    # Route to PROSITE parser if it looks PROSITE-like
    if "-" in text or re.search(r"[xX]\(", text):
        return parse_prosite(text)

    # --- Regex-like parse ---
    segments = []
    i = 0
    n = len(text)

    while i < n:
        ch = text[i]

        # .{n,m} or .{n} — explicit repeat
        if ch in (".", "x", "X") and i + 1 < n and text[i + 1] == "{":
            brace_end = text.find("}", i + 2)
            if brace_end != -1:
                inner = text[i + 2: brace_end]
                parts = [p.strip() for p in inner.split(",")]
                try:
                    mn = int(parts[0])
                    mx = int(parts[1]) if len(parts) > 1 else mn
                    segments.append(_gap(mn, mx))
                    i = brace_end + 1
                    continue
                except ValueError:
                    pass  # fall through to single char

        # x{n,m} — x with explicit repeat (no leading dot)
        if ch in ("x", "X") and i + 1 < n and text[i + 1] == "{":
            brace_end = text.find("}", i + 2)
            if brace_end != -1:
                inner = text[i + 2: brace_end]
                parts = [p.strip() for p in inner.split(",")]
                try:
                    mn = int(parts[0])
                    mx = int(parts[1]) if len(parts) > 1 else mn
                    segments.append(_gap(mn, mx))
                    i = brace_end + 1
                    continue
                except ValueError:
                    pass

        # Plain dot or x — any single residue
        if ch in (".", "x", "X"):
            segments.append(_gap(1, 1))
            i += 1
            continue

        # [ST] style character class
        if ch == "[":
            close = text.find("]", i + 1)
            if close != -1:
                segments.append(_residue(text[i: close + 1]))
                i = close + 1
                continue

        # {P} exclusion group
        if ch == "{":
            close = text.find("}", i + 1)
            if close != -1:
                inner = text[i + 1: close].upper()
                segments.append({"type": "residue", "value": f"[^{inner}]"})
                i = close + 1
                continue

        # Standard amino acid
        if ch.upper() in AMINO_ACIDS:
            segments.append(_residue(ch.upper()))
            i += 1
            continue

        # Unknown / whitespace — skip
        i += 1

    segments = _merge_consecutive_gaps(segments)
    return segments, segments_to_regex(segments)


# ---------------------------------------------------------------------------
# Gap merging helper
# ---------------------------------------------------------------------------

def _merge_consecutive_gaps(segments: list) -> list:
    """
    Merge adjacent gap segments into a single combined gap.
    E.g. [gap(1,1), gap(1,1)] → [gap(2,2)]
         [gap(1,1), gap(0,-1)] → [gap(1,-1)]
    """
    if not segments:
        return segments
    merged = []
    i = 0
    while i < len(segments):
        seg = segments[i]
        if seg["type"] == "gap":
            total_min = seg["min"]
            total_max = seg["max"]
            while i + 1 < len(segments) and segments[i + 1]["type"] == "gap":
                i += 1
                nxt = segments[i]
                total_min += nxt["min"]
                if total_max == _UNLIMITED or nxt["max"] == _UNLIMITED:
                    total_max = _UNLIMITED
                else:
                    total_max += nxt["max"]
            merged.append(_gap(total_min, total_max))
        else:
            merged.append(seg)
        i += 1
    return merged


# ---------------------------------------------------------------------------
# Sequence validation
# ---------------------------------------------------------------------------

def validate_sequence(raw: str) -> tuple[str, str]:
    """
    Clean and validate a protein sequence string.

    Accepts raw sequence or FASTA (strips header lines starting with '>').
    Returns (cleaned_uppercase_sequence, error_message).
    On success error_message is "".
    On failure cleaned_sequence is "" and error_message describes the problem.
    """
    lines = raw.strip().splitlines()
    # Strip all FASTA header lines
    seq_lines = [ln for ln in lines if not ln.startswith(">")]
    seq = "".join(seq_lines).upper()
    # Remove whitespace
    seq = re.sub(r"\s+", "", seq)

    if not seq:
        return "", "No sequence provided."

    # Allow standard AAs plus common ambiguity codes and gap/stop symbols
    # (we strip stop codon '*' and gap '-' before searching)
    valid = AMINO_ACIDS | frozenset("BZOUXJ*-")
    invalid_chars = sorted(set(seq) - valid)
    if invalid_chars:
        return "", (
            f"Sequence contains invalid characters: {', '.join(invalid_chars)}. "
            "Please provide a standard single-letter protein sequence."
        )

    # Strip non-AA symbols for actual searching
    seq_clean = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq)
    if not seq_clean:
        return "", "Sequence contained only non-amino-acid characters after cleaning."

    return seq_clean, ""


# ---------------------------------------------------------------------------
# Search engine
# ---------------------------------------------------------------------------

def search(sequence: str, pattern_input, overlapping: bool = True) -> dict:
    """
    Search a protein sequence for all occurrences of a motif.

    ``pattern_input`` may be:
      - str   : PROSITE-like or regex-like syntax string
      - list  : pre-built segments list
      - dict  : {"segments": [...]} or {"regex": "..."}

    Returns::

        {
            "status":     "ok" | "error",
            "hits":       [{"start": int, "end": int, "match": str}, ...],
            "hit_count":  int,
            "regex_used": str,
            "error":      str,
        }

    All positions are **1-based** inclusive (start=1 is the first residue).
    """
    regex_str = ""
    try:
        if isinstance(pattern_input, str):
            _, regex_str = parse_advanced(pattern_input)
        elif isinstance(pattern_input, list):
            regex_str = segments_to_regex(pattern_input)
        elif isinstance(pattern_input, dict):
            if "regex" in pattern_input:
                regex_str = pattern_input["regex"]
            elif "segments" in pattern_input:
                regex_str = segments_to_regex(pattern_input["segments"])
            else:
                return _err("Invalid pattern_input dict (need 'regex' or 'segments' key)")
        else:
            return _err(f"Invalid pattern_input type: {type(pattern_input)}")
    except Exception as exc:
        return _err(f"Pattern parse error: {exc}")

    if not regex_str:
        return _err("Empty pattern after parsing — nothing to search for.")

    # Compile
    try:
        compiled = re.compile(regex_str, re.IGNORECASE)
    except re.error as exc:
        return _err(f"Invalid regex '{regex_str}': {exc}")

    # Find matches.
    # For overlapping matches we wrap in a lookahead: (?=(PATTERN))
    # The lookahead consumes zero characters so the engine re-tries at pos+1.
    hits = []
    if overlapping:
        try:
            lookahead = re.compile(f"(?=({regex_str}))", re.IGNORECASE)
            for m in lookahead.finditer(sequence):
                matched = m.group(1)
                start   = m.start() + 1        # 1-based
                end     = m.start() + len(matched)  # 1-based inclusive
                hits.append({"start": start, "end": end, "match": matched})
        except re.error as exc:
            return _err(f"Overlapping search regex error: {exc}")
    else:
        for m in compiled.finditer(sequence):
            hits.append({
                "start": m.start() + 1,
                "end":   m.end(),
                "match": m.group(),
            })

    return {
        "status":     "ok",
        "hits":       hits,
        "hit_count":  len(hits),
        "regex_used": regex_str,
        "error":      "",
    }


# ---------------------------------------------------------------------------
# ProtPipe pipeline module entry point
# ---------------------------------------------------------------------------

def run_motif_analysis(sequence: str, motif_queries: list) -> dict:
    """
    ProtPipe pipeline module entry point.

    ``motif_queries`` is a list of dicts::

        [
            {
                "name":     str,             display label
                "pattern":  str,             PROSITE/advanced syntax (optional if segments given)
                "segments": list (optional)  from visual builder
            },
            ...
        ]

    Returns pipeline-compatible result dict::

        {
            "status": "ok",
            "data": {
                "motif_results": [
                    {
                        "name":      str,
                        "pattern":   str,
                        "regex":     str,
                        "hits":      [{"start": int, "end": int, "match": str}, ...],
                        "hit_count": int,
                        "error":     str,
                    },
                    ...
                ]
            },
            "error": ""
        }
    """
    results = []
    for q in (motif_queries or []):
        name     = q.get("name", "Custom motif")
        pattern  = q.get("pattern", "")
        segments = q.get("segments")

        if segments:
            result = search(sequence, segments)
        elif pattern:
            result = search(sequence, pattern)
        else:
            continue  # skip empty queries

        results.append({
            "name":      name,
            "pattern":   pattern,
            "regex":     result.get("regex_used", ""),
            "hits":      result.get("hits", []),
            "hit_count": result.get("hit_count", 0),
            "error":     result.get("error", ""),
        })

    return {
        "status": "ok",
        "data":   {"motif_results": results},
        "error":  "",
    }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _err(msg: str) -> dict:
    return {
        "status":     "error",
        "hits":       [],
        "hit_count":  0,
        "regex_used": "",
        "error":      msg,
    }


def segments_to_human(segments: list) -> str:
    """
    Convert a segments list to a human-readable description string.
    Example: H — directly adjacent — E — exactly 2 any — H
    """
    parts = []
    for seg in segments:
        if seg["type"] == "residue":
            v = seg["value"]
            if v.startswith("["):
                parts.append(f"any of {v}")
            else:
                parts.append(v)
        elif seg["type"] == "gap":
            mn, mx = seg["min"], seg["max"]
            if mn == 0 and mx == 0:
                parts.append("directly adjacent")
            elif mx == _UNLIMITED:
                if mn == 0:
                    parts.append("any number of residues")
                elif mn == 1:
                    parts.append("1 or more residues")
                else:
                    parts.append(f"{mn} or more residues")
            elif mn == mx:
                parts.append(f"exactly {mn} any residue{'s' if mn != 1 else ''}")
            else:
                parts.append(f"{mn}–{mx} any residues")
    return " — ".join(parts)


def segments_to_prosite(segments: list) -> str:
    """
    Render a segments list as a PROSITE-like pattern string.
    Example: H-E-x(2)-H-x(20,80)-E
    """
    parts = []
    for seg in segments:
        if seg["type"] == "residue":
            v = seg["value"]
            if v.startswith("[") or v.startswith("[^"):
                parts.append(v)
            else:
                parts.append(v)
        elif seg["type"] == "gap":
            mn, mx = seg["min"], seg["max"]
            if mn == 0 and mx == 0:
                continue  # adjacent → no token
            elif mx == _UNLIMITED:
                if mn == 0:
                    parts.append("x(0,*)")
                else:
                    parts.append(f"x({mn},*)")
            elif mn == mx:
                parts.append(f"x({mn})")
            else:
                parts.append(f"x({mn},{mx})")
    return "-".join(parts)
