"""
Figure generator — PROSITE MyDomains proxy with internal SVG fallback.

Strategy (in order):
  1. Convert merged annotations to MyDomains data format.
  2. POST to https://prosite.expasy.org/cgi-bin/prosite/mydomains/
  3. If the response is image/png → save as domain_figure.png, status=ok.
  4. If the response is HTML containing an <img> tag → extract the PNG URL,
     fetch it, save as domain_figure.png, status=ok.
  5. If MyDomains is unavailable / returns an error → fall through to the
     internal svgwrite renderer (see _render_svg_fallback).

MyDomains data format:
  Domains: start,stop,shape,color,label     (one per line)
  Ranges:  start,stop,type                  (one per line)
  Sites:   position,type                    (one per line)

  Sections are separated by a blank line.

  Shapes: 1=rectangle, 2=rounded rect, 3=right-arrow, 4=left-arrow, 5=ellipse, 6=diamond
  Colors: 1=orange, 2=green, 3=blue, 4=red
  Range types: 0=dashed, 1=solid
  Site types:  0=small tick, 1=pin marker
"""

import os
import re
import requests

try:
    import svgwrite
    _SVGWRITE = True
except ImportError:
    _SVGWRITE = False

_MYDOMAINS_URL = "https://prosite.expasy.org/cgi-bin/prosite/mydomains/"
_HSCALE        = 0.8
_TIMEOUT       = 30

# MyDomains color cycle (1-4) assigned per unique accession
_COLOR_CYCLE = [3, 2, 1, 4]   # blue, green, orange, red

# MyDomains shape assignments by feature type
_SHAPE_MAP = {
    "domain":         2,   # rounded rectangle
    "family":         2,
    "repeat":         5,   # ellipse
    "motif":          1,   # rectangle
    "coiled_coil":    5,
    "low_complexity": 1,
    "active_site":    1,
    "binding_site":   1,
    "metal_binding":  1,
    "disulfide":      1,
    "site":           1,
    "signal_peptide": 3,   # right-arrow
    "transmembrane":  4,   # left-arrow (shows as block)
}

# Internal SVG constants (fallback)
_SVG_WIDTH  = 900
_MARGIN_L   = 60
_MARGIN_R   = 60
_BAR_Y      = 90
_BAR_H      = 36
_FONT_SANS  = "system-ui, -apple-system, Arial, sans-serif"
_FONT_MONO  = "ui-monospace, 'SF Mono', Consolas, monospace"

_PALETTE = [
    ("#2563eb", "#1e40af"),   # blue
    ("#059669", "#065f46"),   # green
    ("#d97706", "#92400e"),   # amber
    ("#7c3aed", "#4c1d95"),   # violet
    ("#db2777", "#9d174d"),   # pink
    ("#0891b2", "#0e7490"),   # cyan
    ("#65a30d", "#3f6212"),   # lime
    ("#dc2626", "#991b1b"),   # red
]
_SP_COLOR   = "#fef08a"; _SP_BORDER = "#a16207"
_TM_COLOR   = "#fca5a5"; _TM_BORDER = "#b91c1c"
_SITE_COLOR = "#f59e0b"; _SITE_BORDER = "#b45309"


def run(sequence: str, annotations: list, job_dir: str, extra_commands: str = "") -> dict:
    """
    Generate domain architecture figure.

    Args:
        sequence:    full query sequence (for length)
        annotations: list of unified annotation dicts from annotation_merger
        job_dir:     absolute path to job folder

    Returns:
        {
          "status":   "ok" | "error",
          "path":     str  (absolute path to domain_figure.png or .svg),
          "renderer": "mydomains" | "svgwrite" | "none",
          "error":    str
        }
    """
    seq_len = max(len(sequence), 1)

    # ------------------------------------------------------------------
    # Step 1: Try PROSITE MyDomains
    # ------------------------------------------------------------------
    mydomains_result = _try_mydomains(seq_len, annotations, job_dir, extra_commands)
    if mydomains_result["status"] == "ok":
        return mydomains_result

    # MyDomains failed — return the error directly (no SVG fallback).
    print(f"[mydomains] MyDomains unavailable: {mydomains_result['error']}")
    return mydomains_result


# ---------------------------------------------------------------------------
# PROSITE MyDomains integration
# ---------------------------------------------------------------------------

def _try_mydomains(seq_len: int, annotations: list, job_dir: str, extra_commands: str = "") -> dict:
    """
    Attempt to get the domain figure PNG from PROSITE MyDomains.

    Flow:
    1. POST annotation data → receive HTML page
    2. Extract PSImage.cgi URL from the HTML (this is the actual generated figure)
    3. Fetch the PSImage.cgi URL → image/png
    4. Optionally fetch RulerImage.cgi (scale bar) and compose a combined PNG

    Debug artifact: saves raw HTML to <job_dir>/mydomains_debug.html so the
    response can be inspected if image extraction fails.
    """
    data_str = _build_mydomains_data(annotations)
    # Append any user-supplied custom MyDomains commands
    if extra_commands and extra_commands.strip():
        data_str = data_str + "\n\n" + extra_commands.strip() if data_str else extra_commands.strip()

    _headers = {
        "User-Agent": "Mozilla/5.0 (compatible; ProtPipe/1.0)",
        "Referer":    "https://prosite.expasy.org/mydomains/",
    }

    try:
        r = requests.post(
            _MYDOMAINS_URL,
            data={
                "len":    str(seq_len),
                "hscale": str(_HSCALE),
                "data":   data_str,
            },
            headers=_headers,
            timeout=_TIMEOUT,
        )
    except requests.RequestException as e:
        return _err_r(f"MyDomains network error: {e}")

    if not r.ok:
        return _err_r(f"MyDomains HTTP {r.status_code}")

    content_type = r.headers.get("Content-Type", "")

    # Case A: direct image response (unlikely but handle it)
    if "image/" in content_type:
        png_path = os.path.join(job_dir, "domain_figure.png")
        with open(png_path, "wb") as f:
            f.write(r.content)
        print(f"[mydomains] MyDomains direct image saved to {png_path}")
        return {"status": "ok", "path": png_path, "renderer": "mydomains", "error": ""}

    # Case B: HTML page (the normal MyDomains response)
    # Save raw HTML for debugging
    html = r.text
    debug_html_path = os.path.join(job_dir, "mydomains_debug.html")
    try:
        with open(debug_html_path, "w", encoding="utf-8", errors="replace") as fh:
            fh.write(html)
        print(f"[mydomains] debug HTML saved to {debug_html_path}")
    except OSError:
        pass

    # Extract PSImage.cgi URL (the actual generated figure — NOT the logo or icons)
    psimage_url = _extract_img_url(html)
    ruler_url   = _extract_ruler_url(html)

    print(f"[mydomains] PSImage URL: {psimage_url or 'NOT FOUND'}")
    print(f"[mydomains] Ruler URL:   {ruler_url or 'NOT FOUND'}")

    if not psimage_url:
        return _err_r(
            "MyDomains returned HTML but PSImage.cgi URL not found. "
            f"Check {debug_html_path} for the raw response."
        )

    # Fetch the domain figure
    try:
        img_r = requests.get(psimage_url, headers=_headers, timeout=_TIMEOUT)
    except requests.RequestException as e:
        return _err_r(f"MyDomains PSImage fetch failed: {e}")

    if not img_r.ok or "image" not in img_r.headers.get("Content-Type", ""):
        return _err_r(
            f"MyDomains PSImage fetch returned {img_r.status_code} "
            f"Content-Type={img_r.headers.get('Content-Type','?')}"
        )

    figure_png = img_r.content

    # Optionally fetch ruler and stack vertically
    if ruler_url:
        try:
            ruler_r = requests.get(ruler_url, headers=_headers, timeout=_TIMEOUT)
            if ruler_r.ok and "image" in ruler_r.headers.get("Content-Type", ""):
                combined = _stack_images_vertically(figure_png, ruler_r.content)
                if combined:
                    figure_png = combined
                    print("[mydomains] figure + ruler stacked successfully")
        except Exception as e:
            print(f"[mydomains] ruler fetch/compose failed (non-fatal): {e}")

    png_path = os.path.join(job_dir, "domain_figure.png")
    with open(png_path, "wb") as f:
        f.write(figure_png)
    print(f"[mydomains] domain figure saved: {png_path} ({len(figure_png)} bytes)")
    return {"status": "ok", "path": png_path, "renderer": "mydomains", "error": ""}


def _stack_images_vertically(top_bytes: bytes, bottom_bytes: bytes) -> bytes | None:
    """
    Stack two PNG images vertically using only stdlib (no Pillow).
    Returns combined PNG bytes, or None if compositing fails.

    Uses Python's built-in zlib + struct to read/write PNG chunks.
    Falls back gracefully — caller ignores None and uses top image alone.
    """
    try:
        import zlib, struct, io

        def read_png(data):
            """Return (width, height, idat_chunks, bit_depth, color_type)."""
            if data[:8] != b'\x89PNG\r\n\x1a\n':
                return None
            pos = 8
            width = height = bit_depth = color_type = 0
            idat = b''
            while pos < len(data):
                length = struct.unpack('>I', data[pos:pos+4])[0]
                chunk_type = data[pos+4:pos+8]
                chunk_data = data[pos+8:pos+8+length]
                if chunk_type == b'IHDR':
                    width, height = struct.unpack('>II', chunk_data[:8])
                    bit_depth = chunk_data[8]
                    color_type = chunk_data[9]
                elif chunk_type == b'IDAT':
                    idat += chunk_data
                pos += 12 + length
            return width, height, idat, bit_depth, color_type

        top = read_png(top_bytes)
        bot = read_png(bottom_bytes)
        if not top or not bot:
            return None

        tw, th, t_idat, bd, ct = top
        bw, bh, b_idat, _, _  = bot

        # Must have same width, bit depth, color type to composite simply
        if tw != bw or bd != 8:
            return None

        # Decompress both IDAT streams
        t_raw = zlib.decompress(t_idat)
        b_raw = zlib.decompress(b_idat)

        # For RGB (ct=2) each row = 1 filter byte + width*3 bytes
        # For RGBA (ct=6) each row = 1 filter byte + width*4 bytes
        channels = {2: 3, 6: 4}.get(ct)
        if channels is None:
            return None

        row_bytes = 1 + tw * channels
        if len(t_raw) != th * row_bytes or len(b_raw) != bh * row_bytes:
            return None

        combined_raw = t_raw + b_raw
        combined_h   = th + bh
        compressed   = zlib.compress(combined_raw)

        # Build new PNG
        out = io.BytesIO()
        def write_chunk(chunk_type, data):
            out.write(struct.pack('>I', len(data)))
            out.write(chunk_type)
            out.write(data)
            crc = zlib.crc32(chunk_type + data) & 0xFFFFFFFF
            out.write(struct.pack('>I', crc))

        out.write(b'\x89PNG\r\n\x1a\n')
        ihdr = struct.pack('>IIBBBBB', tw, combined_h, bd, ct, 0, 0, 0)
        write_chunk(b'IHDR', ihdr)
        write_chunk(b'IDAT', compressed)
        write_chunk(b'IEND', b'')
        return out.getvalue()

    except Exception:
        return None


def _extract_img_url(html: str) -> str:
    """
    Extract the actual generated domain figure URL from MyDomains HTML response.

    The MyDomains page renders the figure via a CGI endpoint:
      /cgi-bin/prosite/PSImage.cgi?paramfile=/work/expasy/tmp/http/draw.XXXXXX&len=N&hscale=N

    This is NOT a static .png/.gif file — the path contains 'PSImage.cgi'.
    All other images on the page are navigation icons and the PROSITE logo.

    We must specifically target PSImage.cgi and ignore everything else.
    """
    # Primary target: PSImage.cgi (the actual generated domain figure)
    m = re.search(
        r'<img[^>]+src=["\']([^"\']*PSImage\.cgi[^"\']*)["\']',
        html, re.IGNORECASE
    )
    if m:
        url = m.group(1).replace("&amp;", "&")
        if url.startswith("/"):
            url = "https://prosite.expasy.org" + url
        return url

    return ""


def _extract_ruler_url(html: str) -> str:
    """Extract the ruler/scale-bar image URL (RulerImage.cgi) if present."""
    m = re.search(
        r'<img[^>]+src=["\']([^"\']*RulerImage\.cgi[^"\']*)["\']',
        html, re.IGNORECASE
    )
    if m:
        url = m.group(1).replace("&amp;", "&")
        if url.startswith("/"):
            url = "https://prosite.expasy.org" + url
        return url
    return ""


def _build_mydomains_data(annotations: list) -> str:
    """
    Convert unified annotations to MyDomains data format.

    Section 1 (domains): start,stop,shape,color,label
    Section 2 (ranges):  start,stop,type
    Section 3 (sites):   position,type

    Sections are separated by blank lines.

    Feature mapping:
      signal_peptide  → Domain: start,end,2,2,Signal peptide
      transmembrane   → Domain: start,end,4,1,TM helix
      disulfide       → Range (span): start,end,0  +  Site each cysteine: pos,1
      active_site     → Site: pos,1  (pin marker)
      metal_binding   → Site: pos,0  (tick)
      binding_site    → Site: pos,0
      site            → Site: pos,0
      domain/family/… → Domain with per-accession color cycling
    """
    domain_lines = []
    range_lines  = []
    site_lines   = []

    # Assign colors per unique accession (for structural domains only)
    color_map = {}
    color_idx = 0

    for ann in annotations:
        ftype  = ann.get("feature_type", "domain")
        start  = ann.get("start", 1)
        end    = ann.get("end", start)
        label  = ann.get("label", "")[:20]   # MyDomains truncates long labels
        acc    = ann.get("accession", "") or label

        if ftype == "signal_peptide":
            # Domain command: shape 2 (rounded rect), color 2 (green)
            domain_lines.append(f"{start},{end},2,2,Signal peptide")

        elif ftype == "transmembrane":
            # Domain command: shape 4 (left-arrow block), color 1 (orange)
            domain_lines.append(f"{start},{end},4,1,TM helix")

        elif ftype == "disulfide":
            # Range only — shows the bridge span; do NOT add site pins
            # (site pins look like active sites and are visually misleading for disulfides)
            # type 1 = solid line range
            range_lines.append(f"{start},{end},1")

        elif ftype == "active_site":
            # Pin marker at each catalytic residue
            site_lines.append(f"{start},1")
            if end != start:
                for pos in range(start + 1, end + 1):
                    site_lines.append(f"{pos},1")

        elif ftype in ("metal_binding", "binding_site", "site"):
            # Tick marker
            site_lines.append(f"{start},0")

        else:
            # Domain / family / motif / repeat / coiled_coil / low_complexity
            if acc not in color_map:
                color_map[acc] = _COLOR_CYCLE[color_idx % len(_COLOR_CYCLE)]
                color_idx += 1
            color  = color_map[acc]
            shape  = _SHAPE_MAP.get(ftype, 2)
            safe_label = label.replace(",", " ")
            domain_lines.append(f"{start},{end},{shape},{color},{safe_label}")

    parts = []
    if domain_lines:
        parts.append("\n".join(domain_lines))
    if range_lines:
        parts.append("\n".join(range_lines))
    if site_lines:
        parts.append("\n".join(site_lines))

    return "\n\n".join(parts)


def _err_r(msg: str) -> dict:
    return {"status": "error", "path": "", "renderer": "none", "error": msg}


# ---------------------------------------------------------------------------
# Internal SVG fallback (svgwrite)
# ---------------------------------------------------------------------------

def _render_svg_fallback(seq_len: int, annotations: list, job_dir: str) -> str:
    """Render domain architecture as SVG using svgwrite. Returns path to saved SVG."""
    out_path = os.path.join(job_dir, "domain_figure.svg")

    # Separate annotation types
    signal_anns = [a for a in annotations if a["feature_type"] == "signal_peptide"]
    tm_anns     = [a for a in annotations if a["feature_type"] == "transmembrane"]
    domain_anns = [a for a in annotations if a["feature_type"] in
                   ("domain", "family", "repeat", "motif", "coiled_coil", "low_complexity")]
    site_anns   = [a for a in annotations if a["feature_type"] in
                   ("active_site", "binding_site", "metal_binding", "disulfide", "site")]

    # Legend items
    legend_items = []
    color_map = {}
    color_idx = 0

    sp = signal_anns[0] if signal_anns else None
    if sp:
        legend_items.append((sp["label"], _SP_COLOR, _SP_BORDER))

    for ann in tm_anns[:1]:  # just one legend entry for TM
        legend_items.append((f"Transmembrane helix (×{len(tm_anns)})", _TM_COLOR, _TM_BORDER))

    for ann in domain_anns:
        acc = ann["accession"] or ann["label"]
        if acc not in color_map:
            ci = color_idx % len(_PALETTE)
            color_map[acc] = _PALETTE[ci]
            color_idx += 1
            legend_items.append((ann["label"][:60], _PALETTE[ci - 1 + 1 if ci + 1 < len(_PALETTE) else ci][0],
                                   _PALETTE[ci][1]))

    # Recompute color_map cleanly
    color_map = {}
    color_idx = 0
    for ann in domain_anns:
        acc = ann["accession"] or ann["label"]
        if acc not in color_map:
            ci = color_idx % len(_PALETTE)
            color_map[acc] = _PALETTE[ci]
            color_idx += 1
    legend_items = []
    if sp:
        legend_items.append((sp["label"], _SP_COLOR, _SP_BORDER))
    if tm_anns:
        legend_items.append((f"Transmembrane helix (×{len(tm_anns)})", _TM_COLOR, _TM_BORDER))
    seen_labels = set()
    for ann in domain_anns:
        acc = ann["accession"] or ann["label"]
        if acc not in seen_labels:
            seen_labels.add(acc)
            c, b = color_map.get(acc, (_PALETTE[0][0], _PALETTE[0][1]))
            legend_items.append((ann["label"][:60], c, b))

    leg_row_h  = 22
    leg_header = 28 if legend_items else 0
    leg_body_h = len(legend_items) * leg_row_h + 12
    no_data_h  = 36 if not legend_items else 0
    bar_section_h = _BAR_Y + _BAR_H + 32
    svg_height = bar_section_h + leg_header + leg_body_h + no_data_h + 20

    dwg = svgwrite.Drawing(out_path, size=(f"{_SVG_WIDTH}px", f"{svg_height}px"), profile="full")
    dwg.add(dwg.rect(insert=(0, 0), size=(_SVG_WIDTH, svg_height), fill="#ffffff"))
    dwg.add(dwg.rect(insert=(1, 1), size=(_SVG_WIDTH - 2, svg_height - 2),
                     fill="none", stroke="#e2e8f0", stroke_width=1, rx=8, ry=8))

    # Title bar
    dwg.add(dwg.rect(insert=(0, 0), size=(_SVG_WIDTH, 36), fill="#f8fafc", rx=8, ry=8))
    dwg.add(dwg.rect(insert=(0, 20), size=(_SVG_WIDTH, 16), fill="#f8fafc"))
    dwg.add(dwg.line(start=(0, 36), end=(_SVG_WIDTH, 36), stroke="#e2e8f0", stroke_width=1))
    dwg.add(dwg.text("Domain Architecture", insert=(20, 24),
                     font_size="13px", font_family=_FONT_SANS, font_weight="700", fill="#0f172a"))
    dwg.add(dwg.text(f"{seq_len} aa", insert=(_SVG_WIDTH - 20, 24),
                     text_anchor="end", font_size="11px", font_family=_FONT_MONO, fill="#94a3b8"))

    bar_x = _MARGIN_L
    bar_w = _SVG_WIDTH - _MARGIN_L - _MARGIN_R

    def px(pos):
        return bar_x + (pos - 1) / seq_len * bar_w

    def fw(s, e):
        return max(5.0, (e - s + 1) / seq_len * bar_w)

    # Backbone
    dwg.add(dwg.rect(insert=(bar_x, _BAR_Y), size=(bar_w, _BAR_H),
                     rx=5, ry=5, fill="#e2e8f0", stroke="#cbd5e1", stroke_width=1.5))

    # Signal peptide
    if sp:
        sw = fw(sp["start"], sp["end"])
        dwg.add(dwg.rect(insert=(bar_x, _BAR_Y), size=(sw, _BAR_H),
                         rx=5, ry=5, fill=_SP_COLOR, stroke=_SP_BORDER, stroke_width=1.5))
        mid = bar_x + sw / 2
        dwg.add(dwg.text("SP", insert=(mid, _BAR_Y - 7), text_anchor="middle",
                         font_size="10px", font_family=_FONT_SANS, font_weight="700", fill=_SP_BORDER))
        cx = bar_x + sw
        dwg.add(dwg.line(start=(cx, _BAR_Y - 4), end=(cx, _BAR_Y + _BAR_H + 4),
                         stroke=_SP_BORDER, stroke_width=1.5, **{"stroke-dasharray": "3,2"}))

    # TM helices
    for helix in tm_anns:
        hx = px(helix["start"])
        hw = fw(helix["start"], helix["end"])
        dwg.add(dwg.rect(insert=(hx, _BAR_Y), size=(hw, _BAR_H),
                         rx=3, ry=3, fill=_TM_COLOR, stroke=_TM_BORDER, stroke_width=1.5))
        mid = hx + hw / 2
        dwg.add(dwg.text("TM", insert=(mid, _BAR_Y + _BAR_H / 2 + 4), text_anchor="middle",
                         font_size="9px", font_family=_FONT_SANS, font_weight="700", fill=_TM_BORDER))

    # Domains
    inset = 4
    for ann in domain_anns:
        acc = ann["accession"] or ann["label"]
        color, border = color_map.get(acc, ("#6366f1", "#4338ca"))
        dx = px(ann["start"])
        dw = fw(ann["start"], ann["end"])
        dwg.add(dwg.rect(insert=(dx, _BAR_Y + inset), size=(dw, _BAR_H - inset * 2),
                         rx=3, ry=3, fill=color, stroke=border, stroke_width=1.5, opacity="0.92"))
        short = ann["label"].split(".")[0][:12]
        mid = dx + dw / 2
        if dw >= 32:
            dwg.add(dwg.text(short, insert=(mid, _BAR_Y + _BAR_H / 2 + 4), text_anchor="middle",
                             font_size="9px", font_family=_FONT_SANS, font_weight="700", fill="#ffffff"))
        elif dw >= 10:
            dwg.add(dwg.text(short, insert=(mid, _BAR_Y - 6), text_anchor="middle",
                             font_size="8px", font_family=_FONT_SANS, font_weight="700", fill=border))

    # Sites (pin markers)
    for ann in site_anns:
        sx = px(ann["start"])
        pin_color = "#ef4444" if "active" in ann["feature_type"] else _SITE_COLOR
        pin_border = "#991b1b" if "active" in ann["feature_type"] else _SITE_BORDER
        dwg.add(dwg.line(start=(sx, _BAR_Y - 12), end=(sx, _BAR_Y + _BAR_H + 4),
                         stroke=pin_border, stroke_width=1.5))
        dwg.add(dwg.circle(center=(sx, _BAR_Y - 14), r=4,
                           fill=pin_color, stroke=pin_border, stroke_width=1))

    # Ticks
    tick_top    = _BAR_Y + _BAR_H
    tick_bottom = tick_top + 7
    for tick_pos in _tick_positions(seq_len):
        tx = px(tick_pos)
        dwg.add(dwg.line(start=(tx, tick_top), end=(tx, tick_bottom), stroke="#94a3b8", stroke_width=1))
        dwg.add(dwg.text(str(tick_pos), insert=(tx, tick_bottom + 12), text_anchor="middle",
                         font_size="9px", font_family=_FONT_MONO, fill="#64748b"))
    dwg.add(dwg.text("1", insert=(bar_x, tick_bottom + 12), text_anchor="middle",
                     font_size="9px", font_family=_FONT_MONO, fill="#64748b"))

    # Legend
    leg_y_start = bar_section_h
    dwg.add(dwg.line(start=(0, leg_y_start), end=(_SVG_WIDTH, leg_y_start), stroke="#e2e8f0", stroke_width=1))
    if legend_items:
        dwg.add(dwg.text("Legend", insert=(20, leg_y_start + 20),
                         font_size="11px", font_family=_FONT_SANS, font_weight="700", fill="#475569"))
        item_y = leg_y_start + leg_header
        for label, fill, stroke in legend_items:
            dwg.add(dwg.rect(insert=(20, item_y - 11), size=(18, 13), rx=3, ry=3,
                             fill=fill, stroke=stroke, stroke_width=1))
            dwg.add(dwg.text(label[:80], insert=(44, item_y),
                             font_size="11px", font_family=_FONT_SANS, fill="#1e293b"))
            item_y += leg_row_h
    else:
        dwg.add(dwg.text("No annotated features found — run annotation services for domain data.",
                         insert=(_SVG_WIDTH // 2, leg_y_start + 22), text_anchor="middle",
                         font_size="11px", font_family=_FONT_SANS, fill="#94a3b8"))

    dwg.save()
    print(f"[mydomains] SVG fallback written to {out_path}")
    return out_path


def _tick_positions(seq_len: int) -> list:
    if seq_len <= 50:      step = 10
    elif seq_len <= 200:   step = 50
    elif seq_len <= 500:   step = 100
    elif seq_len <= 1000:  step = 200
    else:                  step = 500
    positions = list(range(step, seq_len, step))
    if seq_len not in positions:
        positions.append(seq_len)
    return positions
