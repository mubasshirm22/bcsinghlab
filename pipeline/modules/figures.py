"""
Domain architecture SVG figure generator.

Uses svgwrite (pip install svgwrite) — no external HTTP calls.
Draws a horizontal sequence bar with colored domain boxes,
signal peptide region, and TM helix markers.

Output: writes domain_figure.svg into the job directory.
"""

import os

try:
    import svgwrite
    _SVGWRITE = True
except ImportError:
    _SVGWRITE = False

# Color palette for domains — cycles if more than 12 domains
_DOMAIN_COLORS = [
    "#2563eb", "#059669", "#d97706", "#7c3aed", "#db2777",
    "#0891b2", "#65a30d", "#dc2626", "#9333ea", "#ea580c",
    "#0d9488", "#4f46e5",
]
_DOMAIN_BORDERS = [
    "#1e40af", "#065f46", "#92400e", "#4c1d95", "#9d174d",
    "#0e7490", "#3f6212", "#991b1b", "#581c87", "#c2410c",
    "#115e59", "#312e81",
]

_SP_COLOR  = "#fef08a"   # signal peptide — yellow
_SP_BORDER = "#a16207"
_TM_COLOR  = "#fca5a5"   # TM helices — rose/red
_TM_BORDER = "#b91c1c"

_SVG_WIDTH  = 820
_MARGIN_L   = 50
_MARGIN_R   = 50
_BAR_Y      = 90
_BAR_H      = 36          # backbone bar height
_FONT_SANS  = "system-ui, -apple-system, Arial, sans-serif"
_FONT_MONO  = "ui-monospace, 'SF Mono', Consolas, monospace"


def run(sequence: str, domains: list, phobius: dict, job_dir: str) -> dict:
    """
    Args:
        sequence: full query sequence (used for length scaling)
        domains:  list of dicts from hmmer.py [{name, description, seq_start, seq_end, e_value}]
        phobius:  dict from phobius.py data field (has_signal_peptide, tm_helices, etc.)
        job_dir:  absolute path to the job folder

    Returns:
        {"status": "ok"|"error", "path": str, "error": str}
    """
    if not _SVGWRITE:
        return {"status": "error", "path": "", "error": "svgwrite not installed (pip install svgwrite)"}

    seq_len = max(len(sequence), 1)
    out_path = os.path.join(job_dir, "domain_figure.svg")

    has_sp  = phobius.get("has_signal_peptide", False)
    sp_end  = phobius.get("signal_peptide_end", 0) or 0
    tm_list = phobius.get("tm_helices") or []
    n_dom   = len(domains)

    # --- Legend items ---
    legend_items = []
    if has_sp and sp_end > 0:
        legend_items.append(("Signal Peptide (residues 1–{})".format(sp_end), _SP_COLOR, _SP_BORDER))
    if tm_list:
        legend_items.append(("Transmembrane helix (×{})".format(len(tm_list)), _TM_COLOR, _TM_BORDER))

    color_map  = {}
    color_idx  = 0
    for dom in domains:
        name = dom.get("name", "?")
        if name not in color_map:
            ci = color_idx % len(_DOMAIN_COLORS)
            color_map[name] = (_DOMAIN_COLORS[ci], _DOMAIN_BORDERS[ci])
            color_idx += 1
            desc = dom.get("description", "")
            label = "{} — {}".format(name, desc) if desc else name
            legend_items.append((label, _DOMAIN_COLORS[ci], _DOMAIN_BORDERS[ci]))

    # --- SVG height ---
    leg_row_h  = 22
    leg_header = 28 if legend_items else 0
    leg_body_h = len(legend_items) * leg_row_h + 12
    no_data_h  = 36 if not legend_items else 0
    bar_section_h = _BAR_Y + _BAR_H + 32   # bar + ticks + small gap
    svg_height = bar_section_h + leg_header + leg_body_h + no_data_h + 20

    dwg = svgwrite.Drawing(out_path, size=(f"{_SVG_WIDTH}px", f"{svg_height}px"), profile="full")

    # Background
    dwg.add(dwg.rect(insert=(0, 0), size=(_SVG_WIDTH, svg_height), fill="#ffffff"))

    # Outer border
    dwg.add(dwg.rect(
        insert=(1, 1), size=(_SVG_WIDTH - 2, svg_height - 2),
        fill="none", stroke="#e2e8f0", stroke_width=1, rx=8, ry=8,
    ))

    # --- Title bar ---
    dwg.add(dwg.rect(insert=(0, 0), size=(_SVG_WIDTH, 36), fill="#f8fafc", rx=8, ry=8))
    dwg.add(dwg.rect(insert=(0, 20), size=(_SVG_WIDTH, 16), fill="#f8fafc"))  # square bottom corners
    dwg.add(dwg.line(start=(0, 36), end=(_SVG_WIDTH, 36), stroke="#e2e8f0", stroke_width=1))
    dwg.add(dwg.text(
        "Domain Architecture",
        insert=(20, 24),
        font_size="13px", font_family=_FONT_SANS, font_weight="700",
        fill="#0f172a",
    ))
    # Sequence length in title bar right side
    dwg.add(dwg.text(
        "{} aa".format(seq_len),
        insert=(_SVG_WIDTH - 20, 24),
        text_anchor="end",
        font_size="11px", font_family=_FONT_MONO,
        fill="#94a3b8",
    ))

    # --- Drawing area helpers ---
    bar_x = _MARGIN_L
    bar_w = _SVG_WIDTH - _MARGIN_L - _MARGIN_R

    def px(pos: int) -> float:
        """Convert 1-based sequence position → SVG x coordinate."""
        return bar_x + (pos - 1) / seq_len * bar_w

    def feature_w(start: int, end: int) -> float:
        return max(5.0, (end - start + 1) / seq_len * bar_w)

    # --- Backbone bar ---
    dwg.add(dwg.rect(
        insert=(bar_x, _BAR_Y),
        size=(bar_w, _BAR_H),
        rx=5, ry=5,
        fill="#e2e8f0",
        stroke="#cbd5e1",
        stroke_width=1.5,
    ))

    # --- Signal peptide ---
    if has_sp and sp_end > 0:
        sw = feature_w(1, sp_end)
        dwg.add(dwg.rect(
            insert=(bar_x, _BAR_Y),
            size=(sw, _BAR_H),
            rx=5, ry=5,
            fill=_SP_COLOR,
            stroke=_SP_BORDER,
            stroke_width=1.5,
        ))
        # Label above bar
        mid = bar_x + sw / 2
        dwg.add(dwg.text(
            "SP",
            insert=(mid, _BAR_Y - 7),
            text_anchor="middle",
            font_size="10px", font_family=_FONT_SANS, font_weight="700",
            fill=_SP_BORDER,
        ))
        # Cleavage site marker
        cx = bar_x + sw
        dwg.add(dwg.line(
            start=(cx, _BAR_Y - 4), end=(cx, _BAR_Y + _BAR_H + 4),
            stroke=_SP_BORDER, stroke_width=1.5,
            **{"stroke-dasharray": "3,2"},
        ))

    # --- TM helices ---
    for i, helix in enumerate(tm_list):
        hx = px(helix["start"])
        hw = feature_w(helix["start"], helix["end"])
        dwg.add(dwg.rect(
            insert=(hx, _BAR_Y),
            size=(hw, _BAR_H),
            rx=3, ry=3,
            fill=_TM_COLOR,
            stroke=_TM_BORDER,
            stroke_width=1.5,
        ))
        mid = hx + hw / 2
        dwg.add(dwg.text(
            "TM",
            insert=(mid, _BAR_Y + _BAR_H / 2 + 4),
            text_anchor="middle",
            font_size="9px", font_family=_FONT_SANS, font_weight="700",
            fill=_TM_BORDER,
        ))

    # --- Pfam / HMMER domains ---
    for dom in domains:
        name  = dom.get("name", "?")
        start = int(dom.get("seq_start", 1))
        end   = int(dom.get("seq_end", seq_len))
        color, border = color_map.get(name, ("#6366f1", "#4338ca"))

        dx = px(start)
        dw = feature_w(start, end)
        inset = 4  # domains sit inside the bar vertically

        dwg.add(dwg.rect(
            insert=(dx, _BAR_Y + inset),
            size=(dw, _BAR_H - inset * 2),
            rx=3, ry=3,
            fill=color,
            stroke=border,
            stroke_width=1.5,
            opacity="0.92",
        ))
        # Label inside box if wide enough, otherwise above
        short = name.split(".")[0]
        mid   = dx + dw / 2
        if dw >= 32:
            dwg.add(dwg.text(
                short,
                insert=(mid, _BAR_Y + _BAR_H / 2 + 4),
                text_anchor="middle",
                font_size="9px", font_family=_FONT_SANS, font_weight="700",
                fill="#ffffff",
            ))
        elif dw >= 10:
            # Just a colored block, put label above
            dwg.add(dwg.text(
                short,
                insert=(mid, _BAR_Y - 6),
                text_anchor="middle",
                font_size="8px", font_family=_FONT_SANS, font_weight="700",
                fill=border,
            ))

    # --- Position ticks ---
    tick_top    = _BAR_Y + _BAR_H
    tick_bottom = tick_top + 7
    for tick_pos in _tick_positions(seq_len):
        tx = px(tick_pos)
        dwg.add(dwg.line(
            start=(tx, tick_top), end=(tx, tick_bottom),
            stroke="#94a3b8", stroke_width=1,
        ))
        dwg.add(dwg.text(
            str(tick_pos),
            insert=(tx, tick_bottom + 12),
            text_anchor="middle",
            font_size="9px", font_family=_FONT_MONO,
            fill="#64748b",
        ))

    # Start label
    dwg.add(dwg.text(
        "1",
        insert=(bar_x, tick_bottom + 12),
        text_anchor="middle",
        font_size="9px", font_family=_FONT_MONO,
        fill="#64748b",
    ))

    # --- Legend section ---
    leg_y_start = bar_section_h

    if legend_items:
        dwg.add(dwg.line(
            start=(0, leg_y_start), end=(_SVG_WIDTH, leg_y_start),
            stroke="#e2e8f0", stroke_width=1,
        ))
        dwg.add(dwg.text(
            "Legend",
            insert=(20, leg_y_start + 20),
            font_size="11px", font_family=_FONT_SANS, font_weight="700",
            fill="#475569",
        ))
        item_y = leg_y_start + leg_header
        for label, fill, stroke in legend_items:
            dwg.add(dwg.rect(
                insert=(20, item_y - 11), size=(18, 13),
                rx=3, ry=3, fill=fill, stroke=stroke, stroke_width=1,
            ))
            dwg.add(dwg.text(
                label[:80],
                insert=(44, item_y),
                font_size="11px", font_family=_FONT_SANS,
                fill="#1e293b",
            ))
            item_y += leg_row_h
    else:
        # No features to show
        dwg.add(dwg.line(
            start=(0, leg_y_start), end=(_SVG_WIDTH, leg_y_start),
            stroke="#e2e8f0", stroke_width=1,
        ))
        dwg.add(dwg.text(
            "No annotated features found — run HMMER and Phobius for domain and topology data.",
            insert=(_SVG_WIDTH // 2, leg_y_start + 22),
            text_anchor="middle",
            font_size="11px", font_family=_FONT_SANS,
            fill="#94a3b8",
        ))

    dwg.save()
    print(f"[figures] SVG written to {out_path}")
    return {"status": "ok", "path": out_path, "error": ""}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _tick_positions(seq_len: int) -> list:
    """Return sensible tick mark positions given sequence length."""
    if seq_len <= 50:
        step = 10
    elif seq_len <= 200:
        step = 50
    elif seq_len <= 500:
        step = 100
    elif seq_len <= 1000:
        step = 200
    else:
        step = 500
    positions = list(range(step, seq_len, step))
    if seq_len not in positions:
        positions.append(seq_len)
    return positions
