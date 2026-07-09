import io
import re
from typing import Iterable

try:
    import svgwrite
except ImportError:
    svgwrite = None

try:
    # cairosvg raises OSError (not ImportError) when libcairo is missing, so
    # catch broadly — SVG rendering must keep working even if PNG/PDF can't.
    import cairosvg
except Exception:
    cairosvg = None


STYLE_PRESETS = {
    "journal": {
        "font_size": 11,
        "label_size": 11,
        "char_width": 8.0,
        "row_gap": 18,
        "line_padding": 26,
        "title_size": 16,
        "background": "#ffffff",
        "text": "#111827",
    },
    "presentation": {
        "font_size": 13,
        "label_size": 12,
        "char_width": 9.2,
        "row_gap": 22,
        "line_padding": 30,
        "title_size": 20,
        "background": "#fbfdff",
        "text": "#0f172a",
    },
    "bw": {
        "font_size": 11,
        "label_size": 11,
        "char_width": 8.0,
        "row_gap": 18,
        "line_padding": 24,
        "title_size": 16,
        "background": "#ffffff",
        "text": "#111111",
    },
}

PALETTES = {
    "colorblind": {
        "H": "#0072B2",
        "E": "#009E73",
        "C": "#D55E00",
        "sequence": "#111827",
        "support_low": "#DBEAFE",
        "support_high": "#1D4ED8",
        "region_fill": "#FEF3C7",
        "region_stroke": "#D97706",
        "mismatch": "#DC2626",
    },
    "grayscale": {
        "H": "#111111",
        "E": "#555555",
        "C": "#8A8A8A",
        "sequence": "#111111",
        "support_low": "#E5E7EB",
        "support_high": "#4B5563",
        "region_fill": "#F3F4F6",
        "region_stroke": "#6B7280",
        "mismatch": "#111111",
    },
    "muted": {
        "H": "#1D4ED8",
        "E": "#047857",
        "C": "#B45309",
        "sequence": "#0F172A",
        "support_low": "#E0F2FE",
        "support_high": "#0369A1",
        "region_fill": "#FDF2F8",
        "region_stroke": "#BE185D",
        "mismatch": "#B91C1C",
    },
}


# Distinct fills cycled through for domain-architecture boxes.
DOMAIN_COLORS = [
    "#2563EB", "#059669", "#D97706", "#DB2777",
    "#7C3AED", "#0891B2", "#CA8A04", "#DC2626",
]


def default_options():
    return {
        "row_length": 80,
        "style": "journal",
        "palette": "colorblind",
        "show_sequence": True,
        "show_consensus": True,
        "show_pdb": True,
        "show_confidence": True,
        "show_hydropathy": False,
        "hydropathy_window": 9,
        "predictors": [],
        "legend": True,
        "clean": False,
        "compare": False,
        "title": "SSPred Secondary Structure Summary",
        "regions": [],
        "domains": [],
    }


def render_svg(row: dict, options: dict | None = None) -> str:
    if svgwrite is None:
        raise RuntimeError("svgwrite is required for SSPred figure generation.")
    opts = default_options()
    if options:
        opts.update(options)

    style = STYLE_PRESETS.get(opts["style"], STYLE_PRESETS["journal"])
    palette = PALETTES.get(opts["palette"], PALETTES["colorblind"])
    sequence = (row.get("seq") or "").strip()
    if not sequence:
        raise ValueError("Sequence not available for figure generation.")

    predictors = _resolved_predictors(row, opts.get("predictors") or [])
    row_length = max(20, min(int(opts.get("row_length", 80)), 120))
    pdb_data = _parse_pdb(row.get("pdb"))
    consensus = row.get("majorityvote") or ""
    support = _consensus_support(sequence, predictors, consensus)
    regions = _normalize_regions(opts.get("regions") or [], len(sequence))
    domains = _normalize_domains(opts.get("domains") or [], len(sequence))

    tracks = []
    if opts.get("show_sequence", True):
        tracks.append(("Sequence", sequence, "sequence"))
    if opts.get("show_pdb", True) and pdb_data and pdb_data.get("secondary"):
        tracks.append((f"PDB {pdb_data.get('pdbid', '')}_{pdb_data.get('chain', '')}", pdb_data["secondary"], "structure"))
    for predictor in predictors:
        tracks.append((predictor["label"], predictor["pred"], "predictor"))
    if opts.get("show_consensus", True) and consensus:
        tracks.append(("Consensus", consensus, "consensus"))
    if opts.get("show_confidence", True) and support:
        tracks.append(("Support", support, "confidence"))
    if opts.get("show_hydropathy", False):
        hydro = _hydropathy_values(sequence, opts.get("hydropathy_window", 9))
        if hydro:
            tracks.append(("Hydropathy", hydro, "hydropathy"))

    char_width = style["char_width"]
    left_margin = 120
    right_margin = 32
    top_margin = 32 if opts.get("clean") else 72
    bottom_margin = 36
    block_height = _block_height(tracks, style, bool(regions), bool(domains))
    n_blocks = (len(sequence) + row_length - 1) // row_length
    width = int(left_margin + right_margin + row_length * char_width)
    height = int(top_margin + bottom_margin + n_blocks * block_height + (18 if opts.get("legend") else 0))

    drawing = svgwrite.Drawing(size=(width, height))
    drawing.viewbox(0, 0, width, height)
    drawing.add(drawing.rect(insert=(0, 0), size=(width, height), fill=style["background"]))

    if not opts.get("clean"):
        drawing.add(drawing.text(
            opts.get("title", "SSPred Secondary Structure Summary"),
            insert=(left_margin, 30),
            fill=style["text"],
            font_size=style["title_size"],
            font_family="Helvetica, Arial, sans-serif",
            font_weight="700",
        ))
        subtitle = f"{len(sequence)} aa · {len(predictors)} predictor track(s)"
        drawing.add(drawing.text(
            subtitle,
            insert=(left_margin, 50),
            fill="#64748B",
            font_size=11,
            font_family="Helvetica, Arial, sans-serif",
        ))

    for block_index in range(n_blocks):
        start = block_index * row_length
        end = min(len(sequence), start + row_length)
        base_y = top_margin + block_index * block_height
        _draw_block(
            drawing,
            row,
            tracks,
            palette,
            style,
            left_margin,
            base_y,
            char_width,
            start,
            end,
            regions,
            domains,
            bool(opts.get("compare")),
            pdb_data["secondary"] if pdb_data else "",
        )

    if opts.get("legend"):
        _draw_legend(drawing, tracks, palette, style, left_margin, height - 18, domains)

    return drawing.tostring()


def export_figure(row: dict, options: dict | None = None, output_format: str = "svg", scale: float = 3.0):
    svg_text = render_svg(row, options=options)
    if output_format == "svg":
        return io.BytesIO(svg_text.encode("utf-8")), "image/svg+xml"
    if cairosvg is None:
        raise RuntimeError("CairoSVG is required for PNG/PDF export.")
    if output_format == "png":
        data = cairosvg.svg2png(bytestring=svg_text.encode("utf-8"), scale=max(scale, 1.0))
        return io.BytesIO(data), "image/png"
    if output_format == "pdf":
        data = cairosvg.svg2pdf(bytestring=svg_text.encode("utf-8"))
        return io.BytesIO(data), "application/pdf"
    raise ValueError(f"Unsupported figure format: {output_format}")


def parse_region_text(text: str) -> list[dict]:
    regions = []
    for raw_line in (text or "").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        match = re.match(r"(\d+)\s*-\s*(\d+)(?::|\s+)(.+)$", line)
        if not match:
            continue
        start = int(match.group(1))
        end = int(match.group(2))
        label = match.group(3).strip()
        if end < start:
            start, end = end, start
        regions.append({"start": start, "end": end, "label": label})
    return regions


def parse_domain_text(text: str) -> list[dict]:
    """Parse domain-architecture input.

    Accepts lines of ``start-end: Label`` and optionally ``start-end: Label #hexcolor``.
    """
    domains = []
    for raw_line in (text or "").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        color = ""
        color_match = re.search(r"(#[0-9A-Fa-f]{6})\s*$", line)
        if color_match:
            color = color_match.group(1)
            line = line[:color_match.start()].strip()
        match = re.match(r"(\d+)\s*-\s*(\d+)(?::|\s+)(.+)$", line)
        if not match:
            continue
        start = int(match.group(1))
        end = int(match.group(2))
        label = match.group(3).strip()
        if end < start:
            start, end = end, start
        domains.append({"start": start, "end": end, "label": label, "color": color})
    return domains


def _resolved_predictors(row: dict, requested: Iterable[str]) -> list[dict]:
    available = []
    for key, value in row.items():
        if key.endswith("pred"):
            name = key[:-4]
            if name in {"pdb", "majorityvote"}:
                continue
            status = row.get(name + "stat")
            pred = row.get(name + "pred")
            if status in (1, 3) and pred:
                available.append({
                    "name": name,
                    "label": name.upper(),
                    "pred": pred,
                    "conf": row.get(name + "conf", ""),
                })
    if not requested:
        return available
    requested_set = {item.lower() for item in requested}
    return [item for item in available if item["name"].lower() in requested_set]


def _parse_pdb(value):
    if not value:
        return None
    if isinstance(value, dict):
        return value
    try:
        import json
        data = json.loads(value)
        if isinstance(data, dict):
            return data
    except Exception:
        return None
    return None


def _consensus_support(sequence: str, predictors: list[dict], consensus: str) -> list[float]:
    if not consensus:
        return []
    output = []
    for idx in range(len(sequence)):
        votes = {"H": 0, "E": 0, "C": 0}
        total = 0
        for predictor in predictors:
            if idx >= len(predictor["pred"]):
                continue
            state = predictor["pred"][idx]
            if state in votes:
                votes[state] += 1
                total += 1
        if total == 0:
            output.append(0.0)
            continue
        consensus_state = consensus[idx] if idx < len(consensus) else "X"
        if consensus_state in votes:
            output.append(votes[consensus_state] / total)
        else:
            output.append(max(votes.values()) / total)
    return output


def _normalize_regions(regions: list[dict], seq_len: int) -> list[dict]:
    output = []
    for region in regions:
        try:
            start = max(1, int(region.get("start", 0)))
            end = min(seq_len, int(region.get("end", 0)))
        except (AttributeError, TypeError, ValueError):
            continue
        if start > end:
            continue
        output.append({
            "start": start,
            "end": end,
            "label": str(region.get("label", "")).strip(),
        })
    return output


HYDROPATHY_TRACK_HEIGHT = 46
DOMAIN_BAND_HEIGHT = 20


def _track_height(track_kind: str, style: dict) -> int:
    if track_kind == "hydropathy":
        return HYDROPATHY_TRACK_HEIGHT
    return style["row_gap"]


def _block_height(tracks, style: dict, has_regions: bool, has_domains: bool = False) -> int:
    extra = 18 if has_regions else 0
    extra += DOMAIN_BAND_HEIGHT if has_domains else 0
    tracks_height = sum(_track_height(kind, style) for _, _, kind in tracks)
    return style["line_padding"] + tracks_height + extra + 12


def _draw_block(drawing, row, tracks, palette, style, left_margin, base_y, char_width, start, end, regions, domains, compare, pdb_secondary):
    font_size = style["font_size"]
    label_size = style["label_size"]
    block_width = (end - start) * char_width

    drawing.add(drawing.line(
        start=(left_margin, base_y + 8),
        end=(left_margin + block_width, base_y + 8),
        stroke="#CBD5E1",
        stroke_width=1,
    ))

    for position in range(start, end):
        residue_number = position + 1
        if residue_number % 10 == 0 or residue_number == 1:
            x = left_margin + (position - start) * char_width
            drawing.add(drawing.line(
                start=(x, base_y + 6),
                end=(x, base_y + 12),
                stroke="#94A3B8",
                stroke_width=1,
            ))
            drawing.add(drawing.text(
                str(residue_number),
                insert=(x, base_y),
                fill="#64748B",
                font_size=9,
                font_family="Helvetica, Arial, sans-serif",
            ))

    tracks_span = sum(_track_height(kind, style) for _, _, kind in tracks)
    domain_offset = 0
    if domains:
        _draw_domains(drawing, domains, left_margin, base_y + 14, char_width, start, end)
        domain_offset = DOMAIN_BAND_HEIGHT

    if regions:
        _draw_regions(drawing, regions, palette, left_margin, base_y + 14 + domain_offset, char_width, start, end, tracks_span)

    y = base_y + style["line_padding"] + domain_offset
    for label, track_data, track_kind in tracks:
        drawing.add(drawing.text(
            label,
            insert=(left_margin - 8, y),
            fill=style["text"],
            font_size=label_size,
            font_family="Helvetica, Arial, sans-serif",
            text_anchor="end",
            font_weight="600" if track_kind in {"consensus", "structure"} else "400",
        ))

        if track_kind == "hydropathy":
            _draw_hydropathy(drawing, track_data, palette, style, left_margin, y - 8, char_width, start, end)
            y += _track_height(track_kind, style)
            continue

        if track_kind == "confidence":
            _draw_support_strip(drawing, track_data, palette, left_margin, y - 10, char_width, start, end)
            y += _track_height(track_kind, style)
            continue

        segment = track_data[start:end]
        for offset, char in enumerate(segment):
            x = left_margin + offset * char_width
            if compare and track_kind == "consensus" and pdb_secondary:
                pdb_char = pdb_secondary[start + offset] if start + offset < len(pdb_secondary) else ""
                if pdb_char and pdb_char in {"H", "E", "C"} and char in {"H", "E", "C"} and pdb_char != char:
                    drawing.add(drawing.rect(
                        insert=(x - 0.5, y - font_size + 1),
                        size=(char_width, font_size + 3),
                        fill=palette["mismatch"],
                        fill_opacity=0.15,
                    ))
            drawing.add(drawing.text(
                char,
                insert=(x, y),
                fill=_char_color(char, track_kind, palette),
                font_size=font_size,
                font_family="'Courier New', monospace",
            ))
        y += _track_height(track_kind, style)


def _draw_regions(drawing, regions, palette, left_margin, top_y, char_width, start, end, track_height):
    drawn_labels = set()
    for region in regions:
        if region["end"] < start + 1 or region["start"] > end:
            continue
        region_start = max(region["start"] - 1, start)
        region_end = min(region["end"], end)
        x = left_margin + (region_start - start) * char_width
        width = max(char_width, (region_end - region_start) * char_width)
        drawing.add(drawing.rect(
            insert=(x, top_y),
            size=(width, track_height + 8),
            fill=palette["region_fill"],
            fill_opacity=0.22,
            stroke=palette["region_stroke"],
            stroke_width=0.8,
        ))
        label = region.get("label")
        if label and label not in drawn_labels:
            drawn_labels.add(label)
            drawing.add(drawing.text(
                label,
                insert=(x + 2, top_y - 2),
                fill=palette["region_stroke"],
                font_size=9,
                font_family="Helvetica, Arial, sans-serif",
                font_weight="600",
            ))


def _draw_support_strip(drawing, support_values, palette, left_margin, top_y, char_width, start, end):
    for offset, value in enumerate(support_values[start:end]):
        x = left_margin + offset * char_width
        bar_height = 12 * float(value)
        drawing.add(drawing.rect(
            insert=(x, top_y + (12 - bar_height)),
            size=(max(char_width - 1.0, 1.0), max(bar_height, 1.0)),
            fill=_mix_hex(palette["support_low"], palette["support_high"], float(value)),
            rx=0.8,
            ry=0.8,
        ))


def _draw_legend(drawing, tracks, palette, style, left_margin, y, domains=None):
    x = left_margin
    items = [
        ("Helix", palette["H"]),
        ("Strand", palette["E"]),
        ("Coil", palette["C"]),
    ]
    if any(track[2] == "confidence" for track in tracks):
        items.append(("Consensus support", palette["support_high"]))
    if any(track[2] == "hydropathy" for track in tracks):
        items.append(("Kyte-Doolittle hydropathy", "#0EA5E9"))
    for label, color in items:
        drawing.add(drawing.rect(insert=(x, y - 9), size=(12, 12), fill=color, rx=2, ry=2))
        drawing.add(drawing.text(
            label,
            insert=(x + 18, y + 1),
            fill=style["text"],
            font_size=10,
            font_family="Helvetica, Arial, sans-serif",
        ))
        x += 24 + len(label) * 6.0


def _hydropathy_values(sequence, window):
    """Kyte-Doolittle sliding-window hydropathy, one value per residue."""
    scale = {
        "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
        "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
        "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
        "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
    }
    residues = (sequence or "").upper()
    n = len(residues)
    if n == 0:
        return []
    try:
        window = int(window)
    except (TypeError, ValueError):
        window = 9
    window = max(1, min(window, n))
    if window % 2 == 0:
        window += 1
    half = window // 2
    values = []
    for center in range(n):
        lo = max(0, center - half)
        hi = min(n, center + half + 1)
        chunk = residues[lo:hi]
        values.append(sum(scale.get(aa, 0.0) for aa in chunk) / len(chunk))
    return values


def _draw_hydropathy(drawing, values, palette, style, left_margin, top_y, char_width, start, end):
    """Draw a filled Kyte-Doolittle line plot for the residues in [start, end)."""
    band = HYDROPATHY_TRACK_HEIGHT - 14
    scale_max = 4.5  # KD extremes are +/-4.5
    mid_y = top_y + band / 2.0
    block_width = (end - start) * char_width

    # Baseline (hydropathy = 0) and a faint frame.
    drawing.add(drawing.rect(
        insert=(left_margin, top_y), size=(max(block_width, 1.0), band),
        fill="#F8FAFC", stroke="#E2E8F0", stroke_width=0.6,
    ))
    drawing.add(drawing.line(
        start=(left_margin, mid_y), end=(left_margin + block_width, mid_y),
        stroke="#CBD5E1", stroke_width=0.8, stroke_dasharray="3,3",
    ))

    segment = values[start:end]
    if not segment:
        return

    def point(offset, value):
        x = left_margin + (offset + 0.5) * char_width
        clamped = max(-scale_max, min(scale_max, value))
        y = mid_y - (clamped / scale_max) * (band / 2.0)
        return x, y

    # Filled area to the midline for readability.
    area_points = [(left_margin + 0.5 * char_width, mid_y)]
    line_points = []
    for offset, value in enumerate(segment):
        px, py = point(offset, value)
        area_points.append((px, py))
        line_points.append((px, py))
    area_points.append((point(len(segment) - 1, segment[-1])[0], mid_y))
    drawing.add(drawing.polygon(points=area_points, fill="#0EA5E9", fill_opacity=0.16, stroke="none"))
    drawing.add(drawing.polyline(points=line_points, fill="none", stroke="#0284C7", stroke_width=1.4))


def _normalize_domains(domains, seq_len):
    output = []
    for idx, domain in enumerate(domains or []):
        try:
            d_start = max(1, int(domain.get("start", 0)))
            d_end = min(seq_len, int(domain.get("end", 0)))
        except (AttributeError, TypeError, ValueError):
            continue
        if d_start > d_end:
            continue
        color = str(domain.get("color", "")).strip() or DOMAIN_COLORS[idx % len(DOMAIN_COLORS)]
        output.append({
            "start": d_start,
            "end": d_end,
            "label": str(domain.get("label", "")).strip(),
            "color": color,
        })
    return output


def _draw_domains(drawing, domains, left_margin, top_y, char_width, start, end):
    """Draw domain-architecture boxes (Pfam/InterPro-style) across a block."""
    for domain in domains:
        if domain["end"] < start + 1 or domain["start"] > end:
            continue
        d_start = max(domain["start"] - 1, start)
        d_end = min(domain["end"], end)
        x = left_margin + (d_start - start) * char_width
        width = max(char_width, (d_end - d_start) * char_width)
        drawing.add(drawing.rect(
            insert=(x, top_y), size=(width, 13),
            fill=domain["color"], fill_opacity=0.85,
            stroke=domain["color"], stroke_width=1, rx=3, ry=3,
        ))
        label = domain.get("label")
        if label and width > len(label) * 5.5:
            drawing.add(drawing.text(
                label,
                insert=(x + width / 2.0, top_y + 9.5),
                fill="#ffffff", font_size=8.5,
                font_family="Helvetica, Arial, sans-serif",
                font_weight="700", text_anchor="middle",
            ))


def _char_color(char, track_kind, palette):
    if track_kind == "sequence":
        return palette["sequence"]
    return palette.get(char, palette["sequence"])


def _mix_hex(low_hex, high_hex, weight):
    weight = max(0.0, min(1.0, weight))
    low = _hex_to_rgb(low_hex)
    high = _hex_to_rgb(high_hex)
    mixed = tuple(int(low[i] + (high[i] - low[i]) * weight) for i in range(3))
    return "#%02x%02x%02x" % mixed


def _hex_to_rgb(value):
    value = value.lstrip("#")
    return tuple(int(value[i:i + 2], 16) for i in (0, 2, 4))
