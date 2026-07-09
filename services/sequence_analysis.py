"""Biophysical and secondary-structure analytics for a protein sequence.

Everything here is dependency-light: Biopython is used when available for the
richer ProtParam metrics, but every function degrades gracefully to a pure
Python implementation so the analysis panel never hard-fails a results page.
"""

from __future__ import annotations

import re
from typing import Iterable

try:  # Biopython is in requirements, but never let its absence break a page.
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
except Exception:  # pragma: no cover - defensive import
    ProteinAnalysis = None


# Kyte & Doolittle (1982) hydropathy index.
KD_HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Eisenberg consensus scale — used for the hydrophobic-moment calculation.
EISENBERG = {
    "A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29,
    "Q": -0.85, "E": -0.74, "G": 0.48, "H": -0.40, "I": 1.38,
    "L": 1.06, "K": -1.50, "M": 0.64, "F": 1.19, "P": 0.12,
    "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08,
}

AA_NAMES = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
}

# Physico-chemical groupings used for the composition summary.
AA_GROUPS = {
    "Hydrophobic (A,V,L,I,M,F,W,P)": set("AVLIMFWP"),
    "Polar uncharged (S,T,N,Q,C,G,Y)": set("STNQCGY"),
    "Positively charged (K,R,H)": set("KRH"),
    "Negatively charged (D,E)": set("DE"),
}

_VALID_AA = set(KD_HYDROPATHY.keys())


def clean_sequence(seq: str) -> str:
    """Uppercase and strip anything that is not a standard amino-acid letter."""
    if not seq:
        return ""
    return re.sub(r"[^A-Z]", "", (seq or "").upper())


def hydropathy_profile(seq: str, window: int = 9) -> dict:
    """Kyte-Doolittle sliding-window hydropathy profile.

    Returns per-window values centred on each residue plus summary numbers that
    are useful for spotting transmembrane spans (peaks above ~1.6).
    """
    residues = clean_sequence(seq)
    n = len(residues)
    if n == 0:
        return {"window": window, "values": [], "positions": [], "max": 0.0, "min": 0.0, "gravy": 0.0}

    window = max(1, min(int(window or 9), max(1, n)))
    if window % 2 == 0:  # keep it odd so the window is symmetric around a residue
        window += 1
    half = window // 2

    values: list[float] = []
    positions: list[int] = []
    for center in range(n):
        start = max(0, center - half)
        end = min(n, center + half + 1)
        window_res = residues[start:end]
        score = sum(KD_HYDROPATHY.get(aa, 0.0) for aa in window_res) / len(window_res)
        values.append(round(score, 3))
        positions.append(center + 1)

    gravy = sum(KD_HYDROPATHY.get(aa, 0.0) for aa in residues) / n
    return {
        "window": window,
        "values": values,
        "positions": positions,
        "max": round(max(values), 3),
        "min": round(min(values), 3),
        "gravy": round(gravy, 3),
        # A crude but handy flag: sustained hydrophobic peaks suggest TM helices.
        "likely_tm_spans": _hydrophobic_spans(values, positions, threshold=1.6, min_len=15),
    }


def _hydrophobic_spans(values, positions, threshold=1.6, min_len=15):
    spans = []
    run_start = None
    for idx, value in enumerate(values):
        if value >= threshold:
            if run_start is None:
                run_start = idx
        else:
            if run_start is not None and idx - run_start >= min_len:
                spans.append({"start": positions[run_start], "end": positions[idx - 1]})
            run_start = None
    if run_start is not None and len(values) - run_start >= min_len:
        spans.append({"start": positions[run_start], "end": positions[-1]})
    return spans


def hydrophobic_moment(seq: str, angle: float = 100.0) -> float:
    """Eisenberg mean hydrophobic moment (µH) for the whole peptide.

    100° is the canonical alpha-helix angle; high µH indicates amphipathicity.
    """
    residues = clean_sequence(seq)
    if not residues:
        return 0.0
    import math

    sin_sum = 0.0
    cos_sum = 0.0
    for i, aa in enumerate(residues):
        h = EISENBERG.get(aa, 0.0)
        theta = math.radians(angle * i)
        sin_sum += h * math.sin(theta)
        cos_sum += h * math.cos(theta)
    return round(math.sqrt(sin_sum ** 2 + cos_sum ** 2) / len(residues), 3)


def composition(seq: str) -> dict:
    """Per-residue counts/percentages plus physico-chemical group breakdown."""
    residues = clean_sequence(seq)
    n = len(residues)
    counts = {aa: 0 for aa in sorted(_VALID_AA)}
    for aa in residues:
        if aa in counts:
            counts[aa] += 1

    per_residue = [
        {
            "code": aa,
            "name": AA_NAMES.get(aa, aa),
            "count": counts[aa],
            "percent": round(100.0 * counts[aa] / n, 1) if n else 0.0,
        }
        for aa in sorted(counts.keys())
    ]

    groups = []
    for label, members in AA_GROUPS.items():
        count = sum(counts[aa] for aa in members if aa in counts)
        groups.append({
            "label": label,
            "count": count,
            "percent": round(100.0 * count / n, 1) if n else 0.0,
        })

    return {"length": n, "per_residue": per_residue, "groups": groups}


def biophysical_stats(seq: str) -> dict:
    """Molecular weight, pI, GRAVY, aromaticity, instability, extinction coeff.

    Uses Biopython's ProtParam when available; otherwise fills a minimal subset.
    """
    residues = clean_sequence(seq)
    n = len(residues)
    if n == 0:
        return {"length": 0}

    stats = {"length": n}
    if ProteinAnalysis is not None:
        try:
            analysis = ProteinAnalysis(residues)
            stats["molecular_weight"] = round(analysis.molecular_weight(), 1)
            stats["isoelectric_point"] = round(analysis.isoelectric_point(), 2)
            stats["gravy"] = round(analysis.gravy(), 3)
            stats["aromaticity"] = round(analysis.aromaticity(), 3)
            instability = analysis.instability_index()
            stats["instability_index"] = round(instability, 2)
            stats["stable"] = instability < 40.0
            ext_reduced, ext_cystine = analysis.molar_extinction_coefficient()
            stats["extinction_reduced"] = int(ext_reduced)
            stats["extinction_cystines"] = int(ext_cystine)
            charge = analysis.charge_at_pH(7.0)
            stats["net_charge_ph7"] = round(charge, 2)
            return stats
        except Exception:
            pass  # fall through to the manual estimates

    # Minimal manual fallback (no Biopython or it errored on a rare residue).
    avg_weights = {
        "A": 71.08, "R": 156.19, "N": 114.10, "D": 115.09, "C": 103.14,
        "Q": 128.13, "E": 129.12, "G": 57.05, "H": 137.14, "I": 113.16,
        "L": 113.16, "K": 128.17, "M": 131.19, "F": 147.18, "P": 97.12,
        "S": 87.08, "T": 101.10, "W": 186.21, "Y": 163.18, "V": 99.13,
    }
    weight = sum(avg_weights.get(aa, 110.0) for aa in residues) + 18.02
    stats["molecular_weight"] = round(weight, 1)
    stats["gravy"] = round(sum(KD_HYDROPATHY.get(aa, 0.0) for aa in residues) / n, 3)
    aromatic = sum(1 for aa in residues if aa in "FWY")
    stats["aromaticity"] = round(aromatic / n, 3)
    return stats


def ss_content(track: str) -> dict:
    """Percent helix / strand / coil for a three-state prediction string."""
    if not track:
        return {"helix": 0.0, "strand": 0.0, "coil": 0.0, "counts": {"H": 0, "E": 0, "C": 0}, "total": 0}
    counts = {"H": 0, "E": 0, "C": 0}
    total = 0
    for char in track.upper():
        if char in counts:
            counts[char] += 1
            total += 1
    if total == 0:
        return {"helix": 0.0, "strand": 0.0, "coil": 0.0, "counts": counts, "total": 0}
    return {
        "helix": round(100.0 * counts["H"] / total, 1),
        "strand": round(100.0 * counts["E"] / total, 1),
        "coil": round(100.0 * counts["C"] / total, 1),
        "counts": counts,
        "total": total,
    }


def q3_agreement(prediction: str, reference: str) -> dict:
    """Q3 accuracy of a prediction string against a reference (e.g. PDB DSSP).

    Q3 = fraction of positions where the three-state assignment matches, scored
    only over positions where both strings carry a valid H/E/C state.
    """
    if not prediction or not reference:
        return {"q3": None, "matched": 0, "compared": 0}
    compared = 0
    matched = 0
    for pred_char, ref_char in zip(prediction.upper(), reference.upper()):
        if pred_char in "HEC" and ref_char in "HEC":
            compared += 1
            if pred_char == ref_char:
                matched += 1
    if compared == 0:
        return {"q3": None, "matched": 0, "compared": 0}
    return {"q3": round(100.0 * matched / compared, 1), "matched": matched, "compared": compared}


def analyze_row(row: dict) -> dict:
    """Full analysis bundle for a SSPred database row.

    Pulls the sequence, every completed predictor track, the consensus and any
    known PDB structure, and returns hydropathy + biophysics + SS content + Q3.
    """
    sequence = (row.get("seq") or "").strip()
    predictors = _completed_predictors(row)
    consensus = row.get("majorityvote") or ""
    pdb_secondary = _pdb_secondary(row.get("pdb"))

    ss_summaries = []
    if consensus:
        ss_summaries.append({"name": "Consensus", "content": ss_content(consensus)})
    for predictor in predictors:
        ss_summaries.append({"name": predictor["label"], "content": ss_content(predictor["pred"])})

    accuracy = []
    if pdb_secondary:
        if consensus:
            accuracy.append({"name": "Consensus", "score": q3_agreement(consensus, pdb_secondary)})
        for predictor in predictors:
            accuracy.append({"name": predictor["label"], "score": q3_agreement(predictor["pred"], pdb_secondary)})
        accuracy.sort(key=lambda item: (item["score"]["q3"] is None, -(item["score"]["q3"] or 0)))

    return {
        "length": len(clean_sequence(sequence)),
        "has_sequence": bool(sequence),
        "has_reference": bool(pdb_secondary),
        "hydropathy": hydropathy_profile(sequence),
        "hydrophobic_moment": hydrophobic_moment(sequence),
        "stats": biophysical_stats(sequence),
        "composition": composition(sequence),
        "ss_content": ss_summaries,
        "accuracy": accuracy,
    }


def _completed_predictors(row: dict) -> list[dict]:
    available = []
    for key, value in row.items():
        if key.endswith("pred"):
            name = key[:-4]
            if name in {"pdb", "majorityvote"}:
                continue
            status = row.get(name + "stat")
            pred = row.get(name + "pred")
            if status in (1, 3) and pred:
                available.append({"name": name, "label": name.upper(), "pred": pred})
    return available


def _pdb_secondary(value) -> str:
    if not value:
        return ""
    if isinstance(value, dict):
        return value.get("secondary") or ""
    try:
        import json
        data = json.loads(value)
        if isinstance(data, dict):
            return data.get("secondary") or ""
    except Exception:
        return ""
    return ""
