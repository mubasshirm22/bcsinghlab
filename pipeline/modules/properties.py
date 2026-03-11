"""
Basic sequence properties module.

Uses only Biopython (pip install biopython) — no external HTTP calls.
All calculations are local and instant.

Properties computed:
  - length
  - amino acid composition (count + %)
  - molecular weight estimate (monoisotopic average)
  - isoelectric point estimate (Henderson-Hasselbalch via Biopython)
  - instability index (Guruprasad 1990, via Biopython)
  - GRAVY score (grand average of hydropathicity, Kyte-Doolittle)
  - aromaticity
  - secondary structure fraction estimates (Biopython ProtParam)
"""

try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    _BIOPYTHON = True
except ImportError:
    _BIOPYTHON = False


def run(sequence: str, job_dir: str) -> dict:
    """
    Args:
        sequence: clean uppercase amino acid string
        job_dir:  unused here but kept for uniform module signature

    Returns:
        {"status": "ok"|"error", "data": {...}, "error": str}
    """
    if not _BIOPYTHON:
        return {"status": "error", "data": {}, "error": "Biopython not installed (pip install biopython)"}

    # Biopython ProteinAnalysis does not accept ambiguous residues (B, Z, X, U, O)
    # Replace them with A (alanine) for estimation purposes only
    clean = _sanitise(sequence)

    try:
        analysis = ProteinAnalysis(clean)

        mw = round(analysis.molecular_weight(), 2)
        pi = round(analysis.isoelectric_point(), 2)

        try:
            instability = round(analysis.instability_index(), 2)
            stable = instability < 40
        except Exception:
            instability = None
            stable = None

        try:
            gravy = round(analysis.gravy(), 4)
        except Exception:
            gravy = None

        try:
            aromaticity = round(analysis.aromaticity(), 4)
        except Exception:
            aromaticity = None

        # Secondary structure fraction: helix, turn, sheet
        try:
            ss_fracs = analysis.secondary_structure_fraction()
            helix_frac = round(ss_fracs[0], 4)
            turn_frac  = round(ss_fracs[1], 4)
            sheet_frac = round(ss_fracs[2], 4)
        except Exception:
            helix_frac = turn_frac = sheet_frac = None

        aa_count = analysis.count_amino_acids()
        aa_percent = {aa: round(count / len(clean) * 100, 2)
                      for aa, count in aa_count.items() if count > 0}

        data = {
            "length": len(sequence),
            "molecular_weight_da": mw,
            "isoelectric_point": pi,
            "instability_index": instability,
            "is_stable": stable,
            "gravy": gravy,
            "aromaticity": aromaticity,
            "helix_fraction": helix_frac,
            "turn_fraction": turn_frac,
            "sheet_fraction": sheet_frac,
            "amino_acid_count": aa_count,
            "amino_acid_percent": aa_percent,
            "ambiguous_residues_replaced": len(sequence) - len(clean),
        }
        return {"status": "ok", "data": data, "error": ""}

    except Exception as e:
        return {"status": "error", "data": {}, "error": f"ProtParam calculation failed: {e}"}


def _sanitise(seq: str) -> str:
    """Replace ambiguous IUPAC residues with alanine for ProtParam compatibility."""
    mapping = str.maketrans("BZXUO", "AAAAA")
    return seq.translate(mapping)
