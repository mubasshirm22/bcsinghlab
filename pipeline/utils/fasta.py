"""
FASTA / sequence input validation and normalization.
No external dependencies — uses only stdlib.
"""

import re

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBZXUO")  # standard + ambiguous IUPAC protein


def parse_input(raw: str) -> dict:
    """
    Accept either:
      1. A raw amino-acid string (no > header)
      2. A FASTA block (one or more lines, first line may start with >)

    Returns:
        {
          "ok": True/False,
          "sequence": str,      # uppercased, whitespace stripped (if ok)
          "header": str,        # FASTA header without >, empty string if raw input
          "error": str          # human-readable message (if not ok)
        }
    """
    raw = raw.strip()
    if not raw:
        return {"ok": False, "sequence": "", "header": "", "error": "No input provided."}

    lines = raw.splitlines()
    header = ""
    seq_lines = []

    if lines[0].startswith(">"):
        header = lines[0][1:].strip()
        seq_lines = lines[1:]
    else:
        seq_lines = lines

    sequence = "".join(seq_lines).replace(" ", "").replace("\t", "").upper()

    if not sequence:
        return {"ok": False, "sequence": "", "header": header,
                "error": "No sequence characters found after the header."}

    # Check for DNA/RNA — if >80 % of chars are ACGTU it's almost certainly nucleotide
    nucleotide_chars = set("ACGTU")
    nuc_count = sum(1 for c in sequence if c in nucleotide_chars)
    if len(sequence) > 0 and (nuc_count / len(sequence)) > 0.80:
        return {"ok": False, "sequence": "", "header": header,
                "error": "Input looks like a nucleotide sequence. Please provide a protein sequence."}

    # Check all characters are valid amino acids
    invalid = set(c for c in sequence if c not in VALID_AA)
    if invalid:
        return {"ok": False, "sequence": "", "header": header,
                "error": f"Invalid characters in sequence: {', '.join(sorted(invalid))}"}

    if len(sequence) < 10:
        return {"ok": False, "sequence": "", "header": header,
                "error": "Sequence is too short (minimum 10 residues)."}

    if len(sequence) > 5000:
        return {"ok": False, "sequence": "", "header": header,
                "error": "Sequence is too long (maximum 5000 residues for this tool)."}

    return {"ok": True, "sequence": sequence, "header": header, "error": ""}


def is_accession(text: str) -> bool:
    """Return True if the text looks like a UniProt or NCBI accession, not a FASTA."""
    text = text.strip()
    if "\n" in text or text.startswith(">"):
        return False
    # UniProt: P12345, Q9Y2H9, A0A000XXX, etc.
    if re.match(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$', text):
        return True
    # NCBI protein: NP_123456.1, XP_123456789.2, WP_012345678.1, AAB12345.1
    if re.match(r'^[A-Z]{2}_\d{6,9}(\.\d+)?$', text) or re.match(r'^[A-Z]{2,3}\d{5,8}(\.\d+)?$', text):
        return True
    return False
