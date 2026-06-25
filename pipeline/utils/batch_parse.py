"""
Flexible parsing of batch submission input for ProtSuite.

Accepts pasted text (multi-FASTA, or a comma/newline/whitespace-separated
list of accessions) plus any number of uploaded files (each parsed the same
way), and returns one input string per sequence/accession to submit as its
own pipeline job.
"""

import re

MAX_BATCH_SIZE = 25


class BatchTooLargeError(Exception):
    pass


def _split_fasta_records(text: str) -> list:
    """Split multi-FASTA text into individual '>header\\nSEQUENCE' records."""
    records = []
    chunks = re.split(r'(?=^>)', text, flags=re.MULTILINE)
    for chunk in chunks:
        chunk = chunk.strip()
        if chunk:
            records.append(chunk)
    return records


def _split_accession_list(text: str) -> list:
    """Split a comma/newline/whitespace-separated list of accessions/tokens."""
    tokens = re.split(r'[,\n\r\t ]+', text.strip())
    return [t.strip() for t in tokens if t.strip()]


def parse_text_block(text: str) -> list:
    """Parse one block of raw text (pasted or from a file) into entries."""
    text = (text or "").strip()
    if not text:
        return []
    if '>' in text:
        return _split_fasta_records(text)
    return _split_accession_list(text)


def parse_batch_input(raw_text: str = "", file_contents: list = None) -> list:
    """
    Combine pasted text and uploaded file contents into a flat list of
    submission entries (each becomes one ProtSuite job). Raises
    BatchTooLargeError if the combined entry count exceeds MAX_BATCH_SIZE.
    """
    entries = []
    entries.extend(parse_text_block(raw_text))
    for content in (file_contents or []):
        entries.extend(parse_text_block(content))

    if len(entries) > MAX_BATCH_SIZE:
        raise BatchTooLargeError(
            f"Batch contains {len(entries)} entries, which exceeds the limit of {MAX_BATCH_SIZE}."
        )
    return entries
