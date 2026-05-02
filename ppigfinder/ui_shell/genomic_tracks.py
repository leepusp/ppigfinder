#!/usr/bin/env python3
"""
Genomic signal-track helpers for the guided UI shell.

Tracks that can be computed from a genome sequence:
- GC content
- GC skew: (G - C) / (G + C)

Experimental tracks such as DNA-Seq/RNA-Seq/proteomics coverage require
external files and should be integrated later through bedGraph/WIG/TSV/BAM-style
importers.
"""

from __future__ import annotations

from pathlib import Path
import re


def read_text_any_encoding(path: str | Path) -> str:
    path = Path(path)

    for encoding in ("utf-8", "latin-1", "cp1252"):
        try:
            return path.read_text(encoding=encoding)
        except UnicodeDecodeError:
            continue

    return path.read_text(errors="ignore")


def read_genome_sequence(path: str | Path) -> str:
    """
    Read a main nucleotide sequence from FASTA or GenBank-like files.
    For multi-FASTA, returns the longest sequence.
    """
    path = Path(path)
    text = read_text_any_encoding(path)
    lower = path.suffix.lower()

    if lower in {".gb", ".gbk", ".genbank"} or "ORIGIN" in text[:20000]:
        return _read_genbank_origin(text)

    return _read_fasta_longest(text)


def _read_fasta_longest(text: str) -> str:
    records = []
    chunks = []

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if chunks:
                records.append("".join(chunks))
                chunks = []
            continue

        chunks.append(re.sub(r"[^A-Za-z]", "", line))

    if chunks:
        records.append("".join(chunks))

    if not records:
        seq = re.sub(r"[^A-Za-z]", "", text)
        return seq.upper()

    return max(records, key=len).upper()


def _read_genbank_origin(text: str) -> str:
    idx = text.find("ORIGIN")

    if idx < 0:
        return ""

    origin = text[idx:]
    end_idx = origin.find("//")
    if end_idx >= 0:
        origin = origin[:end_idx]

    seq = re.sub(r"[^A-Za-z]", "", origin.replace("ORIGIN", ""))
    seq = seq.replace("END", "")

    return seq.upper()


def compute_gc_and_skew_bins(sequence: str, start: int, end: int, bins: int = 220) -> list[dict]:
    """
    Compute GC content and GC skew over the visible genomic interval.

    Coordinates are 1-based inclusive.
    """
    sequence = re.sub(r"[^ACGTNacgtn]", "", sequence or "").upper()

    if not sequence:
        return []

    start = max(1, int(start))
    end = min(len(sequence), int(end))

    if end <= start:
        return []

    window = sequence[start - 1:end]
    window_len = len(window)

    bins = max(20, min(int(bins), window_len))
    step = max(1, window_len // bins)

    rows = []

    for offset in range(0, window_len, step):
        chunk = window[offset:offset + step]

        if not chunk:
            continue

        g = chunk.count("G")
        c = chunk.count("C")
        a = chunk.count("A")
        t = chunk.count("T")
        valid = a + c + g + t

        gc = ((g + c) / valid) if valid else 0.0
        skew = ((g - c) / (g + c)) if (g + c) else 0.0

        rows.append(
            {
                "start": start + offset,
                "end": min(end, start + offset + len(chunk) - 1),
                "gc": gc,
                "skew": skew,
            }
        )

    return rows
