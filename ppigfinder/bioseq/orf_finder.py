#!/usr/bin/env python3
"""
ORF prediction routines used by ppigFinder.

This module contains sequence-level ORF prediction logic and is independent
from the PyQt interface.
"""

from __future__ import annotations

from .genetic_code import START_CODONS_BACTERIAL, STOP_CODONS_STANDARD
from .sequence import gc_content, reverse_complement, translate_dna


try:
    import pyrodigal

    PYRODIGAL_AVAILABLE = True
except Exception:
    pyrodigal = None
    PYRODIGAL_AVAILABLE = False


def find_orfs(
    dna_sequence: str,
    min_aa: int = 30,
    start_codons: set[str] | None = None,
) -> list[dict]:
    """
    Find ORFs using a simple six-frame start-to-stop codon scanner.

    Coordinates are returned in genomic 0-based half-open format:
    start inclusive, end exclusive-like as used by the legacy application.
    """
    if start_codons is None:
        start_codons = set(START_CODONS_BACTERIAL)

    stop_codons = set(STOP_CODONS_STANDARD)
    min_len = min_aa * 3
    dna_sequence = dna_sequence.upper()

    orfs: list[dict] = []

    for frame in range(3):
        for strand_seq, strand_name in [
            (dna_sequence, "+"),
            (reverse_complement(dna_sequence), "-"),
        ]:
            i = frame

            while i < len(strand_seq) - 2:
                codon = strand_seq[i:i + 3]

                if codon in start_codons:
                    j = i + 3

                    while j < len(strand_seq):
                        if strand_seq[j:j + 3] in stop_codons:
                            length = j + 3 - i

                            if length >= min_len:
                                dna = strand_seq[i:j + 3]
                                protein = translate_dna(dna)

                                if strand_name == "+":
                                    start = i
                                    end = j + 3
                                else:
                                    start = len(dna_sequence) - (j + 3)
                                    end = len(dna_sequence) - i

                                orfs.append({
                                    "frame": frame + (3 if strand_name == "-" else 0),
                                    "strand": strand_name,
                                    "start": start,
                                    "end": end,
                                    "dna": dna,
                                    "protein": protein,
                                    "length": length,
                                    "gc": gc_content(dna),
                                    "domains": [],
                                    "neighborhood": [],
                                    "candidate_score": 0.0,
                                    "source": "6frame",
                                })

                            i = j
                            break

                        j += 3

                i += 3

    orfs.sort(key=lambda item: item["start"])
    return orfs


def find_orfs_pyrodigal(
    dna_sequence: str,
    meta: bool = True,
    min_aa: int = 30,
    closed_ends: bool = False,
    translation_table: int = 11,
    mask: bool = False,
) -> list[dict]:
    """
    Find ORFs using Pyrodigal/Prodigal.

    Pyrodigal is the preferred bacterial gene caller when available.
    """
    if not PYRODIGAL_AVAILABLE:
        raise ImportError(
            "Pyrodigal is not installed.\n\n"
            "Install with:\n"
            "  pip install pyrodigal\n\n"
            "Or use conda:\n"
            "  conda install -c bioconda pyrodigal"
        )

    seq = dna_sequence.upper()
    encoded_seq = seq.encode() if isinstance(seq, str) else seq

    if meta:
        gene_finder = pyrodigal.GeneFinder(
            meta=True,
            closed=closed_ends,
            min_gene=min_aa * 3,
            mask=mask,
        )
    else:
        gene_finder = pyrodigal.GeneFinder(
            meta=False,
            closed=closed_ends,
            min_gene=min_aa * 3,
            mask=mask,
        )
        gene_finder.train(encoded_seq, translation_table=translation_table)

    genes = gene_finder.find_genes(encoded_seq)
    orfs: list[dict] = []

    for gene in genes:
        start_0 = gene.begin - 1
        end_0 = gene.end
        strand = "+" if gene.strand == 1 else "-"

        if strand == "+":
            frame = start_0 % 3
            dna_sub = seq[start_0:end_0]
        else:
            frame = (len(seq) - end_0) % 3 + 3
            dna_sub = reverse_complement(seq[start_0:end_0])

        protein = translate_dna(dna_sub)

        try:
            cscore = gene.confidence()
        except (AttributeError, TypeError):
            cscore = gene.cscore if hasattr(gene, "cscore") else 0.0

        try:
            rbs_motif = gene.rbs_motif
        except AttributeError:
            rbs_motif = None

        orfs.append({
            "frame": frame,
            "strand": strand,
            "start": start_0,
            "end": end_0,
            "dna": dna_sub,
            "protein": protein,
            "length": end_0 - start_0,
            "gc": gc_content(dna_sub),
            "domains": [],
            "neighborhood": [],
            "candidate_score": 0.0,
            "source": "pyrodigal",
            "pyrodigal_score": round(cscore, 2) if cscore else 0.0,
            "rbs_motif": rbs_motif or "",
            "partial": (
                getattr(gene, "partial_begin", False)
                or getattr(gene, "partial_end", False)
            ),
        })

    orfs.sort(key=lambda item: item["start"])
    return orfs


def find_orfs_hybrid(
    dna_sequence: str,
    min_aa: int = 30,
    start_codons: set[str] | None = None,
    pyro_meta: bool = True,
    pyro_min_aa: int = 30,
    pyro_closed: bool = False,
    pyro_translation_table: int = 11,
    pyro_mask: bool = False,
    pyro_start_filter: dict | None = None,
) -> list[dict]:
    """
    Hybrid gene finder.

    Strategy:
    1. Run Pyrodigal as primary caller.
    2. Identify genome regions not covered by Pyrodigal genes.
    3. Run six-frame scanning only in uncovered gaps.
    4. Merge both sets and sort by genomic coordinate.
    """
    if start_codons is None:
        start_codons = set(START_CODONS_BACTERIAL)

    pyro_orfs = find_orfs_pyrodigal(
        dna_sequence,
        meta=pyro_meta,
        min_aa=pyro_min_aa,
        closed_ends=pyro_closed,
        translation_table=pyro_translation_table,
        mask=pyro_mask,
    )

    if pyro_start_filter and not pyro_start_filter.get("all", True):
        allowed = {
            codon
            for codon in ("ATG", "GTG", "TTG")
            if pyro_start_filter.get(codon)
        }

        if allowed:
            pyro_orfs = [
                orf
                for orf in pyro_orfs
                if orf.get("dna", "")[:3].upper() in allowed
            ]

    seq_len = len(dna_sequence)
    intervals = sorted((orf["start"], orf["end"]) for orf in pyro_orfs)

    merged: list[list[int]] = []

    for seg_start, seg_end in intervals:
        if merged and seg_start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], seg_end)
        else:
            merged.append([seg_start, seg_end])

    gaps: list[tuple[int, int]] = []
    previous_end = 0

    for seg_start, seg_end in merged:
        if seg_start > previous_end:
            gaps.append((previous_end, seg_start))
        previous_end = max(previous_end, seg_end)

    if previous_end < seq_len:
        gaps.append((previous_end, seq_len))

    automatic_orfs: list[dict] = []
    min_gap = min_aa * 3 + 3

    for gap_start, gap_end in gaps:
        if gap_end - gap_start < min_gap:
            continue

        subseq = dna_sequence[gap_start:gap_end]
        gap_orfs = find_orfs(
            subseq,
            min_aa=min_aa,
            start_codons=start_codons,
        )

        for orf in gap_orfs:
            orf["start"] += gap_start
            orf["end"] += gap_start
            orf["source"] = "automatic"
            automatic_orfs.append(orf)

    all_orfs = pyro_orfs + automatic_orfs
    all_orfs.sort(key=lambda item: item["start"])

    return all_orfs
