#!/usr/bin/env python3
"""
Lightweight backend operations for the experimental guided UI shell.

This module intentionally avoids the legacy GUI. It gives the guided shell
real local operations while the full service-layer binding is still evolving.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re


@dataclass(slots=True)
class GuidedORF:
    id: str
    start: int
    end: int
    strand: str
    frame: int
    nt_sequence: str
    protein_sequence: str

    @property
    def aa_length(self) -> int:
        return len(self.protein_sequence)


@dataclass(slots=True)
class ORFPredictionSummary:
    source_file: str
    sequence_name: str
    sequence_length: int
    min_aa: int
    orf_count: int
    longest_orf_aa: int
    shortest_orf_aa: int
    orfs: list[GuidedORF]


GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

START_CODONS = {"ATG", "GTG", "TTG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}


def read_text(path: str | Path) -> str:
    path = Path(path)

    for encoding in ("utf-8", "latin-1", "cp1252"):
        try:
            return path.read_text(encoding=encoding)
        except UnicodeDecodeError:
            continue

    return path.read_text(errors="ignore")


def read_first_fasta_sequence(path: str | Path) -> tuple[str, str]:
    """
    Read the first/longest sequence from a FASTA-like file.
    """
    path = Path(path)
    text = read_text(path)

    records: list[tuple[str, str]] = []
    header = path.stem
    seq_chunks: list[str] = []

    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue

        if line.startswith(">"):
            if seq_chunks:
                records.append((header, "".join(seq_chunks)))
                seq_chunks = []
            header = line[1:].strip() or path.stem
        else:
            seq_chunks.append(re.sub(r"[^A-Za-z]", "", line))

    if seq_chunks:
        records.append((header, "".join(seq_chunks)))

    if not records:
        sequence = re.sub(r"[^A-Za-z]", "", text)
        return path.stem, sequence.upper()

    name, sequence = max(records, key=lambda item: len(item[1]))
    return name, sequence.upper().replace("U", "T")


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(table)[::-1].upper()


def translate_nt(sequence: str) -> str:
    protein = []

    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3].upper()
        protein.append(GENETIC_CODE.get(codon, "X"))

    return "".join(protein)


def _scan_orfs_on_strand(
    sequence: str,
    strand: str,
    original_length: int,
    min_aa: int,
    prefix: str,
) -> list[GuidedORF]:
    orfs: list[GuidedORF] = []
    counter = 1

    for frame in range(3):
        i = frame

        while i <= len(sequence) - 3:
            codon = sequence[i:i + 3]

            if codon not in START_CODONS:
                i += 3
                continue

            j = i + 3

            while j <= len(sequence) - 3:
                stop = sequence[j:j + 3]

                if stop in STOP_CODONS:
                    nt = sequence[i:j + 3]
                    protein = translate_nt(nt).rstrip("*")

                    if len(protein) >= min_aa:
                        if strand == "+":
                            start = i + 1
                            end = j + 3
                        else:
                            start = original_length - (j + 3) + 1
                            end = original_length - i

                        orfs.append(
                            GuidedORF(
                                id=f"{prefix}_{counter:05d}",
                                start=start,
                                end=end,
                                strand=strand,
                                frame=frame,
                                nt_sequence=nt,
                                protein_sequence=protein,
                            )
                        )
                        counter += 1

                    i = j + 3
                    break

                j += 3
            else:
                i += 3

    return orfs


def predict_orfs_from_file(path: str | Path, min_aa: int = 30) -> ORFPredictionSummary:
    """
    Lightweight six-frame ORF scan used by the guided shell.
    """
    name, sequence = read_first_fasta_sequence(path)
    sequence = re.sub(r"[^ACGTN]", "", sequence.upper())

    if not sequence:
        return ORFPredictionSummary(
            source_file=str(path),
            sequence_name=name,
            sequence_length=0,
            min_aa=min_aa,
            orf_count=0,
            longest_orf_aa=0,
            shortest_orf_aa=0,
            orfs=[],
        )

    plus = _scan_orfs_on_strand(sequence, "+", len(sequence), min_aa, "guided_plus")
    minus = _scan_orfs_on_strand(reverse_complement(sequence), "-", len(sequence), min_aa, "guided_minus")

    orfs = sorted(plus + minus, key=lambda item: (item.start, item.end, item.strand))
    lengths = [orf.aa_length for orf in orfs]

    return ORFPredictionSummary(
        source_file=str(path),
        sequence_name=name,
        sequence_length=len(sequence),
        min_aa=min_aa,
        orf_count=len(orfs),
        longest_orf_aa=max(lengths) if lengths else 0,
        shortest_orf_aa=min(lengths) if lengths else 0,
        orfs=orfs,
    )


def write_orfs_fasta(path: str | Path, orfs: list[GuidedORF]) -> None:
    path = Path(path)

    with path.open("w", encoding="utf-8") as handle:
        for orf in orfs:
            header = (
                f">{orf.id} start={orf.start} end={orf.end} "
                f"strand={orf.strand} frame={orf.frame} aa_len={orf.aa_length}"
            )
            handle.write(header + "\n")

            seq = orf.protein_sequence
            for i in range(0, len(seq), 70):
                handle.write(seq[i:i + 70] + "\n")


def build_guided_summary_markdown(state: dict) -> str:
    """
    Build a simple Markdown report for the guided shell state.
    """
    lines = [
        "# ppigFinder Guided Shell Summary",
        "",
        "## Input state",
        "",
        f"- Genome file: {state.get('genome_file') or 'not selected'}",
        f"- Project file: {state.get('project_file') or 'not selected'}",
        f"- Snapshot file: {state.get('snapshot_file') or 'not selected'}",
        f"- AF3 results folder: {state.get('af3_results_folder') or 'not selected'}",
        "",
        "## Genome metadata",
        "",
        f"- Genome name: {state.get('genome_name') or 'N/A'}",
        f"- File type: {state.get('genome_file_type') or 'N/A'}",
        f"- Sequence count: {state.get('genome_sequence_count') or 'N/A'}",
        f"- Total length: {state.get('genome_total_length') or 'N/A'}",
        f"- Longest sequence length: {state.get('genome_longest_length') or 'N/A'}",
        f"- GC%: {state.get('genome_gc_percent') if state.get('genome_gc_percent') is not None else 'N/A'}",
        f"- Validation: {'OK' if state.get('genome_valid') else 'pending/problem'}",
        "",
        "## ORF prediction",
        "",
        f"- ORFs predicted: {state.get('guided_orf_count') or 0}",
        f"- Longest ORF length: {state.get('guided_longest_orf_aa') or 0} aa",
        f"- Shortest ORF length: {state.get('guided_shortest_orf_aa') or 0} aa",
        "",
        "## Annotation",
        "",
        f"- Candidate ORFs available: {state.get('guided_annotation_candidates_count') or state.get('guided_orf_count') or 0}",
        f"- BLAST selected/planned: {'yes' if state.get('guided_blast_planned') else 'no'}",
        f"- HMM selected/planned: {'yes' if state.get('guided_hmm_planned') else 'no'}",
        f"- Neighbourhood selected/planned: {'yes' if state.get('guided_neighborhood_planned') else 'no'}",
        "",
        "## DaVinci / HPC",
        "",
        f"- HPC profile: {state.get('hpc_profile') or 'DaVinci'}",
        f"- HPC host: {state.get('hpc_host') or 'not configured'}",
        f"- HPC status: {state.get('hpc_status') or 'not tested'}",
        f"- HPC mode: {state.get('hpc_mode') or 'optional'}",
        "",
        "## Reports",
        "",
        "This summary was generated from the experimental guided shell.",
        "",
    ]

    return "\n".join(lines)


def write_guided_summary(path: str | Path, state: dict) -> None:
    Path(path).write_text(build_guided_summary_markdown(state), encoding="utf-8")
