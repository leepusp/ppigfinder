#!/usr/bin/env python3
"""
AlphaFold/AF3 table export utilities.
"""

from __future__ import annotations

from pathlib import Path
import csv


AF3_RESULT_COLUMNS = [
    "name",
    "classification",
    "iptm",
    "ptm",
    "ranking_score",
    "pae_inter",
    "pae_min",
    "cp_iptm",
    "contact_percent",
    "result_dir",
    "summary_file",
    "confidence_file",
    "ranking_file",
    "model_file",
]


def write_af3_results_table(
    path: str | Path,
    results: list[dict],
    delimiter: str = "\t",
) -> None:
    """
    Write AF3 results as TSV/CSV.
    """
    path = Path(path)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=AF3_RESULT_COLUMNS,
            delimiter=delimiter,
            extrasaction="ignore",
        )
        writer.writeheader()

        for item in results or []:
            if isinstance(item, dict):
                writer.writerow(item)
