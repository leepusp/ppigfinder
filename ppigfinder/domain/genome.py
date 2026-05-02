#!/usr/bin/env python3
"""
Domain models related to genome and sequence records.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(slots=True)
class GenomeRecord:
    """
    Loaded genome or nucleotide sequence record.
    """

    name: str
    sequence: str
    source_path: Path | None = None
    topology: str = "linear"
    features: list[dict] = field(default_factory=list)
    primers: list[dict] = field(default_factory=list)
    notes: dict = field(default_factory=dict)

    @property
    def length(self) -> int:
        return len(self.sequence)
