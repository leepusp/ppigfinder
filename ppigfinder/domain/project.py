#!/usr/bin/env python3
"""
Versioned project state models for ppigFinder.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from datetime import datetime
from typing import Any


PROJECT_SCHEMA_VERSION = 3


@dataclass(slots=True)
class ProjectMetadata:
    app: str = "ppigFinder"
    schema_version: int = PROJECT_SCHEMA_VERSION
    created_at: str = field(default_factory=lambda: datetime.now().isoformat(timespec="seconds"))
    updated_at: str = field(default_factory=lambda: datetime.now().isoformat(timespec="seconds"))


@dataclass(slots=True)
class ProjectState:
    metadata: ProjectMetadata = field(default_factory=ProjectMetadata)

    genome: dict[str, Any] = field(default_factory=dict)
    orfs: list[dict[str, Any]] = field(default_factory=list)

    blast: dict[str, Any] = field(default_factory=dict)
    hmm: dict[str, Any] = field(default_factory=dict)
    alphafold: dict[str, Any] = field(default_factory=dict)

    hpc: dict[str, Any] = field(default_factory=dict)
    exports: dict[str, Any] = field(default_factory=dict)
    ui: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        self.metadata.updated_at = datetime.now().isoformat(timespec="seconds")
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ProjectState":
        metadata_data = data.get("metadata", {}) or {}

        metadata = ProjectMetadata(
            app=metadata_data.get("app", "ppigFinder"),
            schema_version=int(metadata_data.get("schema_version", 1)),
            created_at=metadata_data.get("created_at", datetime.now().isoformat(timespec="seconds")),
            updated_at=metadata_data.get("updated_at", datetime.now().isoformat(timespec="seconds")),
        )

        return cls(
            metadata=metadata,
            genome=data.get("genome", {}) or {},
            orfs=data.get("orfs", []) or [],
            blast=data.get("blast", {}) or {},
            hmm=data.get("hmm", {}) or {},
            alphafold=data.get("alphafold", {}) or {},
            hpc=data.get("hpc", {}) or {},
            exports=data.get("exports", {}) or {},
            ui=data.get("ui", {}) or {},
        )
