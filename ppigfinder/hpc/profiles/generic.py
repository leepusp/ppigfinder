#!/usr/bin/env python3
"""
Generic HPC / AF3 profile.

This profile intentionally avoids site-specific Slurm assumptions. It is the
safe default for ppigFinder because users may export AlphaFold 3 JSON files,
run local wrappers, or adapt generated scripts to any SSH/Slurm/HPC server.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class GenericAF3Profile:
    """
    Generic AF3 execution defaults.

    Site-specific values such as partition, GRES, memory and wall time should
    be supplied by the user or by a named cluster profile.
    """

    command: str = "af3"
    requires_module_load: bool = False
    module_command: str | None = None

    default_partition: str | None = None
    default_nodes: int | None = None
    default_ntasks: int | None = None
    default_mem: str | None = None
    default_time: str | None = None
    default_gres: str | None = None
    default_resource_mode: str | None = None


GENERIC_AF3 = GenericAF3Profile()
