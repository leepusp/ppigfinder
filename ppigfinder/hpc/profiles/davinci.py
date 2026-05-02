#!/usr/bin/env python3
"""
DaVinci cluster profile.

The current DaVinci AlphaFold 3 integration assumes that the `af3` command
is already available after login. No `module load alphafold3` is required
by default.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class DaVinciAF3Profile:
    """
    Default AF3 submission options for DaVinci.
    """

    command: str = "af3"
    requires_module_load: bool = False
    module_command: str | None = None

    default_partition: str = "max50"
    default_nodes: int = 1
    default_ntasks: int = 16
    default_mem: str = "128G"
    default_time: str = "1:00:00"

    # For shared GPU scheduling, adjust if the wrapper expects another syntax.
    # Common alternatives are "shard:10" or "gres/shard=10".
    default_gres: str = "shard:10"

    default_resource_mode: str = "shared"


DAVINCI_AF3 = DaVinciAF3Profile()
