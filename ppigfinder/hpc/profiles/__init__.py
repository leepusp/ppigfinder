#!/usr/bin/env python3
"""HPC profile registry exports."""

from __future__ import annotations

from ppigfinder.hpc.profiles.generic import GENERIC_AF3, GenericAF3Profile
from ppigfinder.hpc.profiles.davinci import DAVINCI_AF3, DaVinciAF3Profile

__all__ = [
    "GENERIC_AF3",
    "GenericAF3Profile",
    "DAVINCI_AF3",
    "DaVinciAF3Profile",
]
