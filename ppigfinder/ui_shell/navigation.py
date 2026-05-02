#!/usr/bin/env python3
"""
Navigation models for the future ppigFinder workspace shell.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class ModuleRoute:
    id: str
    title: str
    description: str
    data_type: str
    status: str = "Not started"
