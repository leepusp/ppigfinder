#!/usr/bin/env python3
"""
UI shell models for ppigFinder.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class HomeAction:
    """
    One action card on the home/start screen.
    """

    id: str
    title: str
    description: str
    input_data: str
    output_data: str
    action_name: str | None = None


@dataclass
class WorkspaceSection:
    """
    One section in the future workflow workspace.
    """

    id: str
    title: str
    description: str
    data_type: str
