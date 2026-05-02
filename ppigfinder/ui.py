#!/usr/bin/env python3
"""
Convenience module for launching the experimental guided ppigFinder UI.

Usage:
    python -m ppigfinder.ui
"""

from __future__ import annotations

from ppigfinder.ui_shell.launcher import main


if __name__ == "__main__":
    raise SystemExit(main())
