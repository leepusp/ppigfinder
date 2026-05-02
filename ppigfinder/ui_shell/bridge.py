#!/usr/bin/env python3
"""
Bridge between the future UI shell and the current legacy GUI/backend.

This allows new screens to call existing actions without directly depending
on legacy_v20.py internals.
"""

from __future__ import annotations


class LegacyActionBridge:
    """
    Safely call actions on a legacy ppigFinder window.
    """

    def __init__(self, legacy_window=None):
        self.legacy_window = legacy_window

    def available(self, action_name: str | None) -> bool:
        if not action_name or self.legacy_window is None:
            return False

        return callable(getattr(self.legacy_window, action_name, None))

    def call(self, action_name: str | None) -> bool:
        if not self.available(action_name):
            return False

        method = getattr(self.legacy_window, action_name)

        method()
        return True


class PreviewActionBridge:
    """
    Bridge used by standalone preview mode.
    """

    def call(self, action_name: str | None) -> bool:
        return False

    def available(self, action_name: str | None) -> bool:
        return False
