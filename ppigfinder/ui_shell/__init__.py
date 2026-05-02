"""
Future UI shell for ppigFinder.

This package is experimental and does not replace the current legacy GUI yet.
"""

from .home_window import HomeWindow
from .workspace_window import WorkspaceWindow
from .splash import SplashWindow

__all__ = [
    "HomeWindow",
    "WorkspaceWindow",
    "SplashWindow",
]
