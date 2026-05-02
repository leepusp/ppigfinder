#!/usr/bin/env python3
"""
Standalone preview launcher for the future ppigFinder UI shell.
"""

from __future__ import annotations

import sys

from ppigfinder.ui_shell.home_window import HomeWindow
from ppigfinder.ui_shell.qt import QApplication, exec_app
from ppigfinder.ui_shell.theme import shell_stylesheet


def main() -> int:
    app = QApplication(sys.argv)
    app.setApplicationName("ppigFinder UI Shell Preview")
    app.setStyleSheet(shell_stylesheet())

    window = HomeWindow()
    window.show()

    return exec_app(app)


if __name__ == "__main__":
    raise SystemExit(main())
