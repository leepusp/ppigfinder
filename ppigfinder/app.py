#!/usr/bin/env python3
"""
Temporary application launcher.

The modular package structure is already in place, but the active GUI class
is still provided by legacy_v20.py. During the refactor, features should be
moved from legacy_v20.py into dedicated modules.
"""

import sys

from .legacy_v20 import (
    QApplication,
    QT_VERSION,
    _setup_emoji_font,
    _check_dependencies_at_startup,
    ppigFinderApp,
)


def main() -> int:
    app = QApplication(sys.argv)
    app.setApplicationName("ppigFinder")
    app.setApplicationDisplayName(
        "ppigFinder — Protein-Protein Interaction Genomic Finder"
    )
    app.setApplicationVersion("2.0")
    app.setStyle("Fusion")

    _setup_emoji_font(app)
    _check_dependencies_at_startup()

    window = ppigFinderApp()
    window.show()

    return app.exec() if QT_VERSION == 6 else app.exec_()


if __name__ == "__main__":
    raise SystemExit(main())
