#!/usr/bin/env python3
"""
Temporary application launcher.

The modular package structure is already in place, but the active GUI class
is still provided by legacy_v20.py. During the refactor, features should be
moved from legacy_v20.py into dedicated modules.
"""

import sys

from .infrastructure.ipython_runtime import configure_ipython_qt_event_loop

from .legacy_v20 import (
    QApplication,
    QT_VERSION,
    _setup_emoji_font,
    _check_dependencies_at_startup,
    ppigFinderApp,
)

from .ui.icon_provider import set_window_icon
from .ui.text_fallback import apply_text_fallback_to_window
from .ui.window_manager import install_window_management
from .ui.toolbar import polish_toolbars
from .ui.recent_files import install_recent_files_menu
from .ui.file_opening import open_genome_file_into_window


def main() -> int:
    configure_ipython_qt_event_loop()

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

    set_window_icon(window)
    install_window_management(window)
    install_recent_files_menu(window, lambda path: open_genome_file_into_window(window, path))
    polish_toolbars(window)
    apply_text_fallback_to_window(window)

    window.show()

    return app.exec() if QT_VERSION == 6 else app.exec_()


if __name__ == "__main__":
    raise SystemExit(main())
