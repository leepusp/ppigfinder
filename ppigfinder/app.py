#!/usr/bin/env python3
"""
ppigFinder application launcher.

This module is intentionally lightweight at import time. Heavy GUI imports
are loaded only when main() is executed.
"""

from __future__ import annotations

import sys


def main() -> int:
    from ppigfinder.infrastructure.ipython_runtime import configure_ipython_qt_event_loop

    configure_ipython_qt_event_loop()

    from ppigfinder.legacy_v20 import (
        QApplication,
        QT_VERSION,
        QTimer,
        _setup_emoji_font,
        _check_dependencies_at_startup,
        ppigFinderApp,
    )

    from ppigfinder.ui.icon_provider import set_window_icon
    from ppigfinder.ui.text_fallback import apply_text_fallback_to_window
    from ppigfinder.ui.window_manager import install_window_management
    from ppigfinder.ui.toolbar import polish_toolbars
    from ppigfinder.ui.recent_files import install_recent_files_menu
    from ppigfinder.ui.file_opening import open_genome_file_into_window
    from ppigfinder.ui.menu_installers import install_modular_gui_actions
    from ppigfinder.ui.tab_compactor import compact_tab_labels

    app = QApplication(sys.argv)
    app.setApplicationName("ppigFinder")
    app.setApplicationDisplayName(
        "ppigFinder — Protein-Protein Interaction Genomic Finder"
    )
    app.setApplicationVersion("2.0")
    app.setStyle("Fusion")

    _setup_emoji_font(app)

    window = ppigFinderApp()

    set_window_icon(window)
    install_window_management(window, add_menu=False)

    window.show()

    def post_startup_ui_polish():
        _check_dependencies_at_startup()
        install_recent_files_menu(
            window,
            lambda path: open_genome_file_into_window(window, path),
        )
        polish_toolbars(window)
        install_modular_gui_actions(window)
        compact_tab_labels(window)
        apply_text_fallback_to_window(window)

    QTimer.singleShot(100, post_startup_ui_polish)

    return app.exec() if QT_VERSION == 6 else app.exec_()


if __name__ == "__main__":
    raise SystemExit(main())
