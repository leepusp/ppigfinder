#!/usr/bin/env python3
"""
ppigFinder application launcher.

This module is intentionally lightweight at import time. Heavy GUI imports
are loaded only when main() is executed.
"""

from __future__ import annotations

import os as _ppig_os
import time as _ppig_time
from pathlib import Path as _ppig_Path


def _ppig_startup_profiler_enabled() -> bool:
    return _ppig_os.environ.get("PPIG_PROFILE_STARTUP") in {"1", "true", "TRUE", "yes", "YES"}


def _ppig_startup_log_path() -> _ppig_Path:
    return _ppig_Path(
        _ppig_os.environ.get(
            "PPIG_APP_STARTUP_LOG",
            "docs/developer/startup_profile/app_startup_phases.tsv",
        )
    )


def _ppig_startup_log(event: str, elapsed_ms: float = 0.0, note: str = "") -> None:
    if not _ppig_startup_profiler_enabled():
        return

    path = _ppig_startup_log_path()
    path.parent.mkdir(parents=True, exist_ok=True)

    new_file = not path.exists()

    with path.open("a", encoding="utf-8") as handle:
        if new_file:
            handle.write("event\telapsed_ms\tnote\n")
        handle.write(f"{event}\t{elapsed_ms:.3f}\t{note}\n")
        handle.flush()


class _PpigStartupTimer:
    def __init__(self, event: str, note: str = ""):
        self.event = event
        self.note = note
        self.t0 = 0.0

    def __enter__(self):
        self.t0 = _ppig_time.perf_counter()
        return self

    def __exit__(self, exc_type, exc, tb):
        elapsed = (_ppig_time.perf_counter() - self.t0) * 1000.0
        suffix = self.note
        if exc is not None:
            suffix = f"{suffix} ERROR={exc}" if suffix else f"ERROR={exc}"
        _ppig_startup_log(self.event, elapsed, suffix)
        return False



import sys



def _ppig_signal_window_ready() -> None:
    """Signal the external splash launcher that the main window is visible."""
    ready_file = _ppig_os.environ.get("PPIG_READY_FILE")
    if not ready_file:
        return

    try:
        path = _ppig_Path(ready_file)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("ready\n", encoding="utf-8")
    except Exception:
        pass


def main() -> int:
    with _PpigStartupTimer("import.ipython_runtime"):
        from ppigfinder.infrastructure.ipython_runtime import configure_ipython_qt_event_loop

    with _PpigStartupTimer("configure_ipython_qt_event_loop"):
        configure_ipython_qt_event_loop()

    with _PpigStartupTimer("import.legacy_gui"):
        from ppigfinder.legacy_v29_11_7 import (
            QApplication,
            QT_VERSION,
            QTimer,
            _setup_emoji_font,
            _check_dependencies_at_startup,
            ppigFinderApp,
        )

    with _PpigStartupTimer("import.ui_helpers"):
        from ppigfinder.ui.icon_provider import set_window_icon
        from ppigfinder.ui.text_fallback import apply_text_fallback_to_window
        from ppigfinder.ui.window_manager import install_window_management
        from ppigfinder.ui.toolbar import polish_toolbars
        from ppigfinder.ui.recent_files import install_recent_files_menu
        from ppigfinder.ui.file_opening import open_genome_file_into_window
        from ppigfinder.ui.menu_installers import install_modular_gui_actions
        from ppigfinder.ui.tab_compactor import compact_tab_labels

    with _PpigStartupTimer("QApplication"):
        app = QApplication(sys.argv)
    app.setApplicationName("ppigFinder")
    app.setApplicationDisplayName(
        "ppigFinder — Protein-Protein Interaction Genomic Finder"
    )
    app.setApplicationVersion("2.0")
    app.setStyle("Fusion")

    _setup_emoji_font(app)

    with _PpigStartupTimer("ppigFinderApp.__init__"):
        window = ppigFinderApp()

    with _PpigStartupTimer("set_window_icon"):
        set_window_icon(window)
    with _PpigStartupTimer("install_window_management"):
        install_window_management(window, add_menu=False)

    with _PpigStartupTimer("window.show"):
        window.show()

    # Let Qt process the first paint/layout pass, then notify the splash.
    QTimer.singleShot(100, _ppig_signal_window_ready)

    def post_startup_recent_files():
        with _PpigStartupTimer("post_startup.install_recent_files_menu"):
            install_recent_files_menu(
                window,
                lambda path: open_genome_file_into_window(window, path),
            )

    def post_startup_toolbar_polish():
        with _PpigStartupTimer("post_startup.polish_toolbars"):
            polish_toolbars(window)

    def post_startup_modular_actions():
        with _PpigStartupTimer("post_startup.install_modular_gui_actions"):
            install_modular_gui_actions(window)

    def post_startup_compact_tabs():
        with _PpigStartupTimer("post_startup.compact_tab_labels"):
            compact_tab_labels(window)

    def post_startup_text_fallback():
        with _PpigStartupTimer("post_startup.apply_text_fallback_to_window"):
            apply_text_fallback_to_window(window)

    def post_startup_dependency_check():
        with _PpigStartupTimer("post_startup._check_dependencies_at_startup"):
            _check_dependencies_at_startup()

    # Keep first paint responsive: run non-critical UI polish after the window
    # has already entered the Qt event loop. Backend checks are intentionally
    # delayed because they may touch optional dependencies or filesystem paths.
    QTimer.singleShot(350, post_startup_recent_files)
    QTimer.singleShot(700, post_startup_toolbar_polish)
    QTimer.singleShot(1000, post_startup_modular_actions)
    QTimer.singleShot(1300, post_startup_compact_tabs)
    QTimer.singleShot(1700, post_startup_text_fallback)
    QTimer.singleShot(2600, post_startup_dependency_check)

    _ppig_startup_log("app.exec.enter", 0.0, "entering Qt event loop")
    return app.exec() if QT_VERSION == 6 else app.exec_()


if __name__ == "__main__":
    raise SystemExit(main())
