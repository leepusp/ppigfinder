#!/usr/bin/env python3
"""
Experimental ppigFinder shell launcher.

This launcher does not replace main.py yet. It provides a future entry point:

Splash
  -> Home
      -> current legacy interface
      -> guided workspace preview
      -> direct calls into the current interface
"""

from __future__ import annotations

import sys

from ppigfinder.ui_shell.qt import QApplication, QTimer, exec_app
from ppigfinder.ui_shell.splash import SplashWindow
from ppigfinder.ui_shell.home_window import HomeWindow
from ppigfinder.ui_shell.workspace_window import WorkspaceWindow
from ppigfinder.ui_shell.bridge import LegacyActionBridge
from ppigfinder.ui_shell.models import HomeAction
from ppigfinder.ui_shell.theme import shell_stylesheet


class ShellController:
    """
    Controls the experimental startup flow.
    """

    def __init__(self):
        self.legacy_window = None
        self.home = None
        self.splash = None
        self.workspace = None

    def start(self) -> None:
        self.splash = SplashWindow()
        self.splash.show()
        self.splash.start()

        QTimer.singleShot(1600, self.show_home)

    def show_home(self) -> None:
        if self.splash is not None:
            self.splash.close()

        self.home = HomeWindow(
            bridge=ShellBridge(self),
            actions=self._home_actions(),
        )
        self.home.show()

    def _home_actions(self) -> list[HomeAction]:
        """
        Home actions for the experimental shell.

        action_name values are resolved by ShellBridge. Existing legacy
        actions are opened on demand.
        """
        return [
            HomeAction(
                id="open_genome",
                title="Start new analysis",
                description="Begin by adding genome data or selecting an input file.",
                input_data="DNA / genome sequence",
                output_data="Genome workspace",
                action_name="open_workspace:data",
            ),
            HomeAction(
                id="open_project",
                title="Open project",
                description="Resume a previous ppigFinder session.",
                input_data="Project file",
                output_data="Restored session",
                action_name="open_workspace:data",
            ),
            HomeAction(
                id="predict_orfs",
                title="Predict ORFs",
                description="Identify protein-coding regions in the loaded genome.",
                input_data="Loaded DNA / genome",
                output_data="ORF and protein table",
                action_name="open_workspace:orfs",
            ),
            HomeAction(
                id="guided_workspace",
                title="Guided workspace preview",
                description="Open the future workflow-oriented workspace preview.",
                input_data="Experimental interface shell",
                output_data="Stepwise analysis workspace",
                action_name="open_workspace:overview",
            ),
            HomeAction(
                id="reports",
                title="Reports",
                description="Generate HTML reports, snapshots and tabular exports.",
                input_data="Current project state",
                output_data="HTML / JSON / TSV",
                action_name="open_workspace:reports",
            ),
        ]

    def open_legacy_interface(self) -> bool:
        """
        Open the current legacy ppigFinder interface.
        """
        if self.legacy_window is not None:
            self.legacy_window.show()
            self.legacy_window.raise_()
            return True

        from ppigfinder.legacy_v20 import (
            _setup_emoji_font,
            _check_dependencies_at_startup,
            ppigFinderApp,
        )
        from ppigfinder.ui.icon_provider import set_window_icon
        from ppigfinder.ui.window_manager import install_window_management
        from ppigfinder.ui.toolbar import polish_toolbars
        from ppigfinder.ui.text_fallback import apply_text_fallback_to_window
        from ppigfinder.ui.menu_installers import install_modular_gui_actions
        from ppigfinder.ui.tab_compactor import compact_tab_labels

        app = QApplication.instance()
        if app is not None:
            _setup_emoji_font(app)

        self.legacy_window = ppigFinderApp()

        set_window_icon(self.legacy_window)
        install_window_management(self.legacy_window, add_menu=False)
        polish_toolbars(self.legacy_window)
        install_modular_gui_actions(self.legacy_window)
        compact_tab_labels(self.legacy_window)
        apply_text_fallback_to_window(self.legacy_window)

        self.legacy_window.show()
        self.legacy_window.raise_()

        QTimer.singleShot(100, _check_dependencies_at_startup)

        return True

    def open_guided_workspace(self, route_id: str = "overview") -> bool:
        """
        Open the future guided workspace preview and optionally switch route.
        """
        if self.workspace is None:
            self.workspace = WorkspaceWindow(bridge=ShellBridge(self))

        self.workspace.show()
        self.workspace.raise_()

        try:
            self.workspace.show_route(route_id)
        except Exception:
            pass

        return True

    def call_legacy_action(self, action_name: str | None) -> bool:
        """
        Open legacy interface if needed and call one of its existing methods.
        """
        if not action_name:
            return False

        if action_name == "open_legacy_interface":
            return self.open_legacy_interface()

        if action_name == "open_guided_workspace":
            return self.open_guided_workspace()

        if action_name.startswith("open_workspace:"):
            route_id = action_name.split(":", 1)[1] or "overview"
            return self.open_guided_workspace(route_id)

        self.open_legacy_interface()

        if self.legacy_window is None:
            return False

        method = getattr(self.legacy_window, action_name, None)

        if not callable(method):
            try:
                self.legacy_window._status.showMessage(
                    f"Action not available yet: {action_name}"
                )
            except Exception:
                pass
            return False

        method()
        return True


class ShellBridge(LegacyActionBridge):
    """
    Bridge used by the experimental shell.

    It marks actions as available because it can open the legacy interface
    on demand before calling them.
    """

    def __init__(self, controller: ShellController):
        super().__init__(None)
        self.controller = controller

    def available(self, action_name: str | None) -> bool:
        return action_name is not None

    def call(self, action_name: str | None) -> bool:
        return self.controller.call_legacy_action(action_name)


def main() -> int:
    app = QApplication(sys.argv)
    app.setApplicationName("ppigFinder Experimental Shell")
    app.setStyleSheet(shell_stylesheet())

    controller = ShellController()
    controller.start()

    return exec_app(app)


if __name__ == "__main__":
    raise SystemExit(main())
