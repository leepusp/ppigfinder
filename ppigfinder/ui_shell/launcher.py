#!/usr/bin/env python3
"""
Experimental ppigFinder shell launcher.

This launcher does not replace main.py yet. It provides a future entry point:
Splash -> Home -> current legacy interface.
"""

from __future__ import annotations

import sys

from ppigfinder.ui_shell.qt import QApplication, QTimer, exec_app
from ppigfinder.ui_shell.splash import SplashWindow
from ppigfinder.ui_shell.home_window import HomeWindow
from ppigfinder.ui_shell.bridge import LegacyActionBridge
from ppigfinder.ui_shell.theme import shell_stylesheet


class ShellController:
    """
    Controls the experimental startup flow.
    """

    def __init__(self):
        self.legacy_window = None
        self.home = None
        self.splash = None

    def start(self) -> None:
        self.splash = SplashWindow()
        self.splash.show()
        self.splash.start()

        QTimer.singleShot(1700, self.show_home)

    def show_home(self) -> None:
        if self.splash is not None:
            self.splash.close()

        self.home = HomeWindow()
        self._add_open_legacy_button()
        self.home.show()

    def _add_open_legacy_button(self) -> None:
        """
        Add an explicit action to open the current working ppigFinder interface.
        """
        from ppigfinder.ui_shell.models import HomeAction

        action = HomeAction(
            id="open_current_interface",
            title="Open current ppigFinder interface",
            description="Open the current full ppigFinder interface while the guided interface is under development.",
            input_data="Current application",
            output_data="Current ppigFinder workspace",
            action_name="open_legacy_interface",
        )

        self.home.actions.insert(0, action)

        # Rebuild home after inserting the action.
        self.home.close()
        self.home = HomeWindow(
            bridge=ShellBridge(self),
            actions=self.home.actions,
        )

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

        QTimer.singleShot(100, _check_dependencies_at_startup)

        return True


class ShellBridge(LegacyActionBridge):
    """
    Bridge used by the experimental shell.
    """

    def __init__(self, controller: ShellController):
        super().__init__(None)
        self.controller = controller

    def available(self, action_name: str | None) -> bool:
        if action_name == "open_legacy_interface":
            return True

        if self.controller.legacy_window is None:
            return False

        return callable(getattr(self.controller.legacy_window, action_name, None))

    def call(self, action_name: str | None) -> bool:
        if action_name == "open_legacy_interface":
            return self.controller.open_legacy_interface()

        if self.controller.legacy_window is None:
            self.controller.open_legacy_interface()

        if self.controller.legacy_window is None or action_name is None:
            return False

        method = getattr(self.controller.legacy_window, action_name, None)

        if callable(method):
            method()
            return True

        return False


def main() -> int:
    app = QApplication(sys.argv)
    app.setApplicationName("ppigFinder Experimental Shell")
    app.setStyleSheet(shell_stylesheet())

    controller = ShellController()
    controller.start()

    return exec_app(app)


if __name__ == "__main__":
    raise SystemExit(main())
