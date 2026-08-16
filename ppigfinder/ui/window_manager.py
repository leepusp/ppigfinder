#!/usr/bin/env python3
"""
Responsive window management for ppigFinder.

This module avoids fixed-size windows and makes the desktop UI behave better
across local desktops, X11 forwarding, VNC and different monitor sizes.
"""

from __future__ import annotations

try:
    from PyQt6.QtCore import QSettings
    from PyQt6.QtGui import QAction
    from PyQt6.QtWidgets import QApplication, QInputDialog
except Exception:
    from PyQt5.QtCore import QSettings
    from PyQt5.QtWidgets import QApplication, QAction, QInputDialog

from ppigfinder.ui.text_fallback import clean_ui_text


WINDOW_PRESETS = {
    "Small 60%": (0.60, 0.65),
    "Medium 75%": (0.75, 0.80),
    "Large 85%": (0.85, 0.90),
    "X-Large 95%": (0.95, 0.95),
}


def _primary_available_geometry(window):
    try:
        screen = window.screen() if hasattr(window, "screen") else None
    except Exception:
        screen = None

    if screen is None:
        screen = QApplication.primaryScreen()

    if screen is None:
        return None

    return screen.availableGeometry()


def resize_relative(window, width_fraction: float = 0.85, height_fraction: float = 0.90) -> None:
    """
    Resize and center a window relative to the available screen.
    """
    rect = _primary_available_geometry(window)

    if rect is None:
        window.resize(1280, 800)
        return

    width = max(900, int(rect.width() * width_fraction))
    height = max(650, int(rect.height() * height_fraction))

    window.resize(width, height)

    frame = window.frameGeometry()
    frame.moveCenter(rect.center())
    window.move(frame.topLeft())


def _get_or_create_menu(menu_bar, title: str):
    """
    Return an existing menu by title, or create it if it does not exist.
    """
    wanted = clean_ui_text(title).replace("&", "")

    for action in menu_bar.actions():
        try:
            current = clean_ui_text(action.text()).replace("&", "")
            if current == wanted and action.menu() is not None:
                return action.menu()
        except Exception:
            pass

    return menu_bar.addMenu(title)


def add_window_size_menu(window) -> None:
    """
    Add responsive size presets to the existing Window menu.
    """
    if getattr(window, "_ppig_window_size_menu_installed", False):
        return

    try:
        menu_bar = window.menuBar()
    except Exception:
        return

    window_menu = _get_or_create_menu(menu_bar, "Window")
    window_menu.addSeparator()

    for label, fractions in WINDOW_PRESETS.items():
        action = QAction(label, window)
        action.triggered.connect(
            lambda checked=False, f=fractions: resize_relative(window, f[0], f[1])
        )
        window_menu.addAction(action)

    window_menu.addSeparator()

    maximize_action = QAction("Maximize", window)
    maximize_action.triggered.connect(window.showMaximized)
    window_menu.addAction(maximize_action)

    normal_action = QAction("Restore", window)
    normal_action.triggered.connect(window.showNormal)
    window_menu.addAction(normal_action)

    custom_action = QAction("Custom size...", window)

    def _custom_size():
        width, ok_w = QInputDialog.getInt(
            window,
            "Custom width",
            "Width:",
            max(900, window.width()),
            600,
            10000,
            50,
        )
        if not ok_w:
            return

        height, ok_h = QInputDialog.getInt(
            window,
            "Custom height",
            "Height:",
            max(650, window.height()),
            400,
            10000,
            50,
        )
        if not ok_h:
            return

        window.resize(width, height)

    custom_action.triggered.connect(_custom_size)
    window_menu.addAction(custom_action)

    window._ppig_window_size_menu_installed = True


def install_window_management(
    window,
    organization: str = "LEEPBioinfo",
    application: str = "ppigFinder",
    default_width_fraction: float = 0.85,
    default_height_fraction: float = 0.90,
    add_menu: bool = True,
) -> None:
    """
    Install responsive sizing and persistent geometry on a main window.
    """
    settings = QSettings(organization, application)

    try:
        window.setMinimumSize(900, 650)
    except Exception:
        pass

    restored = False

    try:
        geometry = settings.value("main_window/geometry")
        if geometry:
            restored = bool(window.restoreGeometry(geometry))
    except Exception:
        restored = False

    if not restored:
        resize_relative(window, default_width_fraction, default_height_fraction)

    try:
        state = settings.value("main_window/state")
        if state and hasattr(window, "restoreState"):
            window.restoreState(state)
    except Exception:
        pass

    if add_menu:
        add_window_size_menu(window)

    original_close_event = getattr(window, "closeEvent", None)

    def _close_event(event):
        try:
            settings.setValue("main_window/geometry", window.saveGeometry())
            if hasattr(window, "saveState"):
                settings.setValue("main_window/state", window.saveState())
        except Exception:
            pass

        if callable(original_close_event):
            original_close_event(event)
        else:
            event.accept()

    window.closeEvent = _close_event
