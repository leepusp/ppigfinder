#!/usr/bin/env python3
"""
Central installer for modular GUI actions.

This module adds new actions to the legacy GUI without editing legacy_v20.py.
It keeps the refactor safer by installing extensions after the main window
has been created.
"""

from __future__ import annotations

try:
    from PyQt6.QtGui import QAction
except Exception:
    from PyQt5.QtWidgets import QAction

from ppigfinder.ui.text_fallback import clean_ui_text


def _normalized_menu_title(text: str) -> str:
    return clean_ui_text(text).replace("&", "").strip().lower()


def get_or_create_menu(window, title: str):
    """
    Return an existing menu by title or create it.
    """
    menu_bar = window.menuBar()
    wanted = _normalized_menu_title(title)

    for action in menu_bar.actions():
        current = _normalized_menu_title(action.text())
        if current == wanted and action.menu() is not None:
            return action.menu()

    return menu_bar.addMenu(title)


def install_action_once(
    window,
    menu_title: str,
    action_text: str,
    attribute_name: str,
    callback,
    separator_before: bool = False,
) -> None:
    """
    Install a QAction only once on a menu.
    """
    if getattr(window, attribute_name, False):
        return

    menu = get_or_create_menu(window, menu_title)

    if separator_before:
        menu.addSeparator()

    action = QAction(clean_ui_text(action_text), window)
    action.triggered.connect(lambda checked=False: callback(window))
    menu.addAction(action)

    setattr(window, attribute_name, True)


def install_modular_gui_actions(window) -> None:
    """
    Install ppigFinder refactor actions into the legacy GUI.
    """
    # Imports are intentionally local to keep startup light and avoid breaking
    # the GUI if one optional extension has a problem.

    try:
        from ppigfinder.ui.backend_status import install_backend_status_action

        install_backend_status_action(window)
    except Exception:
        pass

    try:
        from ppigfinder.ui.af3_export import export_selected_orfs_as_server_json

        install_action_once(
            window,
            menu_title="AlphaFold",
            action_text="Export AF3 Server JSON",
            attribute_name="_ppig_action_export_af3_server_json_installed",
            callback=export_selected_orfs_as_server_json,
            separator_before=False,
        )
    except Exception:
        pass

    try:
        from ppigfinder.ui.af3_results_import import import_af3_results_folder

        install_action_once(
            window,
            menu_title="AlphaFold",
            action_text="Import AF3 Results Folder",
            attribute_name="_ppig_action_import_af3_results_folder_installed",
            callback=import_af3_results_folder,
            separator_before=False,
        )
    except Exception:
        pass

    try:
        from ppigfinder.ui.af3_table_export import export_af3_results_table

        install_action_once(
            window,
            menu_title="AlphaFold",
            action_text="Export AF3 Results Table",
            attribute_name="_ppig_action_export_af3_results_table_installed",
            callback=export_af3_results_table,
            separator_before=False,
        )
    except Exception:
        pass

    try:
        from ppigfinder.ui.project_bridge import (
            export_project_snapshot_from_window,
            import_project_snapshot_into_window,
        )

        install_action_once(
            window,
            menu_title="Project",
            action_text="Export Project Snapshot v3",
            attribute_name="_ppig_action_export_project_snapshot_installed",
            callback=export_project_snapshot_from_window,
            separator_before=False,
        )

        install_action_once(
            window,
            menu_title="Project",
            action_text="Import Project Snapshot v3",
            attribute_name="_ppig_action_import_project_snapshot_installed",
            callback=import_project_snapshot_into_window,
            separator_before=False,
        )
    except Exception:
        pass
