#!/usr/bin/env python3
"""
Optional IPython/Jupyter runtime helpers.

IPython is not required by ppigFinder. This module only improves behavior
when the application is launched from an IPython shell or notebook-like
environment.
"""

from __future__ import annotations


def get_ipython_shell():
    """
    Return the active IPython shell if available.
    """
    try:
        from IPython import get_ipython
    except Exception:
        return None

    try:
        return get_ipython()
    except Exception:
        return None


def is_running_in_ipython() -> bool:
    """
    True when ppigFinder is being launched from IPython/Jupyter.
    """
    return get_ipython_shell() is not None


def configure_ipython_qt_event_loop() -> bool:
    """
    Enable Qt event loop integration when running under IPython.

    Returns True if an IPython shell was detected and a Qt GUI loop was
    requested successfully.
    """
    shell = get_ipython_shell()

    if shell is None:
        return False

    enable_gui = getattr(shell, "enable_gui", None)

    if not callable(enable_gui):
        return False

    try:
        enable_gui("qt")
        return True
    except Exception:
        return False
