#!/usr/bin/env python3
"""
Generic Qt worker infrastructure for long-running tasks.

This is used to move heavy backend operations away from the GUI thread.
"""

from __future__ import annotations

import traceback


try:
    from PyQt6.QtCore import QObject, QThread, pyqtSignal
except Exception:
    from PyQt5.QtCore import QObject, QThread, pyqtSignal


class WorkerSignals(QObject):
    """
    Signals emitted by a background worker.
    """

    started = pyqtSignal()
    finished = pyqtSignal(object)
    failed = pyqtSignal(str)
    progress = pyqtSignal(str)


class FunctionWorker(QObject):
    """
    Run a Python function in a QThread.
    """

    def __init__(self, function, *args, **kwargs):
        super().__init__()
        self.function = function
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

    def run(self):
        self.signals.started.emit()

        try:
            result = self.function(*self.args, **self.kwargs)
        except Exception:
            self.signals.failed.emit(traceback.format_exc())
            return

        self.signals.finished.emit(result)


def run_in_thread(
    parent,
    function,
    on_finished=None,
    on_failed=None,
    on_started=None,
    *args,
    **kwargs,
):
    """
    Run a function in a background QThread.

    The thread/worker references are stored on parent._active_workers to avoid
    premature garbage collection.
    """
    thread = QThread(parent)
    worker = FunctionWorker(function, *args, **kwargs)
    worker.moveToThread(thread)

    if not hasattr(parent, "_active_workers"):
        parent._active_workers = []

    parent._active_workers.append((thread, worker))

    thread.started.connect(worker.run)

    if on_started:
        worker.signals.started.connect(on_started)

    if on_finished:
        worker.signals.finished.connect(on_finished)

    if on_failed:
        worker.signals.failed.connect(on_failed)

    def _cleanup():
        try:
            parent._active_workers.remove((thread, worker))
        except Exception:
            pass

        worker.deleteLater()
        thread.quit()
        thread.wait()
        thread.deleteLater()

    worker.signals.finished.connect(lambda _: _cleanup())
    worker.signals.failed.connect(lambda _: _cleanup())

    thread.start()
    return thread, worker
