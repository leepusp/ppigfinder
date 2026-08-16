"""
Infrastructure modules for ppigFinder.
"""

from .cache import JsonCache, directory_signature, file_signature
from .ipython_runtime import (
    configure_ipython_qt_event_loop,
    get_ipython_shell,
    is_running_in_ipython,
)

__all__ = [
    "JsonCache",
    "directory_signature",
    "file_signature",
    "configure_ipython_qt_event_loop",
    "get_ipython_shell",
    "is_running_in_ipython",
]

from .parallel import (
    recommended_workers,
    chunked,
    parallel_map,
)

__all__ += [
    "recommended_workers",
    "chunked",
    "parallel_map",
]
