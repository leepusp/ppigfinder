"""
AlphaFold/AF3 backend modules for ppigFinder.
"""

from .classifier import classify_metrics, classification_color
from .metrics import (
    InteractionMetrics,
    classify_interaction,
    contact_percent_from_pae,
    metrics_from_af3_summary,
    off_diagonal_blocks,
    pae_inter_from_blocks,
    pae_min_from_blocks,
)
from .result_parser import (
    AF3ParsedResult,
    parse_af3_result_directory,
    parse_af3_results_root,
)

__all__ = [
    "InteractionMetrics",
    "classify_interaction",
    "contact_percent_from_pae",
    "metrics_from_af3_summary",
    "off_diagonal_blocks",
    "pae_inter_from_blocks",
    "pae_min_from_blocks",
    "classify_metrics",
    "classification_color",
    "AF3ParsedResult",
    "parse_af3_result_directory",
    "parse_af3_results_root",
]
