"""
Core domain models for ppigFinder.
"""

from .genome import GenomeRecord
from .orf import ORF, ORFSet, DomainHit
from .annotation import (
    domain_hit_from_legacy_dict,
    orf_from_legacy_dict,
    orf_to_legacy_dict,
)

__all__ = [
    "GenomeRecord",
    "ORF",
    "ORFSet",
    "DomainHit",
    "domain_hit_from_legacy_dict",
    "orf_from_legacy_dict",
    "orf_to_legacy_dict",
]
