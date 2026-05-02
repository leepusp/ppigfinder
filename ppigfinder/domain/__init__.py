"""
Core domain models for ppigFinder.
"""

from .genome import GenomeRecord
from .project import ProjectState, ProjectMetadata, PROJECT_SCHEMA_VERSION
from .blast import BlastHit, blast_hit_from_legacy_dict
from .hmm import HMMHit, hmm_hit_from_legacy_dict
from .orf import ORF, ORFSet, DomainHit
from .annotation import (
    domain_hit_from_legacy_dict,
    orf_from_legacy_dict,
    orf_to_legacy_dict,
)

__all__ = [
    "ProjectState",
    "ProjectMetadata",
    "PROJECT_SCHEMA_VERSION",
    "HMMHit",
    "hmm_hit_from_legacy_dict",
    "BlastHit",
    "blast_hit_from_legacy_dict",
    "GenomeRecord",
    "ORF",
    "ORFSet",
    "DomainHit",
    "domain_hit_from_legacy_dict",
    "orf_from_legacy_dict",
    "orf_to_legacy_dict",
]
