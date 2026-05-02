#!/usr/bin/env python3
"""
HMM/domain search service for ppigFinder.

This service separates HMM search orchestration from the GUI. The current
implementation delegates algorithm execution to the legacy analyzer while
providing a stable service API for the new architecture.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from ppigfinder.infrastructure.backends import BACKENDS


@dataclass(slots=True)
class HMMSearchParams:
    """
    Parameters for HMM/domain search.
    """

    evalue: float = 1e-3
    score_threshold: float = 0.0
    max_hits: int = 100
    use_hmmer: bool = True
    fallback_enabled: bool = True

    def to_legacy_dict(self) -> dict[str, Any]:
        return {
            "evalue": self.evalue,
            "score_threshold": self.score_threshold,
            "max_hits": self.max_hits,
            "use_hmmer": self.use_hmmer,
            "fallback_enabled": self.fallback_enabled,
        }


@dataclass(slots=True)
class HMMSearchResult:
    """
    Result returned by HMMSearchService.
    """

    orfs: list[dict] = field(default_factory=list)
    hits: list[dict] = field(default_factory=list)
    algorithm_used: str = ""
    backend_error: str = ""


class HMMSearchService:
    """
    Select and run HMMER/PSSM-like domain search backends.
    """

    def __init__(self, analyzer):
        self.analyzer = analyzer

    def search(
        self,
        orfs: list[dict],
        hmm_file: str | None = None,
        algorithm: str = "Auto",
        params: HMMSearchParams | None = None,
    ) -> HMMSearchResult:
        """
        Run HMM/domain search.

        This method is intentionally permissive because older ppigFinder
        versions used different analyzer method names. It tries known method
        names and falls back gracefully.
        """
        params = params or HMMSearchParams()
        backend_error = ""
        algorithm_used = algorithm
        hits: list[dict] = []

        self.analyzer._last_hmm_error = ""

        hmmer_available = BACKENDS.get("hmmer3", {}).get("available", False)

        should_try_hmmer = (
            algorithm.startswith("Auto")
            or algorithm.upper().startswith("HMMER")
            or algorithm.upper().startswith("HMM")
        ) and hmmer_available and hmm_file

        if should_try_hmmer:
            for method_name in [
                "run_hmmer_search",
                "run_hmm_search",
                "scan_hmm_profiles",
                "run_hmmer",
            ]:
                method = getattr(self.analyzer, method_name, None)
                if callable(method):
                    try:
                        result = method(orfs, hmm_file, params.to_legacy_dict())
                        algorithm_used = "HMMER3"
                        if isinstance(result, tuple):
                            orfs, hits = result
                        elif isinstance(result, list):
                            hits = result
                        return HMMSearchResult(
                            orfs=orfs,
                            hits=hits or [],
                            algorithm_used=algorithm_used,
                            backend_error="",
                        )
                    except Exception as exc:
                        backend_error = str(exc)

        # Fallback to built-in domain scanner if available.
        if params.fallback_enabled:
            for method_name in [
                "scan_domains",
                "scan_pssm_profiles",
                "run_pssm_scanner",
                "annotate_domains_builtin",
                "annotate_domains",
            ]:
                method = getattr(self.analyzer, method_name, None)
                if callable(method):
                    try:
                        result = method(orfs)
                        algorithm_used = "Built-in domain scanner"
                        if isinstance(result, tuple):
                            orfs, hits = result
                        elif isinstance(result, list):
                            # Some legacy methods return updated ORFs only.
                            orfs = result
                            hits = []
                        return HMMSearchResult(
                            orfs=orfs,
                            hits=hits or [],
                            algorithm_used=algorithm_used,
                            backend_error=backend_error,
                        )
                    except Exception as exc:
                        backend_error = str(exc)

        return HMMSearchResult(
            orfs=orfs,
            hits=hits,
            algorithm_used=algorithm_used,
            backend_error=backend_error or self.analyzer._last_hmm_error,
        )
