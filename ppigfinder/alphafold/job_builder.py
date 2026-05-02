#!/usr/bin/env python3
"""
AlphaFold job builders.

This module exposes builders for AlphaFold Server JSON and, progressively,
will also expose local AlphaFold 3 JSON builders.
"""

from ppigfinder.alphafold.server_json import (
    ServerJob,
    ServerProteinChain,
    build_server_job,
    build_pair_job,
    build_pair_jobs_from_sequences,
    build_pair_jobs_from_legacy_orfs,
    jobs_to_json_payload,
    write_server_json,
    validate_server_payload,
)

__all__ = [
    "ServerJob",
    "ServerProteinChain",
    "build_server_job",
    "build_pair_job",
    "build_pair_jobs_from_sequences",
    "build_pair_jobs_from_legacy_orfs",
    "jobs_to_json_payload",
    "write_server_json",
    "validate_server_payload",
]
