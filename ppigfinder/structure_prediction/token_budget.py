from __future__ import annotations

from typing import Iterable, List

from ppigfinder.structure_prediction.models import PredictionJobSpec, SequenceTarget


def estimate_job_tokens(job: PredictionJobSpec) -> int:
    return job.estimated_tokens()


def classify_token_load(token_count: int) -> str:
    if token_count <= 1000:
        return "small"
    if token_count <= 2500:
        return "medium"
    if token_count <= 5000:
        return "large"
    return "xlarge"


def partition_targets_by_token_budget(
    targets: Iterable[SequenceTarget],
    max_tokens: int = 5000,
) -> List[List[SequenceTarget]]:
    batches: List[List[SequenceTarget]] = []
    current: List[SequenceTarget] = []
    current_tokens = 0

    for target in targets:
        tokens = target.token_length()

        if current and current_tokens + tokens > max_tokens:
            batches.append(current)
            current = []
            current_tokens = 0

        current.append(target)
        current_tokens += tokens

    if current:
        batches.append(current)

    return batches
