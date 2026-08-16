from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from ppigfinder.structure_prediction.hpc_planner import HPCResourcePlan


@dataclass(frozen=True)
class SlurmDirectiveOverrides:
    partition: Optional[str] = None
    gres: Optional[str] = None
    notes: List[str] = field(default_factory=list)


@dataclass(frozen=True)
class SlurmClusterProfile:
    name: str
    partition_by_hint: Dict[str, str]
    gres_by_gpu_hint: Dict[str, str]
    notes: List[str] = field(default_factory=list)

    def resolve(self, plan: HPCResourcePlan) -> SlurmDirectiveOverrides:
        partition = self.partition_by_hint.get(plan.partition_hint)
        gres = self.gres_by_gpu_hint.get(plan.gpu_hint)

        notes = list(self.notes)

        if partition is None:
            notes.append(
                f"No partition mapping for hint '{plan.partition_hint}'. Review manually."
            )

        if plan.gpu_hint != "none_by_default" and gres is None:
            notes.append(
                f"No GRES mapping for GPU hint '{plan.gpu_hint}'. Review manually."
            )

        return SlurmDirectiveOverrides(
            partition=partition,
            gres=gres,
            notes=notes,
        )



GENERIC_SLURM_PROFILE = SlurmClusterProfile(
    name="generic_slurm",
    partition_by_hint={},
    gres_by_gpu_hint={
        "none_by_default": "",
    },
    notes=[
        "Generic Slurm profile: review partition and GPU/GRES directives for the target cluster.",
        "DaVinci is an optional preconfigured institutional profile, not a ppigFinder requirement.",
    ],
)

DAVINCI_PROFILE = SlurmClusterProfile(
    name="davinci",
    partition_by_hint={
        "gpu_shared_or_standard": "max50",
        "gpu_standard": "max50",
        "gpu_high_memory_preferred": "max90",
        "gpu_high_memory_required": "unrestricted",
        "cpu_or_backend_specific": "basic",
        "review_required": "unrestricted",
    },
    gres_by_gpu_hint={
        "standard_gpu_or_slice": "shard:10",
        "one_standard_gpu": "shard:32",
        "one_high_memory_gpu": "shard:50",
        "largest_available_gpu_memory": "shard:50",
        "none_by_default": "",
    },
    notes=[
        "DaVinci profile maps generic ppigFinder hints to local Slurm partitions.",
        "Review GRES values against current sinfo/scontrol before production submission.",
    ],
)


PROFILES: Dict[str, SlurmClusterProfile] = {
    "generic_slurm": GENERIC_SLURM_PROFILE,
    "davinci": DAVINCI_PROFILE,
}


def get_cluster_profile(name: str = "generic_slurm") -> SlurmClusterProfile:
    key = name.lower()
    if key not in PROFILES:
        raise KeyError(f"Unknown cluster profile: {name}")
    return PROFILES[key]


def resolve_slurm_overrides(
    plan: HPCResourcePlan,
    cluster: str = "generic_slurm",
) -> SlurmDirectiveOverrides:
    return get_cluster_profile(cluster).resolve(plan)
