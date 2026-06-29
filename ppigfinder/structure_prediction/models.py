from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass(frozen=True)
class StructuralBackend:
    identifier: str
    display_name: str
    supports_monomer: bool
    supports_complex: bool
    supports_multicomponent: bool
    supports_gpu: bool
    notes: str = ""


AF3_BACKEND = StructuralBackend(
    identifier="af3",
    display_name="AlphaFold 3",
    supports_monomer=True,
    supports_complex=True,
    supports_multicomponent=True,
    supports_gpu=True,
    notes="AF3 backend with token-aware JSON generation and HPC scheduling.",
)

BOLTZ2_BACKEND = StructuralBackend(
    identifier="boltz2",
    display_name="Boltz-2",
    supports_monomer=True,
    supports_complex=True,
    supports_multicomponent=True,
    supports_gpu=True,
    notes="Alternative structural prediction backend for model comparison.",
)

FOLDCP_BACKEND = StructuralBackend(
    identifier="foldcp",
    display_name="FoldCP",
    supports_monomer=True,
    supports_complex=True,
    supports_multicomponent=True,
    supports_gpu=False,
    notes="Optional structural comparison or prediction-support backend.",
)

BACKENDS: Dict[str, StructuralBackend] = {
    "af3": AF3_BACKEND,
    "boltz2": BOLTZ2_BACKEND,
    "foldcp": FOLDCP_BACKEND,
}


@dataclass
class SequenceTarget:
    target_id: str
    sequence: str
    molecule_type: str = "protein"
    chain_id: Optional[str] = None
    role: str = "target"
    metadata: Dict[str, str] = field(default_factory=dict)

    def token_length(self) -> int:
        return len((self.sequence or "").strip())


@dataclass
class PredictionJobSpec:
    job_id: str
    backend_id: str
    targets: List[SequenceTarget]
    model_mode: str = "complex"
    priority: str = "normal"
    metadata: Dict[str, str] = field(default_factory=dict)

    def estimated_tokens(self) -> int:
        return sum(target.token_length() for target in self.targets)

    def target_count(self) -> int:
        return len(self.targets)

    def is_multicomponent(self) -> bool:
        return len(self.targets) > 2
