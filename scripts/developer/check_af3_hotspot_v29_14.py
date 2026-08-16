#!/usr/bin/env python3
"""Self-check for the ppigFinder v29.14 AF3 PAE hotspot algorithm."""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import numpy as np

from ppigfinder.alphafold.hotspot import compute_pae_hotspot


def main() -> int:
    sub_ab = np.full((20, 25), 20.0, dtype=float)
    sub_ab[5:10, 8:14] = 3.0
    sub_ab[6, 9] = 2.0

    sub_ba = sub_ab.T + 0.5

    contact_ab = np.zeros_like(sub_ab)
    contact_ab[5:10, 8:14] = 0.8

    hotspot = compute_pae_hotspot(
        sub_AB=sub_ab,
        sub_BA=sub_ba,
        contact_AB=contact_ab,
        radius=10,
    )

    assert hotspot["min_row"] == 6, hotspot
    assert hotspot["min_col"] == 9, hotspot
    assert hotspot["hotspot_min_pae"] == 2.0, hotspot
    assert hotspot["hotspot_core_cells"] >= 25, hotspot
    assert hotspot["hotspot_score"] > 0.0, hotspot
    assert hotspot["hotspot_mean_cp"] is not None, hotspot

    print("OK: AF3 hotspot v29.14 self-check passed.")
    print("hotspot:", hotspot)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
