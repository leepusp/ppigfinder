"""AlphaFold PAE hotspot metrics.

This module contains the ppigFinder v29.14 hotspot algorithm extracted from
the archived legacy GUI source so it can be tested and reused independently.
"""

from __future__ import annotations

def compute_pae_hotspot(sub_AB, sub_BA, contact_AB=None, radius: int = 10) -> dict:
    """v29.14 — anchor at the PAE minimum, then alternate grow/trim until the
    window is stable.

    Changes vs v29.11.x:
      1. The window may now CONTRACT as well as expand.  Previously it was
         anchored at argmin and could only grow, so whenever the PAE minimum
         sat at the rim of an interface (or on a plateau, where argmin returns
         the first index = top-left corner) the initial +/-radius box hung half
         outside the interface, the border means were polluted by high PAE, no
         expansion ever triggered, and the window froze at (2*radius+1)^2.
      2. Border tests are evaluated against the SAME window and applied
         together, so the result no longer depends on the order in which the
         four edges happen to be tested.
      3. recip_factor is 1.0 when reciprocity is UNKNOWN.  Previously a missing
         sub_BA produced the same 0.7 penalty as measurably poor reciprocity,
         silently docking 30% off the score for absent input.
    """
    if sub_AB is None or getattr(sub_AB, 'size', 0) == 0:
        return {}
    # v29.14 COHERENCE FIX: PAE_min is reported as the minimum over BOTH
    # directions, min(min(sub_AB), min(sub_BA)), but the window used to be
    # anchored only on argmin(sub_AB).  Whenever the true minimum lay in
    # the B->A block the hotspot was seeded at a different, higher-PAE
    # point, so hotspot_min_pae disagreed with the reported PAE_min and the
    # Seq A / Seq B anchor columns pointed at the wrong residues.  Anchor
    # on whichever direction actually holds the minimum, keeping chain A on
    # the rows, and use the other direction for the reciprocity term.
    if sub_BA is not None:
        try:
            if (sub_BA.shape == (sub_AB.shape[1], sub_AB.shape[0])
                    and float(sub_BA.min()) < float(sub_AB.min())):
                sub_AB, sub_BA = sub_BA.T, sub_AB.T
                contact_AB = None   # orientation no longer matches
        except Exception:
            pass
    nA, nB = sub_AB.shape
    if nA == 0 or nB == 0:
        return {}
    min_idx = int(sub_AB.argmin())
    min_r = min_idx // nB
    min_c = min_idx % nB
    min_val = float(sub_AB[min_r, min_c])
    soft_thresh = min(min_val * 2.0 + 2.0, 12.0)

    r0 = max(0, min_r - radius); r1 = min(nA, min_r + radius + 1)
    c0 = max(0, min_c - radius); c1 = min(nB, min_c + radius + 1)

    def _grow(r0, r1, c0, c1):
        while True:
            m_r0 = float(sub_AB[r0 - 1, c0:c1].mean()) if r0 > 0 else None
            m_r1 = float(sub_AB[r1, c0:c1].mean()) if r1 < nA else None
            m_c0 = float(sub_AB[r0:r1, c0 - 1].mean()) if c0 > 0 else None
            m_c1 = float(sub_AB[r0:r1, c1].mean()) if c1 < nB else None
            nr0 = r0 - 1 if (m_r0 is not None and m_r0 < soft_thresh) else r0
            nr1 = r1 + 1 if (m_r1 is not None and m_r1 < soft_thresh) else r1
            nc0 = c0 - 1 if (m_c0 is not None and m_c0 < soft_thresh) else c0
            nc1 = c1 + 1 if (m_c1 is not None and m_c1 < soft_thresh) else c1
            if (nr0, nr1, nc0, nc1) == (r0, r1, c0, c1):
                return r0, r1, c0, c1
            r0, r1, c0, c1 = nr0, nr1, nc0, nc1

    def _trim(r0, r1, c0, c1):
        while True:
            m_r0 = float(sub_AB[r0, c0:c1].mean())
            m_r1 = float(sub_AB[r1 - 1, c0:c1].mean())
            m_c0 = float(sub_AB[r0:r1, c0].mean())
            m_c1 = float(sub_AB[r0:r1, c1 - 1].mean())
            # never trim past the anchor, never empty the window
            nr0 = r0 + 1 if (r0 < min_r and m_r0 >= soft_thresh) else r0
            nr1 = r1 - 1 if (r1 - 1 > min_r and m_r1 >= soft_thresh) else r1
            nc0 = c0 + 1 if (c0 < min_c and m_c0 >= soft_thresh) else c0
            nc1 = c1 - 1 if (c1 - 1 > min_c and m_c1 >= soft_thresh) else c1
            if (nr0, nr1, nc0, nc1) == (r0, r1, c0, c1):
                return r0, r1, c0, c1
            r0, r1, c0, c1 = nr0, nr1, nc0, nc1

    for _ in range(20):                      # bounded; converges in 2-3 rounds
        prev = (r0, r1, c0, c1)
        r0, r1, c0, c1 = _grow(r0, r1, c0, c1)
        r0, r1, c0, c1 = _trim(r0, r1, c0, c1)
        if (r0, r1, c0, c1) == prev:
            break

    window = sub_AB[r0:r1, c0:c1]
    hmean = float(window.mean())
    hmin = float(window.min())
    hdens = float((window < 5.0).mean())
    hmcp = None
    if contact_AB is not None:
        try:
            hmcp = float(contact_AB[r0:r1, c0:c1].mean())
        except Exception:
            pass
    recip_mean = None
    if sub_BA is not None and sub_BA.shape[0] >= c1 and sub_BA.shape[1] >= r1:
        try:
            recip_mean = float(sub_BA[c0:c1, r0:r1].mean())
        except Exception:
            pass
    # Unknown reciprocity is neutral (1.0); only measurably poor reciprocity
    # is penalised.
    recip_factor = 0.7 if (recip_mean is not None
                           and recip_mean >= hmean * 1.5) else 1.0
    # 4. Size guard: a window that survives grow/trim can legitimately be very
    #    small, and a lone confident token pair would otherwise score ~1.0.
    #    Full credit requires at least MIN_CORE sub-5 A cells (a 5x5 patch);
    #    smaller cores are damped by sqrt of the shortfall.
    MIN_CORE = 25.0
    n_core = float((window < 5.0).sum())
    size_factor = min(1.0, (n_core / MIN_CORE) ** 0.5)
    hscore = max(0.0, min(1.0, (1.0 - hmean / 30.0) * (hdens ** 0.5)
                          * recip_factor * size_factor))
    return {
        'min_row': min_r, 'min_col': min_c,
        'hotspot_r0': r0, 'hotspot_r1': r1,
        'hotspot_c0': c0, 'hotspot_c1': c1,
        'hotspot_mean_pae': round(hmean, 2),
        'hotspot_min_pae': round(hmin, 2),
        'hotspot_contact_density': round(hdens, 3),
        'hotspot_mean_cp': round(hmcp, 3) if hmcp is not None else None,
        'hotspot_size_A': r1 - r0,
        'hotspot_size_B': c1 - c0,
        'hotspot_recip_mean': round(recip_mean, 2) if recip_mean is not None else None,
        'hotspot_core_cells': int(n_core),
        'hotspot_score': round(hscore, 3),
    }
