"""
tests.test_filter_sensitivity — §5.5 lower-cutoff filter-sensitivity scan.

Audit-requested artifact for the paper's §5.5 filter-sensitivity table.
Scans the A1 prior hit fraction as a function of mu_star > mu_min for
mu_min in {0.7, 1.0, 1.5, 2.0} GeV.

Pre-registered (matches paper §5.5):
    seed = 2026 (same as A1 baseline in test_a_random_spectrum.py)
    N = 10000 valid Koide triples
    Prior = A1 Donoghue log-uniform
    Cutoffs = [0.7, 1.0, 1.5, 2.0] GeV

Usage:
    python -m tests.test_filter_sensitivity
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from engine import algorithm as alg
from engine import koide_manifold as km
from engine import priors
from engine import running
from engine import statistics as stats


N_SAMPLES = 10000
SEED = 2026
MAX_ATTEMPTS_MULT = 50
CUTOFFS_GEV = [0.7, 1.0, 1.5, 2.0]

M_S_2GEV_GEV = 0.09344
SIGMA_M_S_2GEV_GEV = 0.00068


def sample_one_raw(rng, draw_heavy, draw_ratio):
    """Sample one Koide-valid triple and return F^2, mu_star, residual
    without applying any mu_star cutoff."""
    triple = km.sample_koide_triple(rng, draw_heavy, draw_ratio)
    m1, m2, m3 = triple
    if math.isnan(m1):
        return None
    mu = alg.mu_star(m1, m2, m3)
    # Require mu_star > 0.5 GeV to avoid unphysical running. We defer the
    # per-cutoff filter to the count loop; only reject triples that are
    # structurally unusable.
    if mu < 0.5:
        return None
    F2 = alg.outer_soddy_sq(m1, m2, m3)
    if not math.isfinite(F2):
        return None
    try:
        m_s, sig_m_s = running.m_s_at(mu, M_S_2GEV_GEV, SIGMA_M_S_2GEV_GEV)
    except Exception:
        return None
    if not math.isfinite(m_s) or not math.isfinite(sig_m_s):
        return None
    return {
        "mu_star": mu,
        "F2": F2,
        "m_s_mu": m_s,
        "sigma_m_s_mu": sig_m_s,
    }


def main():
    print("=" * 70)
    print("Filter sensitivity scan (paper §5.5)")
    print("=" * 70)
    print(f"N target:           {N_SAMPLES}")
    print(f"Seed:               {SEED}")
    print(f"Prior:              A1 Donoghue log-uniform")
    print(f"Cutoffs (GeV):      {CUTOFFS_GEV}")
    print()

    # Draw samples with a relaxed internal cutoff (0.5 GeV) for the raw pool,
    # then re-filter per cutoff. For each requested cutoff we run a
    # dedicated sampling loop until we reach N_SAMPLES valid triples under
    # that cutoff, using the same seed (so the 1.0 GeV row reproduces
    # test_a's A1 hit count exactly).

    draw_h, draw_r = priors.PRIORS["A1_donoghue_loguniform"]

    results = {}
    for cut in CUTOFFS_GEV:
        rng = np.random.default_rng(SEED)
        samples = []
        attempts = 0
        max_attempts = N_SAMPLES * MAX_ATTEMPTS_MULT
        while len(samples) < N_SAMPLES and attempts < max_attempts:
            attempts += 1
            s = sample_one_raw(rng, draw_h, draw_r)
            if s is None:
                continue
            if s["mu_star"] < cut:
                continue
            samples.append(s)
        if len(samples) < N_SAMPLES:
            print(f"  cutoff {cut} GeV: only {len(samples)} valid samples in {attempts} attempts")

        F2s = np.array([s["F2"] for s in samples])
        m_ss = np.array([s["m_s_mu"] for s in samples])
        sigs = np.array([s["sigma_m_s_mu"] for s in samples])
        hit_res = stats.report_null_a_hits(F2s, m_ss, sigs)

        results[f"mu_min_{cut:.1f}_GeV"] = {
            "cutoff_GeV": cut,
            "seed": SEED,
            "n_attempts": attempts,
            "n_samples": len(samples),
            **hit_res,
        }
        print(f"  cutoff {cut:.1f} GeV: "
              f"{hit_res['n_hits']}/{hit_res['n_valid']} = "
              f"{hit_res['hit_fraction']*100:.3f}%  "
              f"CI95 [{hit_res['ci95_lo']*100:.3f}, {hit_res['ci95_hi']*100:.3f}]")

    out_path = Path(__file__).resolve().parent.parent / "results" / "test_filter_sensitivity.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults written to {out_path}")


if __name__ == "__main__":
    main()
