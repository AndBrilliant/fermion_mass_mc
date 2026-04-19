"""
tests.test_a_random_spectrum — Null A: random-spectrum Koide-conditioned MC.

For each prior shape (A1-A4), sample N Koide-valid triples, compute F^2
and mu_star per draw, run m_s^MSbar from 2 GeV to mu_star, count hits
within 1 sigma of m_s(mu_star).  Report tail fractions with exact
Clopper-Pearson CI and full residual histograms.

Pre-registered:
    seed_base = 2026
    N_sample = 10000 valid Koide triples per prior
    Cutoff: mu_star > 1.0 GeV (matches paper §5.4)
    Hit criterion: |F^2 - m_s(mu_star)| < sigma_m_s(mu_star)

Usage:
    python -m tests.test_a_random_spectrum
    (from the v4 directory)
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


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
N_SAMPLES_PER_PRIOR = 10000      # valid Koide triples per prior
MU_STAR_MIN_GEV = 1.0            # perturbative cutoff (paper §5.4)
MAX_ATTEMPTS_MULT = 50           # oversample factor; A3 needs ~42x
SEED_BASE = 2026                 # pre-registered


# Lattice input (FLAG 2024, 2+1+1)
M_S_2GEV_GEV = 0.09344
SIGMA_M_S_2GEV_GEV = 0.00068


def sample_one(rng, prior_draw_heaviest, prior_draw_ratio):
    """One MC draw.  Returns dict with triple, F^2, mu_star, m_s(mu_star),
    sigma_m_s(mu_star), residual."""
    triple = km.sample_koide_triple(rng, prior_draw_heaviest, prior_draw_ratio)
    m1, m2, m3 = triple
    if math.isnan(m1):
        return None
    mu = alg.mu_star(m1, m2, m3)
    if mu < MU_STAR_MIN_GEV:
        return None
    F2 = alg.outer_soddy_sq(m1, m2, m3)
    if not math.isfinite(F2):
        return None
    m_s, sig_m_s = running.m_s_at(mu, M_S_2GEV_GEV, SIGMA_M_S_2GEV_GEV)
    if not math.isfinite(m_s) or not math.isfinite(sig_m_s):
        return None
    residual = F2 - m_s
    return {
        "m1": m1, "m2": m2, "m3": m3,
        "mu_star": mu,
        "F2": F2,
        "m_s_mu": m_s,
        "sigma_m_s_mu": sig_m_s,
        "residual": residual,
    }


def run_one_prior(prior_name: str, draw_heavy, draw_ratio,
                   n_target: int, seed: int) -> dict:
    """Run Null A for a single prior."""
    rng = np.random.default_rng(seed)
    samples = []
    max_attempts = n_target * MAX_ATTEMPTS_MULT
    attempts = 0
    while len(samples) < n_target and attempts < max_attempts:
        attempts += 1
        s = sample_one(rng, draw_heavy, draw_ratio)
        if s is not None:
            samples.append(s)
    if len(samples) < n_target:
        print(f"  WARNING: {prior_name} only produced {len(samples)} valid "
              f"samples in {attempts} attempts (accept rate "
              f"{len(samples)/attempts*100:.2f}%)")
    F2s = np.array([s["F2"] for s in samples])
    m_ss = np.array([s["m_s_mu"] for s in samples])
    sigs = np.array([s["sigma_m_s_mu"] for s in samples])
    residuals = np.abs(F2s - m_ss)

    hit_result = stats.report_null_a_hits(F2s, m_ss, sigs)
    resid_summary = stats.residual_distribution_summary(residuals,
                                                         observed=0.038e-3)
    return {
        "prior": prior_name,
        "seed": seed,
        "n_attempts": attempts,
        "n_samples": len(samples),
        **hit_result,
        "residual_summary_MeV": resid_summary,
    }


def main():
    results = {}
    print("=" * 70)
    print("Null A: random-spectrum Koide-conditioned MC")
    print("=" * 70)
    print(f"N target per prior: {N_SAMPLES_PER_PRIOR}")
    print(f"mu_star cutoff:     {MU_STAR_MIN_GEV} GeV")
    print(f"m_s(2 GeV) input:   {M_S_2GEV_GEV*1000:.2f} +/- "
          f"{SIGMA_M_S_2GEV_GEV*1000:.2f} MeV")
    print(f"Seed base:          {SEED_BASE}")
    print(f"Max attempts mult:  {MAX_ATTEMPTS_MULT}")
    print()

    for i, (name, (dh, dr)) in enumerate(priors.PRIORS.items()):
        seed = SEED_BASE + i
        print(f"Prior {name} (seed {seed})...", flush=True)
        result = run_one_prior(name, dh, dr,
                                N_SAMPLES_PER_PRIOR, seed)
        results[name] = result
        print(f"  samples:      {result['n_samples']} "
              f"(from {result['n_attempts']} attempts)")
        print(f"  hit fraction: {result.get('hit_fraction', float('nan'))*100:.3f}%")
        print(f"  95% CI:       [{result.get('ci95_lo', 0)*100:.3f}, "
              f"{result.get('ci95_hi', 0)*100:.3f}]")
        rs = result["residual_summary_MeV"]
        print(f"  residual median: {rs.get('median_MeV', float('nan')):.2f} MeV")
        print(f"  P(residual <= observed +0.038 MeV): "
              f"{rs.get('p_leq_observed', float('nan')):.4f}")
        print(flush=True)

    # Save JSON
    out_dir = Path(__file__).resolve().parent.parent / "results"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "test_a_random_spectrum.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"Results written to {out_path}")

    # Fold-variation summary across priors
    fracs = [r["hit_fraction"] for r in results.values()
             if "hit_fraction" in r and r["hit_fraction"] is not None]
    if fracs:
        print(f"\nSummary across priors:")
        print(f"  min hit fraction: {min(fracs)*100:.3f}%")
        print(f"  max hit fraction: {max(fracs)*100:.3f}%")
        if min(fracs) > 0:
            print(f"  fold variation:   {max(fracs)/min(fracs):.2f}x")
        else:
            print(f"  (min is zero; fold variation undefined)")


if __name__ == "__main__":
    main()
