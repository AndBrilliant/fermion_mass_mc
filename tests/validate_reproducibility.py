"""
tests.validate_reproducibility — smoke test.

Runs all engine modules' self-tests plus a small-N version of each
test_a, test_b, test_c to verify the pipeline works end-to-end before
committing to a full production run.

Expected output:
    - All engine self-tests pass (observed-lepton F^2 ~ 95.11 MeV, etc.)
    - Null A produces some hits (rates ~0.1%-1%) across priors
    - Null B produces moderate tail fractions (~0.3-0.6)
    - Algorithm-space enumeration runs without exceptions

Usage:
    python -m tests.validate_reproducibility
"""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from engine import algorithm as alg
from engine import algorithm_variants as av
from engine import koide_manifold as km
from engine import priors
from engine import running
from engine import statistics as stats


def section(title):
    print()
    print("=" * 70)
    print(title)
    print("=" * 70)


def main():
    section("1. Engine: Soddy algorithm on observed leptons")
    m_e, m_mu, m_tau = 0.51099895e-3, 105.6583755e-3, 1.77693
    p = alg.soddy_pipeline(m_e, m_mu, m_tau)
    F2_MeV = p["F2"] * 1000
    mu_star_MeV = p["mu_star"] * 1000
    print(f"F^2 (outer Soddy squared): {F2_MeV:.4f} MeV")
    print(f"Expected (paper):          95.1134 MeV")
    assert abs(F2_MeV - 95.11) < 0.01, f"F^2 mismatch: {F2_MeV}"
    print(f"mu_star: {mu_star_MeV:.2f} MeV   (expected ~1883.10)")
    assert abs(mu_star_MeV - 1883.1) < 1.0

    section("2. Engine: m_s RG running, 2 GeV -> mu_star")
    m_s_mu, sig_mu = running.m_s_at(p["mu_star"])
    print(f"m_s({mu_star_MeV:.2f} MeV) = {m_s_mu*1000:.3f} +/- "
          f"{sig_mu*1000:.3f} MeV")
    print(f"Expected (paper):          95.07 +/- 0.69 MeV")
    assert abs(m_s_mu * 1000 - 95.07) < 2.0

    section("3. Engine: Koide manifold direct sampler")
    rng = np.random.default_rng(42)
    n_valid = 0
    for _ in range(500):
        t = km.sample_koide_triple(rng,
                                    priors.a1_draw_heaviest,
                                    priors.a1_draw_ratio)
        if math.isnan(t[0]):
            continue
        res = km.koide_residual(*t)
        if res < 1e-10:
            n_valid += 1
    print(f"Koide-valid triples: {n_valid} / 500   (res < 1e-10)")
    # ~58% expected due to m_3/m_1 > 67.9 constraint on the "minus" branch.
    assert n_valid > 200, f"Direct sampler produced only {n_valid}/500"

    section("4. Null A: quick run, N=500 per prior")
    from tests import test_a_random_spectrum as ta
    ta.N_SAMPLES_PER_PRIOR = 500
    results_a = {}
    for i, (name, (dh, dr)) in enumerate(priors.PRIORS.items()):
        print(f"  {name}...")
        r = ta.run_one_prior(name, dh, dr, 500, ta.SEED_BASE + i)
        results_a[name] = r
        print(f"    hit_fraction = {r.get('hit_fraction', 0)*100:.2f}% "
              f"({r['n_hits']}/{r['n_samples']})")

    section("5. Null B: quick run, N=200 per variant")
    from tests import test_b_bootstrap as tb
    variants = [
        ("B1_baseline", 1.0, 0.0, tb.SEED_BASE),
        ("B3_tighter",  0.5, 0.0, tb.SEED_BASE + 2),
    ]
    for name, scale, rho, seed in variants:
        print(f"  {name}...")
        r = tb.run_variant(name, scale, rho, seed, n=200)
        print(f"    P(perturbed residual <= observed) = "
              f"{r['p_leq_observed_residual']:.3f}")
        print(f"    hit_fraction = {r['hit_fraction']*100:.2f}%")

    section("6. Algorithm-space: observed-lepton hit count")
    from tests import test_c_algorithm_space as tc
    obs = tc.observed_hit_check()
    n_hits = sum(1 for e in obs["variants"] if e["hit"])
    print(f"Algorithm variants hitting within 1 sigma: {n_hits} of "
          f"{len(obs['variants'])}")

    section("SMOKE TEST PASSED")


if __name__ == "__main__":
    main()
