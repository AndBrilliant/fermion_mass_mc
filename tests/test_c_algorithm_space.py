"""
tests.test_c_algorithm_space — algorithm-space check for paper Appendix.

Defends the paper's §3 claim that F^2 = (k_4^-)^2 is not selected from
a catalog of symmetric polynomials of the lepton triple.  Enumerates
~25 similar-complexity functions of (m_e, m_mu, m_tau), evaluates each
on the observed leptons, and counts how many land within 1 sigma of
m_s^MSbar(mu_star).

If only outer_soddy_sq (and its surrogate e1_minus_sqrt_p2) hit, the
paper's argument is strengthened.  If multiple functions hit, that is
useful information about the algorithm-selection freedom.

We also run each variant across a random-spectrum pseudo-universe
sample (Null A style) to compare the hit rate ACROSS algorithms
against the hit rate for outer_soddy_sq alone.  A variant that hits
on observed leptons but has a high hit rate under Null A would be
"algorithm-space noise" rather than a signal.

Pre-registered seed: 2226

Usage:
    python -m tests.test_c_algorithm_space
"""
import json
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


# Observed leptons
M_E = 0.51099895e-3
M_MU = 105.6583755e-3
M_TAU = 1.77693

# FLAG 2024 m_s(2 GeV)
M_S_2GEV = 0.09344
SIGMA_M_S_2GEV = 0.00068

# Perturbative cutoff
MU_STAR_MIN_GEV = 1.0

# Pre-registered
SEED = 2226
N_NULLA_SAMPLES = 2000  # per-variant for cross-algorithm null (smaller N
                         # is fine here — we're comparing algorithms, not
                         # pinning one to <1% precision)


def observed_hit_check() -> dict:
    """Part 1: evaluate every variant on the observed lepton triple.
    Count how many land within 1 sigma of m_s(mu_star)."""
    mu = alg.mu_star(M_E, M_MU, M_TAU)
    m_s_mu, sig_m_s_mu = running.m_s_at(mu, M_S_2GEV, SIGMA_M_S_2GEV)
    print(f"Observed leptons: mu_star = {mu*1000:.2f} MeV")
    print(f"m_s(mu_star) = {m_s_mu*1000:.3f} +/- {sig_m_s_mu*1000:.3f} MeV")
    print(f"1-sigma hit band: "
          f"[{(m_s_mu-sig_m_s_mu)*1000:.3f}, "
          f"{(m_s_mu+sig_m_s_mu)*1000:.3f}] MeV")
    print()

    vals = av.evaluate_all_variants(M_E, M_MU, M_TAU)
    entries = []
    for name, v in vals.items():
        if not math.isfinite(v):
            entries.append({"name": name, "value_MeV": None,
                            "n_sigma": None, "hit": False})
            continue
        n_sig = (v - m_s_mu) / sig_m_s_mu
        hit = abs(v - m_s_mu) < sig_m_s_mu
        entries.append({
            "name": name,
            "value_MeV": v * 1000,
            "n_sigma": n_sig,
            "hit": bool(hit),
        })
    entries.sort(key=lambda e: (e["n_sigma"] is None, abs(e["n_sigma"] or 1e18)))
    return {
        "mu_star_MeV": mu * 1000,
        "m_s_at_mu_MeV": m_s_mu * 1000,
        "sigma_m_s_MeV": sig_m_s_mu * 1000,
        "variants": entries,
    }


def per_algorithm_null_rate(seed: int = SEED,
                             n_samples: int = N_NULLA_SAMPLES) -> dict:
    """Part 2: for each variant, how often does it produce a hit on a
    random-spectrum pseudo-universe?  High rate = variant is noisy,
    landing near m_s doesn't mean much for that variant.  Low rate =
    variant is genuinely discriminating."""
    rng = np.random.default_rng(seed)
    # Use A1 Donoghue log-uniform as reference prior
    dh, dr = priors.PRIORS["A1_donoghue_loguniform"]

    # Collect triples first
    triples = []
    attempts = 0
    max_attempts = n_samples * 30
    while len(triples) < n_samples and attempts < max_attempts:
        attempts += 1
        t = km.sample_koide_triple(rng, dh, dr)
        if math.isnan(t[0]):
            continue
        mu = alg.mu_star(*t)
        if mu < MU_STAR_MIN_GEV:
            continue
        m_s, sig = running.m_s_at(mu, M_S_2GEV, SIGMA_M_S_2GEV)
        if not math.isfinite(m_s):
            continue
        triples.append((t, mu, m_s, sig))

    # Hit rate per algorithm variant
    variant_stats = {}
    for name, func, _note in av.ALGORITHM_VARIANTS:
        hits = 0
        finite = 0
        for (t, mu, m_s, sig) in triples:
            try:
                v = func(*t)
            except (ValueError, ZeroDivisionError, OverflowError):
                continue
            if not math.isfinite(v):
                continue
            finite += 1
            if abs(v - m_s) < sig:
                hits += 1
        if finite == 0:
            variant_stats[name] = {"hit_rate": None, "n_finite": 0}
            continue
        ci_lo, ci_hi = stats.clopper_pearson_ci(hits, finite)
        variant_stats[name] = {
            "n_finite": finite,
            "n_hits": hits,
            "hit_rate": hits / finite,
            "ci95_lo": ci_lo,
            "ci95_hi": ci_hi,
        }
    return {
        "seed": seed,
        "n_triples_used": len(triples),
        "n_triples_target": n_samples,
        "per_variant": variant_stats,
    }


def main():
    print("=" * 70)
    print("Algorithm-space check (paper Appendix)")
    print("=" * 70)
    print()

    print("Part 1: observed-lepton hit count")
    print("-" * 70)
    observed = observed_hit_check()
    n_hits = sum(1 for e in observed["variants"] if e["hit"])
    n_finite = sum(1 for e in observed["variants"]
                    if e["value_MeV"] is not None)
    n_total = len(observed["variants"])
    print(f"Variants evaluated: {n_total} (finite: {n_finite})")
    print(f"Variants landing within 1 sigma of m_s(mu_star): {n_hits}")
    print()
    print(f"{'Variant':<30} {'Value (MeV)':>12} {'n-sigma':>10} {'Hit?':>6}")
    print("-" * 60)
    for e in observed["variants"]:
        v = e["value_MeV"]
        ns = e["n_sigma"]
        vstr = f"{v:>12.3f}" if v is not None else f"{'(nan)':>12}"
        nsstr = f"{ns:>+10.2f}" if ns is not None else f"{'(nan)':>10}"
        hit_str = "YES" if e["hit"] else "no"
        print(f"{e['name']:<30} {vstr} {nsstr} {hit_str:>6}")
    print()

    print("Part 2: per-algorithm Null A hit rate (random universes)")
    print("-" * 70)
    per_algo = per_algorithm_null_rate()
    print(f"Random-spectrum triples used: {per_algo['n_triples_used']}")
    print()
    print(f"{'Variant':<30} {'Hit rate':>10} {'95% CI':>20}")
    print("-" * 60)
    for name in sorted(per_algo["per_variant"].keys()):
        s = per_algo["per_variant"][name]
        if s.get("hit_rate") is None:
            print(f"{name:<30} {'(nan)':>10}")
            continue
        rate = s["hit_rate"] * 100
        lo = s["ci95_lo"] * 100
        hi = s["ci95_hi"] * 100
        print(f"{name:<30} {rate:>9.3f}% [{lo:>5.3f}, {hi:>5.3f}]")
    print()

    # Save
    out_dir = Path(__file__).resolve().parent.parent / "results"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "test_c_algorithm_space.json"
    with open(out_path, "w") as f:
        json.dump({
            "observed": observed,
            "per_algorithm_null_a": per_algo,
        }, f, indent=2, default=str)
    print(f"Results written to {out_path}")


if __name__ == "__main__":
    main()
