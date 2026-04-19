"""
tests.test_b_bootstrap — Null B: measurement-noise parametric bootstrap.

Fix the real-world setup (actual PDG+FLAG central values and sigmas).
Perturb each input within its uncertainty.  Ask: is the observed
residual (+0.038 MeV) distinguishable from what measurement noise alone
would produce?

This is the null the Soddy paper's §5.4-5.6 MC is MISSING.  The v3 LEE
paper calls the corresponding null "Model B" for the quark-Koide claims.

Perturbation sources (pre-registered):
    m_e, m_mu, m_tau  -- PDG 2024 pole mass sigmas (independent Gaussian)
    m_s(2 GeV)         -- FLAG 2024 lattice sigma
    alpha_s(M_Z)       -- PDG 2024 sigma

Variants:
    B1 uncorrelated Gaussians, 1 sigma   (baseline)
    B2 correlated lattice inputs (rho=0.5 among lattice-based sigmas)
    B3 tighter (0.5 sigma)
    B4 looser (2 sigma)

Test statistics reported:
    - Fraction of draws with perturbed residual <= observed (+0.038 MeV)
    - Fraction of draws producing a "hit" under paper's |F^2 - m_s| < sigma

Usage:
    python -m tests.test_b_bootstrap
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from engine import algorithm as alg
from engine import running
from engine import statistics as stats


# ---------------------------------------------------------------------------
# Real-world central values and uncertainties
# ---------------------------------------------------------------------------
# PDG 2024 charged lepton pole masses (GeV) and sigmas
M_E = 0.51099895e-3
SIG_M_E = 0.00000015e-3

M_MU = 105.6583755e-3
SIG_M_MU = 0.0000023e-3

M_TAU = 1.77693           # PDG central
SIG_M_TAU = 0.00009       # PDG sigma (0.09 MeV)

# FLAG 2024 N_f = 2+1+1
M_S_2GEV = 0.09344
SIG_M_S_2GEV = 0.00068

# PDG 2024
ALPHA_S_MZ = 0.1180
SIG_ALPHA_S = 0.0009

# Observed residual (paper, corrected per audit A.1: F^2 - m_s at observed leptons)
OBSERVED_RESIDUAL_GEV = 0.038e-3

# Pre-registered
N_BOOTSTRAP = 2000
SEED_BASE = 2126


def sample_one_bootstrap(rng, sigma_scale: float = 1.0,
                         correlation_rho: float = 0.0) -> dict:
    """Draw one perturbed measurement.

    correlation_rho: common correlation among the three lattice-sensitive
    inputs (m_s, alpha_s, and we arbitrarily include m_tau via the
    "shared experimental systematics" proxy).  For the baseline rho=0.
    """
    if correlation_rho == 0.0:
        m_e = rng.normal(M_E, SIG_M_E * sigma_scale)
        m_mu = rng.normal(M_MU, SIG_M_MU * sigma_scale)
        m_tau = rng.normal(M_TAU, SIG_M_TAU * sigma_scale)
        m_s = rng.normal(M_S_2GEV, SIG_M_S_2GEV * sigma_scale)
        as_mz = rng.normal(ALPHA_S_MZ, SIG_ALPHA_S * sigma_scale)
    else:
        # Draw the three lattice-sensitive inputs from a correlated
        # Gaussian, with leptons independent (leptons are not
        # lattice-determined).
        r = correlation_rho
        cov = np.array([
            [1.0, r, r],
            [r, 1.0, r],
            [r, r, 1.0],
        ]) * (sigma_scale ** 2)
        latt = rng.multivariate_normal(mean=[0, 0, 0], cov=cov)
        m_s = M_S_2GEV + SIG_M_S_2GEV * latt[0]
        as_mz = ALPHA_S_MZ + SIG_ALPHA_S * latt[1]
        # (Third draw unused for now; reserved for future lattice sys)
        m_e = rng.normal(M_E, SIG_M_E * sigma_scale)
        m_mu = rng.normal(M_MU, SIG_M_MU * sigma_scale)
        m_tau = rng.normal(M_TAU, SIG_M_TAU * sigma_scale)

    # Guard against negative/zero draws (extreme tails)
    if min(m_e, m_mu, m_tau, m_s, as_mz) <= 0:
        return None

    # Compute F^2 from perturbed leptons
    F2 = alg.outer_soddy_sq(m_e, m_mu, m_tau)
    if not math.isfinite(F2):
        return None

    # Compute mu_star from perturbed leptons
    mu = alg.mu_star(m_e, m_mu, m_tau)

    # Run perturbed m_s to perturbed mu_star using perturbed alpha_s
    m_s_at_mu = running.run_mass(m_s, 2.0, mu, as_mz=as_mz)
    if not math.isfinite(m_s_at_mu):
        return None

    # Paper-style sigma (propagated lattice sigma to mu_star)
    _, sig_m_s_at_mu = running.m_s_at(mu, M_S_2GEV, SIG_M_S_2GEV, as_mz=as_mz)

    residual_signed = F2 - m_s_at_mu
    residual_abs = abs(residual_signed)
    return {
        "F2": F2,
        "mu_star": mu,
        "m_s_at_mu": m_s_at_mu,
        "sigma_m_s_at_mu": sig_m_s_at_mu,
        "residual_signed": residual_signed,
        "residual_abs": residual_abs,
    }


def run_variant(name: str, sigma_scale: float, correlation_rho: float,
                 seed: int, n: int = N_BOOTSTRAP) -> dict:
    """Run a single bootstrap variant."""
    rng = np.random.default_rng(seed)
    draws = []
    while len(draws) < n:
        d = sample_one_bootstrap(rng, sigma_scale, correlation_rho)
        if d is not None:
            draws.append(d)

    residuals_abs = np.array([d["residual_abs"] for d in draws])
    residuals_signed = np.array([d["residual_signed"] for d in draws])

    # P(perturbed residual <= observed)
    p_leq_observed = stats.tail_fraction(residuals_abs,
                                          OBSERVED_RESIDUAL_GEV, "leq")

    # P(hit under paper's |F^2 - m_s| < sigma)
    F2s = np.array([d["F2"] for d in draws])
    mss = np.array([d["m_s_at_mu"] for d in draws])
    sigs = np.array([d["sigma_m_s_at_mu"] for d in draws])
    hit_result = stats.report_null_a_hits(F2s, mss, sigs)

    return {
        "variant": name,
        "sigma_scale": sigma_scale,
        "correlation_rho": correlation_rho,
        "seed": seed,
        "n": len(draws),
        "p_leq_observed_residual": p_leq_observed,
        "hit_fraction": hit_result.get("hit_fraction"),
        "ci95_lo": hit_result.get("ci95_lo"),
        "ci95_hi": hit_result.get("ci95_hi"),
        "residual_abs_summary_MeV": {
            "median": float(np.median(residuals_abs)) * 1000,
            "p10":    float(np.percentile(residuals_abs, 10)) * 1000,
            "p90":    float(np.percentile(residuals_abs, 90)) * 1000,
        },
        "residual_signed_summary_MeV": {
            "mean":   float(np.mean(residuals_signed)) * 1000,
            "median": float(np.median(residuals_signed)) * 1000,
            "std":    float(np.std(residuals_signed)) * 1000,
        },
    }


def main():
    print("=" * 70)
    print("Null B: measurement-noise parametric bootstrap")
    print("=" * 70)
    print(f"N per variant: {N_BOOTSTRAP}")
    print(f"Observed residual: +{OBSERVED_RESIDUAL_GEV*1000:.3f} MeV")
    print(f"Seed base: {SEED_BASE}")
    print()

    variants = [
        ("B1_baseline_uncorrelated", 1.0, 0.0, SEED_BASE),
        ("B2_correlated_rho0.5",      1.0, 0.5, SEED_BASE + 1),
        ("B3_tighter_0.5sigma",       0.5, 0.0, SEED_BASE + 2),
        ("B4_looser_2sigma",          2.0, 0.0, SEED_BASE + 3),
    ]

    results = {}
    for name, scale, rho, seed in variants:
        print(f"Variant {name} (sigma_scale={scale}, rho={rho}, seed={seed})...")
        r = run_variant(name, scale, rho, seed)
        results[name] = r
        print(f"  P(perturbed residual <= observed):   "
              f"{r['p_leq_observed_residual']:.4f}")
        print(f"  Fraction producing paper-style hit:   "
              f"{r['hit_fraction']*100:.2f}%")
        print(f"  Residual |.| median: "
              f"{r['residual_abs_summary_MeV']['median']:.3f} MeV")
        print(f"  Residual signed mean: "
              f"{r['residual_signed_summary_MeV']['mean']:+.3f} +/- "
              f"{r['residual_signed_summary_MeV']['std']:.3f} MeV")
        print()

    out_dir = Path(__file__).resolve().parent.parent / "results"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "test_b_bootstrap.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"Results written to {out_path}")


if __name__ == "__main__":
    main()
