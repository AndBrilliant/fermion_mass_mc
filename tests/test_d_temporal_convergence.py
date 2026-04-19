"""
tests.test_d_temporal_convergence -- temporal convergence MC.

Rationale: as experimental precision improves, a real regularity should
produce a residual that shrinks in MeV while staying roughly constant
in units of sigma (i.e., F^2 tracks m_s as m_s gets pinned down).  A
coincidence at one epoch should produce a residual that is STABLE in
MeV while sigma narrows around it -- i.e. residual/sigma grows over
time as the measurement tightens onto a central value different from
F^2.

This is the third discriminator in the v3 LEE paper section 6.3 and in
the Soddy paper section 1 closing paragraph ("the discriminatory power
to notice it has only recently arrived"). This test makes that
argument quantitative.

Inputs: historical m_s(2 GeV) central values and sigmas (FLAG 2013 -
FLAG 2024) from results/historical_masses.json.  Lepton masses are
essentially stationary across the same epoch (sigma_{F^2} ~ 0.01 MeV),
so F^2 is treated as fixed at the PDG 2024 central value.

Test statistics (all descriptive):

1. Trajectory table.  For each epoch: residual (MeV), sigma_m_s (MeV
   at that epoch, evolved to mu_star), residual/sigma.

2. Weighted-mean residual across epochs under the "all epochs sample
   a common true m_s^true" model.

3. Chi-square against F^2 = 95.113 MeV hypothesis.  Under the
   hypothesis that m_s^true = F^2, chi^2 = sum of (r_i/sigma_i)^2
   follows the chi-square distribution with df = N_indep.

4. MC under the coincidence null.  Draw m_s^true uniformly from
   [60, 130] MeV (historical range of all determinations).  Report
   fraction of MC draws where F^2 lands within 1 sigma at EVERY epoch.

Pre-registered seed: 2326.  N_MC = 50,000.

Usage:
    python -m tests.test_d_temporal_convergence
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import chi2, norm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from engine import algorithm as alg
from engine import running
from engine import statistics as stats


# ---- Observed lepton inputs (PDG 2024) ------------------------------
M_E, M_MU, M_TAU = 0.51099895e-3, 105.6583755e-3, 1.77693

F2_FIXED_GEV = alg.outer_soddy_sq(M_E, M_MU, M_TAU)
MU_STAR_FIXED = alg.mu_star(M_E, M_MU, M_TAU)

M_S_PRIOR_LO_MEV, M_S_PRIOR_HI_MEV = 60.0, 130.0

N_MC = 50000
SEED = 2326

HISTORY_PATH = Path(__file__).resolve().parent.parent / "results" / \
               "historical_masses.json"


def load_history():
    with open(HISTORY_PATH) as f:
        hist = json.load(f)
    epochs = []
    for key, entry in hist["flag_ms_history"].items():
        if key == "note":
            continue
        epochs.append({
            "label": key,
            "year": entry["year"],
            "m_s_2gev_MeV": entry["m_s_MeV"],
            "sigma_2gev_MeV": entry["m_s_sigma_MeV"],
            "nf": entry.get("nf", "2+1+1"),
        })
    epochs.sort(key=lambda e: e["year"])
    return epochs


def evolve_epoch_to_mustar(epoch):
    m_s_mustar_GeV, sigma_mustar_GeV = running.m_s_at(
        MU_STAR_FIXED,
        epoch["m_s_2gev_MeV"] * 1e-3,
        epoch["sigma_2gev_MeV"] * 1e-3,
    )
    residual_GeV = F2_FIXED_GEV - m_s_mustar_GeV
    return {
        **epoch,
        "m_s_mustar_MeV": m_s_mustar_GeV * 1000,
        "sigma_mustar_MeV": sigma_mustar_GeV * 1000,
        "residual_MeV": residual_GeV * 1000,
        "residual_over_sigma": residual_GeV / sigma_mustar_GeV,
        "within_1sigma": abs(residual_GeV) < sigma_mustar_GeV,
    }


def chi2_test(epochs):
    residuals = [e["residual_over_sigma"] for e in epochs]
    chi2_stat = sum(r * r for r in residuals)
    N = len(residuals)
    independent_indices = [i for i, e in enumerate(epochs)
                           if e["year"] in (2013, 2016, 2024)]
    N_indep = len(independent_indices)
    chi2_indep = sum(residuals[i] ** 2 for i in independent_indices)
    return {
        "chi2_all_epochs": chi2_stat,
        "N_all_epochs": N,
        "p_value_all_epochs": float(1 - chi2.cdf(chi2_stat, df=N)),
        "chi2_independent": chi2_indep,
        "N_independent": N_indep,
        "p_value_independent": float(1 - chi2.cdf(chi2_indep, df=N_indep)),
        "interpretation": (
            "Small chi^2 means residuals are consistent with or smaller "
            "than expected under m_s^true = F^2.  Large chi^2 means "
            "residuals are inconsistent with that hypothesis."
        ),
    }


def weighted_mean_residual(epochs):
    weights = np.array([1.0 / (e["sigma_mustar_MeV"] ** 2) for e in epochs])
    residuals = np.array([e["residual_MeV"] for e in epochs])
    w_mean = float(np.sum(weights * residuals) / np.sum(weights))
    w_sigma = float(1.0 / math.sqrt(np.sum(weights)))
    return {
        "weighted_mean_MeV": w_mean,
        "weighted_sigma_MeV": w_sigma,
        "n_sigma_from_zero": w_mean / w_sigma if w_sigma > 0 else float("nan"),
    }


def mc_coincidence_null(epochs, n_mc=N_MC, seed=SEED):
    rng = np.random.default_rng(seed)
    F2_MeV = F2_FIXED_GEV * 1000
    sigmas = np.array([e["sigma_mustar_MeV"] for e in epochs])
    n_epochs = len(sigmas)
    hits_all_epochs = 0
    hits_per_epoch = np.zeros(n_epochs, dtype=int)
    n_observed_hits = sum(1 for e in epochs if e["within_1sigma"])
    hits_at_least_observed = 0

    for _ in range(n_mc):
        m_s_true = rng.uniform(M_S_PRIOR_LO_MEV, M_S_PRIOR_HI_MEV)
        measured = rng.normal(m_s_true, sigmas)
        residuals = F2_MeV - measured
        within = np.abs(residuals) < sigmas
        hits_per_epoch += within.astype(int)
        n_hits_this = int(np.sum(within))
        if n_hits_this == n_epochs:
            hits_all_epochs += 1
        if n_hits_this >= n_observed_hits:
            hits_at_least_observed += 1

    return {
        "n_mc": n_mc,
        "n_observed_hits": n_observed_hits,
        "n_epochs": n_epochs,
        "p_all_epochs_hit_under_null": hits_all_epochs / n_mc,
        "p_at_least_observed_hits_under_null": hits_at_least_observed / n_mc,
        "per_epoch_hit_rate_under_null": (hits_per_epoch / n_mc).tolist(),
    }


def main():
    print("=" * 70)
    print("Null D: temporal convergence of F^2 vs historical m_s(2 GeV)")
    print("=" * 70)
    print(f"F^2 (PDG 2024 leptons): {F2_FIXED_GEV*1000:.4f} MeV")
    print(f"mu_star: {MU_STAR_FIXED*1000:.2f} MeV")
    print()

    epochs_raw = load_history()
    epochs = [evolve_epoch_to_mustar(e) for e in epochs_raw]

    print("Trajectory (FLAG m_s(2 GeV) evolved to mu_star):")
    print(f"{'Epoch':<14} {'m_s(mu*)':>10} {'sigma':>8} {'resid':>8} "
          f"{'r/sig':>7} {'<1sig?':>8}")
    print("-" * 60)
    for e in epochs:
        hit = "YES" if e["within_1sigma"] else "no"
        print(f"{e['label']:<14} {e['m_s_mustar_MeV']:>10.3f} "
              f"{e['sigma_mustar_MeV']:>8.3f} {e['residual_MeV']:>+8.3f} "
              f"{e['residual_over_sigma']:>+7.3f} {hit:>8}")
    print()

    chi = chi2_test(epochs)
    print(f"Chi-square against m_s^true = F^2 hypothesis:")
    print(f"  All {chi['N_all_epochs']} epochs:  chi^2 = "
          f"{chi['chi2_all_epochs']:.3f},  "
          f"P(chi^2 >= obs) = {chi['p_value_all_epochs']:.3f}")
    print(f"  {chi['N_independent']} independent: chi^2 = "
          f"{chi['chi2_independent']:.3f},  "
          f"P(chi^2 >= obs) = {chi['p_value_independent']:.3f}")
    print()

    wm = weighted_mean_residual(epochs)
    print(f"Inverse-variance weighted mean residual:")
    print(f"  <residual> = {wm['weighted_mean_MeV']:+.4f} +/- "
          f"{wm['weighted_sigma_MeV']:.4f} MeV "
          f"({wm['n_sigma_from_zero']:+.2f} sigma from zero)")
    print()

    mc = mc_coincidence_null(epochs)
    print(f"MC coincidence null (N={mc['n_mc']}, m_s^true ~ U(60,130) MeV):")
    print(f"  P(F^2 lands in 1-sig band at ALL {mc['n_epochs']} epochs) = "
          f"{mc['p_all_epochs_hit_under_null']:.4f}")
    print(f"  Observed: F^2 lands in 1-sig band at "
          f"{mc['n_observed_hits']}/{mc['n_epochs']} epochs")
    print(f"  P(at-least-observed-hits under null) = "
          f"{mc['p_at_least_observed_hits_under_null']:.4f}")
    print(f"  Per-epoch hit rate under null: "
          f"{[f'{r:.3f}' for r in mc['per_epoch_hit_rate_under_null']]}")
    print()

    out = {
        "F2_MeV": F2_FIXED_GEV * 1000,
        "mu_star_MeV": MU_STAR_FIXED * 1000,
        "trajectory": epochs,
        "chi2_test": chi,
        "weighted_mean": wm,
        "mc_coincidence_null": mc,
    }
    out_path = Path(__file__).resolve().parent.parent / "results" / \
               "test_d_temporal_convergence.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"Results written to {out_path}")


if __name__ == "__main__":
    main()
