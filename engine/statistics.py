"""
engine.statistics — tail probability, Clopper-Pearson CI, histograms.

All p-values reported here are empirical tail probabilities under the
specified null.  They are NOT severity assessments in the Mayo sense,
which require an alternative pdf and are distinct (Cousins 2018,
arXiv:1807.05996).  We use "tail probability" or "p-value" exclusively.
"""
import math
import numpy as np
from scipy.stats import beta as beta_dist


def clopper_pearson_ci(n_hits: int, n_trials: int,
                       confidence: float = 0.95) -> tuple:
    """Exact binomial confidence interval.  Returns (lo, hi) on the hit
    fraction."""
    if n_trials == 0:
        return (float("nan"), float("nan"))
    alpha = 1.0 - confidence
    lo = beta_dist.ppf(alpha / 2.0, n_hits, n_trials - n_hits + 1) \
         if n_hits > 0 else 0.0
    hi = beta_dist.ppf(1.0 - alpha / 2.0, n_hits + 1, n_trials - n_hits) \
         if n_hits < n_trials else 1.0
    return (lo, hi)


def tail_fraction(values: np.ndarray, threshold: float,
                  direction: str = "leq") -> float:
    """Fraction of `values` that satisfy the tail condition vs threshold.

    Args:
        values: array of per-draw test statistic values.
        threshold: the observed value (or band edge).
        direction: "leq" for values <= threshold (i.e. at least as
                   small as observed), "geq" for values >= threshold.

    Returns:
        fraction in [0, 1].  Ignores NaN values.
    """
    v = np.asarray(values)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return float("nan")
    if direction == "leq":
        return float(np.sum(v <= threshold)) / float(v.size)
    else:
        return float(np.sum(v >= threshold)) / float(v.size)


def hit_within_sigma(F2_sample: float, m_s_at_mu: float,
                     sigma_m_s_at_mu: float) -> bool:
    """Paper's hit criterion: |F^2 - m_s(mu_star)| < sigma_m_s(mu_star)."""
    if not all(map(math.isfinite, [F2_sample, m_s_at_mu, sigma_m_s_at_mu])):
        return False
    return abs(F2_sample - m_s_at_mu) < sigma_m_s_at_mu


def residual_distribution_summary(residuals: np.ndarray,
                                   observed: float = 0.038e-3) -> dict:
    """Summary statistics for a residual distribution.

    Args:
        residuals: per-draw |F^2 - m_s(mu_star)| values in GeV.
        observed: the observed residual in GeV (default +0.038 MeV).

    Returns:
        dict of summary stats and tail probabilities.
    """
    r = np.asarray(residuals)
    r = r[np.isfinite(r)]
    if r.size == 0:
        return {"n": 0}
    return {
        "n": int(r.size),
        "median_MeV": float(np.median(r)) * 1000,
        "p10_MeV": float(np.percentile(r, 10)) * 1000,
        "p90_MeV": float(np.percentile(r, 90)) * 1000,
        "p_leq_observed": tail_fraction(r, observed, "leq"),
        "observed_MeV": observed * 1000,
    }


def report_null_a_hits(F2_samples: np.ndarray,
                       m_s_samples: np.ndarray,
                       sigma_samples: np.ndarray) -> dict:
    """Counting hits under Null A (paper's §5.4 convention).

    Hit = |F^2 - m_s(mu_star)| < sigma_m_s(mu_star).
    """
    F = np.asarray(F2_samples)
    M = np.asarray(m_s_samples)
    S = np.asarray(sigma_samples)
    good = np.isfinite(F) & np.isfinite(M) & np.isfinite(S)
    F, M, S = F[good], M[good], S[good]
    if F.size == 0:
        return {"n_valid": 0}
    hits = np.abs(F - M) < S
    n_hits = int(np.sum(hits))
    n_valid = int(F.size)
    frac = n_hits / n_valid
    ci_lo, ci_hi = clopper_pearson_ci(n_hits, n_valid)
    return {
        "n_valid": n_valid,
        "n_hits": n_hits,
        "hit_fraction": frac,
        "ci95_lo": ci_lo,
        "ci95_hi": ci_hi,
    }


if __name__ == "__main__":
    # Self-test: paper reports 121 hits / 25485 valid = 0.4748%
    result = {
        "n_valid": 25485,
        "n_hits": 121,
        "hit_fraction": 121 / 25485,
    }
    lo, hi = clopper_pearson_ci(121, 25485)
    print(f"Paper figures: 121/25485 = {121/25485*100:.4f}% "
          f"[95% CI: {lo*100:.3f}, {hi*100:.3f}]")
    print(f"Paper reports: 0.4748% [95% CI: 0.394%, 0.567%]")
