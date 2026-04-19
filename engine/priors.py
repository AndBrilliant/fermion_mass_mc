"""
engine.priors — prior-shape grid for the Koide-manifold Null A MC.

Each prior is a pair of callables (draw_heaviest, draw_ratio) that feed
koide_manifold.sample_koide_triple.  All priors operate on the same
parameterisation (m_3, r = m_1/m_3); what varies is the shape of the
distribution, not the parameterisation.  This lets us attribute any
observed sensitivity to the prior *shape* rather than to parameterisation
artifacts.

The four shapes are chosen from the HEP literature on Yukawa / mass
priors.  Each has a citation justifying it as a legitimate "random
universe" distribution for the purposes of a null test.

A1  Donoghue log-uniform  [PRIMARY]
    m_3 ~ log-U[1 MeV, 200 GeV], r ~ log-U[1e-3, 1e-1]
    Donoghue 1998, PRD 57 5499, hep-ph/9712333.  Scale-invariant
    distribution shown to best fit observed charged-fermion spectrum
    under an anthropic-free null.  Taken as the default Yukawa-landscape
    prior in the literature.

A2  Hall-Salem-Watari log-normal
    m_3 ~ log-N(geom mean of SM charged fermion masses, sigma=2.5),
    r ~ log-N(median SM ratio, sigma=1.5).
    Hall, Salem, Watari 2007, PRD 76 093001, arXiv:0707.3446.
    Gaussian-landscape alternative in log space.  Gives extra weight
    to "central" Yukawa values and suppresses extreme tails.

A3  Yukawa-anarchy SVD (Hall-Murayama-Weiner inspired)
    Draw a random 3x3 Yukawa matrix with entries ~ N(0, 1).
    Compute singular values, sort, identify min and max as m_1 and m_3.
    Rescale globally so m_3 sits at a Donoghue-drawn scale.
    Hall, Murayama, Weiner 2000, PRL 84 2572, hep-ph/9911341;
    methodology from de Gouvea-Murayama 2003, hep-ph/0301050.
    Natural "no-structure" prior on mass matrices.

A4  Linear-uniform stress test
    m_3 ~ U[1 MeV, 200 GeV], r ~ U[1e-3, 1e-1]
    Same ranges as A1 but uniform on the linear scale.  Included as a
    stress test: strongly disfavours the light-mass region and is the
    worst-case against scale-invariance.  Matches the paper's §5.5 "D"
    prior for continuity.
"""
import math
import numpy as np


# ---------------------------------------------------------------------------
# A1: Donoghue log-uniform (primary)
# ---------------------------------------------------------------------------
A1_M3_RANGE_GEV = (1e-3, 200.0)      # [1 MeV, 200 GeV]
A1_RATIO_RANGE = (1e-3, 1e-1)


def a1_draw_heaviest(rng: np.random.Generator) -> float:
    lo, hi = map(math.log, A1_M3_RANGE_GEV)
    return math.exp(rng.uniform(lo, hi))


def a1_draw_ratio(rng: np.random.Generator) -> float:
    lo, hi = map(math.log, A1_RATIO_RANGE)
    return math.exp(rng.uniform(lo, hi))


# ---------------------------------------------------------------------------
# A2: Hall-Salem-Watari log-normal
# ---------------------------------------------------------------------------
# SM charged-lepton+quark reference: geometric mean ~ 0.1 GeV, log-sigma ~2.5
A2_M3_LOGMEAN = math.log(0.1)       # ~100 MeV (geom. mean of SM masses)
A2_M3_LOGSIGMA = 2.5                # covers [1 MeV, 10 GeV] within ~1 sigma
A2_R_LOGMEAN = math.log(0.01)       # SM lepton median ratio m_e/m_tau ~ 3e-4
A2_R_LOGSIGMA = 1.5


def a2_draw_heaviest(rng: np.random.Generator) -> float:
    return math.exp(rng.normal(A2_M3_LOGMEAN, A2_M3_LOGSIGMA))


def a2_draw_ratio(rng: np.random.Generator) -> float:
    # Truncate to (0, 1) to avoid unphysical ratio
    for _ in range(100):
        r = math.exp(rng.normal(A2_R_LOGMEAN, A2_R_LOGSIGMA))
        if 0 < r < 1:
            return r
    return 0.01  # fallback after truncation failures


# ---------------------------------------------------------------------------
# A3: Yukawa-anarchy SVD (Hall-Murayama-Weiner inspired)
# ---------------------------------------------------------------------------
# Draw 3x3 Yukawa matrix entries from N(0, 1), take singular values as
# proxy for mass ratios, then globally scale so m_3 lies at a Donoghue
# draw.  This gives Koide-compatible (m_3, r) pairs whose internal ratio
# structure comes from SVD statistics, not from a user-chosen range.
def _draw_svd_ratio(rng: np.random.Generator) -> float:
    """Return s_min/s_max of a 3x3 N(0,1) Yukawa matrix."""
    for _ in range(10):
        Y = rng.normal(0.0, 1.0, size=(3, 3))
        s = np.linalg.svd(Y, compute_uv=False)
        s_min, s_max = s.min(), s.max()
        if s_max > 0:
            r = s_min / s_max
            # Square because m ~ Y^2 (Yukawa to mass is linear in vev;
            # we treat singular values as mass ratios directly here).
            if 1e-6 < r < 1.0:
                return r
    return 1e-3


def a3_draw_heaviest(rng: np.random.Generator) -> float:
    # Same scale as A1 so A3 is directly comparable
    return a1_draw_heaviest(rng)


def a3_draw_ratio(rng: np.random.Generator) -> float:
    return _draw_svd_ratio(rng)


# ---------------------------------------------------------------------------
# A4: Linear-uniform stress test (paper's prior D)
# ---------------------------------------------------------------------------
A4_M3_RANGE_GEV = (1e-3, 200.0)
A4_RATIO_RANGE = (1e-3, 1e-1)


def a4_draw_heaviest(rng: np.random.Generator) -> float:
    return rng.uniform(*A4_M3_RANGE_GEV)


def a4_draw_ratio(rng: np.random.Generator) -> float:
    return rng.uniform(*A4_RATIO_RANGE)


# ---------------------------------------------------------------------------
# Registry for iterating over all priors
# ---------------------------------------------------------------------------
PRIORS = {
    "A1_donoghue_loguniform": (a1_draw_heaviest, a1_draw_ratio),
    "A2_hsw_lognormal":       (a2_draw_heaviest, a2_draw_ratio),
    "A3_anarchy_svd":         (a3_draw_heaviest, a3_draw_ratio),
    "A4_linear_uniform":      (a4_draw_heaviest, a4_draw_ratio),
}


if __name__ == "__main__":
    rng = np.random.default_rng(2026)
    for name, (dm, dr) in PRIORS.items():
        m3_samples = [dm(rng) for _ in range(10000)]
        r_samples = [dr(rng) for _ in range(10000)]
        print(f"{name}:")
        print(f"  m_3: median={np.median(m3_samples):.3e} "
              f"GeV, 5-95% = [{np.percentile(m3_samples,5):.2e}, "
              f"{np.percentile(m3_samples,95):.2e}]")
        print(f"  r:   median={np.median(r_samples):.3e}, "
              f"5-95% = [{np.percentile(r_samples,5):.2e}, "
              f"{np.percentile(r_samples,95):.2e}]")
