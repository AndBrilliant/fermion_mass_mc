"""
engine.koide_manifold — analytic direct sampling on the Koide 2/3 manifold.

The Koide condition Q(m_1, m_2, m_3) = 2/3 is one equation in three
unknowns, so the compatible set is a 2-dimensional manifold in mass
space. Rejection sampling with tight tolerance on |Q - 2/3| is ~10^-5
efficient; direct sampling via the Koide quadratic is exact to machine
precision and ~10^5 times faster.

Parameterisation used here (matches Brilliant 2026 Soddy paper §5.4):
    draw m_3  (the heaviest, tau-like)
    draw r = m_1 / m_3  (mass ratio, electron-like over heaviest)
    solve Koide quadratic for m_2  (the middle mass, muon-like)

Math derivation. Let u_i = sqrt(m_i). Koide Q = 2/3 is equivalent to
    3(u_1^2 + u_2^2 + u_3^2) = 2(u_1 + u_2 + u_3)^2
which rearranges to the symmetric form
    u_1^2 + u_2^2 + u_3^2 = 4(u_1 u_2 + u_2 u_3 + u_1 u_3).
Given u_1 and u_3, this is a quadratic in u_2:
    u_2^2 - 4(u_1 + u_3) u_2 + (u_1^2 + u_3^2 - 4 u_1 u_3) = 0
with discriminant 3(u_1^2 + 4 u_1 u_3 + u_3^2) (always positive for
positive u_1, u_3), giving two roots
    u_2 = 2(u_1 + u_3) +/- sqrt(3) * sqrt(u_1^2 + 4 u_1 u_3 + u_3^2).
The "-" root lies between u_1 and u_3 and corresponds to the physical
charged-lepton pattern (m_mu between m_e and m_tau). The "+" root
lies outside and corresponds to Koide's alternate quadratic solution
(~3.3 MeV for the real leptons; see Brilliant 2026 Soddy paper).
"""
import math
import numpy as np


def koide_solve_middle(m1: float, m3: float, branch: str = "minus") -> float:
    """Given m_1 and m_3, solve the Koide quadratic for m_2.

    Args:
        m1: lightest mass (GeV), e.g. electron-like.
        m3: heaviest mass (GeV), e.g. tau-like.  Must satisfy m3 > m1.
        branch: "minus" for the physical middle root (between m1 and m3),
                "plus" for the alternate (Koide's second root).

    Returns:
        m_2 in GeV on the Koide manifold, or nan if the discriminant
        yields a non-physical (m_2 <= 0 or complex) result.
    """
    if m1 <= 0 or m3 <= 0 or m3 <= m1:
        return float("nan")
    u1 = math.sqrt(m1)
    u3 = math.sqrt(m3)
    b = 2.0 * (u1 + u3)
    disc = u1 * u1 + 4.0 * u1 * u3 + u3 * u3
    if disc < 0:
        return float("nan")
    root = math.sqrt(3.0 * disc)
    u2 = b - root if branch == "minus" else b + root
    if u2 <= 0 or not math.isfinite(u2):
        return float("nan")
    return u2 * u2


def sample_koide_triple(rng: np.random.Generator,
                         prior_draw_heaviest,
                         prior_draw_ratio) -> tuple:
    """Draw a Koide-satisfying triple (m1, m2, m3) sorted ascending.

    Args:
        rng: numpy Generator (seeded by caller).
        prior_draw_heaviest: callable(rng) -> m_3 in GeV.
        prior_draw_ratio: callable(rng) -> m_1/m_3, expected < 1.

    Returns:
        (m_1, m_2, m_3) sorted ascending, in GeV. Returns (nan, nan, nan)
        if the draw is unphysical (e.g. ratio >= 1, Koide root negative).
    """
    m3 = float(prior_draw_heaviest(rng))
    r = float(prior_draw_ratio(rng))
    if not (0.0 < r < 1.0) or not (m3 > 0.0) or not math.isfinite(m3):
        return (float("nan"),) * 3
    m1 = r * m3
    m2 = koide_solve_middle(m1, m3, branch="minus")
    if math.isnan(m2):
        return (float("nan"),) * 3
    # Enforce ascending order (m2 should already satisfy m1 < m2 < m3 for
    # the minus branch, but numerics can fail near the boundary)
    if not (m1 < m2 < m3):
        return (float("nan"),) * 3
    return (m1, m2, m3)


def koide_residual(m1: float, m2: float, m3: float) -> float:
    """|Q - 2/3| for diagnostic checks.  Should be ~1e-15 for direct-
    sampled triples."""
    u = math.sqrt(m1) + math.sqrt(m2) + math.sqrt(m3)
    if u <= 0:
        return float("nan")
    Q = (m1 + m2 + m3) / (u * u)
    return abs(Q - 2.0 / 3.0)


if __name__ == "__main__":
    # Self-test: known leptons should produce near-exact m_mu when we
    # solve from (m_e, m_tau).
    m_e, m_tau = 0.00051099895, 1.77693
    m_mu_solved = koide_solve_middle(m_e, m_tau, branch="minus")
    print(f"Koide-solved m_mu from (m_e, m_tau): "
          f"{m_mu_solved * 1000:.4f} MeV   (PDG: 105.6584 MeV)")
    # Self-test: random draws.
    rng = np.random.default_rng(2026)
    draw_m3 = lambda r: np.exp(r.uniform(np.log(1e-3), np.log(200.0)))
    draw_rat = lambda r: np.exp(r.uniform(np.log(1e-3), np.log(1e-1)))
    n_good = 0
    for _ in range(1000):
        t = sample_koide_triple(rng, draw_m3, draw_rat)
        if not math.isnan(t[0]):
            res = koide_residual(*t)
            if res < 1e-12:
                n_good += 1
    print(f"Koide-compatible fraction of 1000 draws (res<1e-12): {n_good}/1000")
