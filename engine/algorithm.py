"""
engine.algorithm â the Soddy pipeline.

Given a Koide-saturating triple (m_1, m_2, m_3), compute:
  - e_1 = sum sqrt(m_i)
  - e_2 = sum_{i<j} sqrt(m_i m_j)
  - outer Soddy curvature k_4^- = e_1 - 2 sqrt(e_2)
  - F^2 = (k_4^-)^2   [the comparison quantity, in units of mass]
  - mu_star = m_1 + m_2 + m_3  [the natural lepton-sum evaluation scale]

Reference: Brilliant 2026 "A Descartes–Soddy completion of the Koide lepton triple", Proposition 1.
"""
import math


def e1(m1: float, m2: float, m3: float) -> float:
    """Sum of mass square-roots."""
    return math.sqrt(m1) + math.sqrt(m2) + math.sqrt(m3)


def e2(m1: float, m2: float, m3: float) -> float:
    """Sum of pairwise sqrt products."""
    return (math.sqrt(m1 * m2) + math.sqrt(m2 * m3) + math.sqrt(m1 * m3))


def outer_soddy_curvature(m1: float, m2: float, m3: float) -> float:
    """k_4^- = e_1 - 2 sqrt(e_2) in units of sqrt(mass)."""
    E2 = e2(m1, m2, m3)
    if E2 < 0:
        return float("nan")
    return e1(m1, m2, m3) - 2.0 * math.sqrt(E2)


def outer_soddy_sq(m1: float, m2: float, m3: float) -> float:
    """F^2 = (k_4^-)^2 in units of mass (GeV).  The comparison quantity."""
    k = outer_soddy_curvature(m1, m2, m3)
    if not math.isfinite(k):
        return float("nan")
    return k * k


def inner_soddy_curvature(m1: float, m2: float, m3: float) -> float:
    """k_4^+ = e_1 + 2 sqrt(e_2).  The other Descartes root."""
    E2 = e2(m1, m2, m3)
    if E2 < 0:
        return float("nan")
    return e1(m1, m2, m3) + 2.0 * math.sqrt(E2)


def inner_soddy_sq(m1: float, m2: float, m3: float) -> float:
    k = inner_soddy_curvature(m1, m2, m3)
    if not math.isfinite(k):
        return float("nan")
    return k * k


def mu_star(m1: float, m2: float, m3: float) -> float:
    """The lepton-sum evaluation scale.  For the observed leptons this
    equals m_e + m_mu + m_tau = 1.883 GeV."""
    return m1 + m2 + m3


def soddy_pipeline(m1: float, m2: float, m3: float) -> dict:
    """Compute all Soddy quantities at once.  Returns a dict."""
    return {
        "e_1": e1(m1, m2, m3),
        "e_2": e2(m1, m2, m3),
        "k4_minus": outer_soddy_curvature(m1, m2, m3),
        "k4_plus": inner_soddy_curvature(m1, m2, m3),
        "F2": outer_soddy_sq(m1, m2, m3),
        "mu_star": mu_star(m1, m2, m3),
    }


if __name__ == "__main__":
    # Self-test: observed leptons should reproduce paper's F^2 = 95.113 MeV
    m_e = 0.51099895e-3       # GeV
    m_mu = 105.6583755e-3     # GeV
    m_tau = 1.77693           # GeV
    p = soddy_pipeline(m_e, m_mu, m_tau)
    print(f"Observed lepton Soddy pipeline:")
    print(f"  F^2     = {p['F2']*1000:.4f} MeV  (paper: 95.113 MeV)")
    print(f"  k4_minus = {p['k4_minus']*1000**0.5:.4f} MeV^(1/2)  (paper: 9.7526)")
    print(f"  mu_star = {p['mu_star']*1000:.2f} MeV  (paper: 1883.10)")
