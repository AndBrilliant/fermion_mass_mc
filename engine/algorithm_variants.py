"""
engine.algorithm_variants — enumeration of similar-complexity functions
of the Koide lepton triple.

Purpose: defend the paper's §3 claim that F = e_1 - 2 sqrt(e_2) is not
selected from a catalog of polynomials.  We enumerate ~25 algebraically
simple functions of (m_1, m_2, m_3), compute each on the observed
charged-lepton triple, and check how many land within 1 sigma of
m_s^MSbar(mu_star).

The list is organised by complexity tier.  Functions of dimension
[mass^(1/2)] are squared before comparison.  Functions of dimension
[mass^k] for k != 1 are raised to the appropriate power or omitted.

If only the outer Soddy curvature lands, the paper's §3 argument is
strengthened.  If multiple functions land, that is relevant
information.
"""
import math
from typing import Callable, List, Tuple


def _sqrt_safe(x: float) -> float:
    return math.sqrt(x) if x >= 0 else float("nan")


# Each entry: (name, callable(m1, m2, m3) -> value_in_GeV, dimension_note)
# All return a quantity of dimension [mass] = [GeV].
# For functions naturally of dimension [mass^(1/2)] we square; for
# dimension [mass^2] we take sqrt; etc.  Dimensional consistency is
# noted in `dimension_note`.
ALGORITHM_VARIANTS: List[Tuple[str, Callable, str]] = [

    # --- Soddy / Descartes roots (the paper's family) ----------------
    ("outer_soddy_sq",
     lambda m1, m2, m3: (_sqrt_safe(m1) + _sqrt_safe(m2) + _sqrt_safe(m3)
                         - 2*_sqrt_safe(_sqrt_safe(m1*m2) + _sqrt_safe(m2*m3)
                                        + _sqrt_safe(m1*m3)))**2,
     "F^2 = (k_4^-)^2, the paper's choice"),

    ("inner_soddy_sq",
     lambda m1, m2, m3: (_sqrt_safe(m1) + _sqrt_safe(m2) + _sqrt_safe(m3)
                         + 2*_sqrt_safe(_sqrt_safe(m1*m2) + _sqrt_safe(m2*m3)
                                        + _sqrt_safe(m1*m3)))**2,
     "(k_4^+)^2, inner Soddy"),

    # --- Elementary symmetric polynomials, dimensionally squared to mass
    ("e1_squared",
     lambda m1, m2, m3: (_sqrt_safe(m1) + _sqrt_safe(m2) + _sqrt_safe(m3))**2,
     "e_1^2 = p_2 + 2 e_2 = (3/2) p_2 under Koide"),

    ("e1_squared_over_3",
     lambda m1, m2, m3: (_sqrt_safe(m1) + _sqrt_safe(m2) + _sqrt_safe(m3))**2 / 3.0,
     "e_1^2 / 3"),

    ("e2_raw",
     lambda m1, m2, m3: _sqrt_safe(m1*m2) + _sqrt_safe(m2*m3) + _sqrt_safe(m1*m3),
     "e_2 (pairwise sqrt products)"),

    ("sqrt_e2_sq",
     lambda m1, m2, m3: _sqrt_safe(m1*m2) + _sqrt_safe(m2*m3) + _sqrt_safe(m1*m3),
     "e_2 again, under different naming"),

    # --- Power mean functions ---------------------------------------
    ("arithmetic_mean",
     lambda m1, m2, m3: (m1 + m2 + m3) / 3.0,
     "(m1+m2+m3)/3"),

    ("geometric_mean",
     lambda m1, m2, m3: (m1 * m2 * m3) ** (1.0/3.0) if m1*m2*m3 > 0 else float("nan"),
     "(m1 m2 m3)^(1/3)"),

    ("harmonic_mean",
     lambda m1, m2, m3: 3.0 / (1.0/m1 + 1.0/m2 + 1.0/m3)
                        if min(m1,m2,m3) > 0 else float("nan"),
     "harmonic mean"),

    ("sum_all",
     lambda m1, m2, m3: m1 + m2 + m3,
     "p_2 = sum m_i = mu_star"),

    # --- Koide quadratic "second root" squared ----------------------
    # u_2 = 2(u_1+u_3) + sqrt(3)*sqrt(u_1^2 + 4 u_1 u_3 + u_3^2)
    ("koide_plus_root_sq",
     lambda m1, m2, m3: (
         2*(_sqrt_safe(m1) + _sqrt_safe(m3))
         + math.sqrt(3) * _sqrt_safe(m1 + 4*_sqrt_safe(m1*m3) + m3)
     )**2,
     "square of Koide quadratic '+' root using m1, m3"),

    # --- Functions of the heaviest two -------------------------------
    ("sqrt_m2_m3_prod",
     lambda m1, m2, m3: _sqrt_safe(m2 * m3),
     "sqrt(m_mu * m_tau)"),

    ("m3_minus_m2",
     lambda m1, m2, m3: m3 - m2,
     "m_tau - m_mu"),

    ("m3_over_m2_times_m2",
     lambda m1, m2, m3: m3 * m2 / (m2 + m3),
     "reduced-mass-like combo of m_mu and m_tau"),

    # --- Functions of lightest + heaviest ----------------------------
    ("sqrt_m1_m3_prod",
     lambda m1, m2, m3: _sqrt_safe(m1 * m3),
     "sqrt(m_e * m_tau)"),

    ("m1_plus_m2",
     lambda m1, m2, m3: m1 + m2,
     "m_e + m_mu"),

    # --- Double-sqrt / quartic combinations --------------------------
    ("sqrt_sum_sqrt",
     lambda m1, m2, m3: (_sqrt_safe(m1) + _sqrt_safe(m2) + _sqrt_safe(m3))**2 / 6.0,
     "e_1^2 / 6 = e_2 under Koide"),

    ("fourth_root_prod_sq",
     lambda m1, m2, m3: ((m1 * m2 * m3) ** (1.0/4.0))**2 if m1*m2*m3 > 0 else float("nan"),
     "sqrt((m1 m2 m3)^(1/2)) = (m1 m2 m3)^(1/2)"),

    # --- Ratios and differences of elementary symmetric pieces --------
    ("e1_minus_sqrt_p2",
     lambda m1, m2, m3: ((_sqrt_safe(m1) + _sqrt_safe(m2) + _sqrt_safe(m3))
                         - _sqrt_safe(m1+m2+m3))**2,
     "F surrogate = (e_1 - sqrt(p_2))^2, equals outer Soddy under Koide"),

    # --- Simple differences ------------------------------------------
    ("e1_sq_minus_p2",
     lambda m1, m2, m3: (_sqrt_safe(m1) + _sqrt_safe(m2) + _sqrt_safe(m3))**2
                        - (m1 + m2 + m3),
     "e_1^2 - p_2 = 2 e_2"),

    # --- Determinant-like -------------------------------------------
    ("m_mid_squared",
     lambda m1, m2, m3: m2**2 / (m1 + m3) if (m1 + m3) > 0 else float("nan"),
     "m_mu^2 / (m_e + m_tau)"),

    ("m_heavy_sqrt_times_light",
     lambda m1, m2, m3: _sqrt_safe(m3) * _sqrt_safe(m1),
     "same as sqrt_m1_m3_prod, duplicate for completeness"),

    # --- Functions raised to different powers (dimensional test) -----
    ("p2_minus_e2",
     lambda m1, m2, m3: (m1 + m2 + m3)
                        - (_sqrt_safe(m1*m2) + _sqrt_safe(m2*m3) + _sqrt_safe(m1*m3)),
     "p_2 - e_2"),

    ("two_e2_minus_p2",
     lambda m1, m2, m3: 2 * (_sqrt_safe(m1*m2) + _sqrt_safe(m2*m3) + _sqrt_safe(m1*m3))
                        - (m1 + m2 + m3),
     "2 e_2 - p_2 = 0 under Koide (control: should give ~0 for leptons)"),
]


def evaluate_all_variants(m1: float, m2: float, m3: float) -> dict:
    """Evaluate every algorithm variant on (m1, m2, m3).  Returns
    {name: value_in_GeV}."""
    out = {}
    for name, func, _note in ALGORITHM_VARIANTS:
        try:
            v = func(m1, m2, m3)
        except (ValueError, ZeroDivisionError, OverflowError):
            v = float("nan")
        out[name] = v
    return out


def count_hits_on_observed(m1: float, m2: float, m3: float,
                           m_s_target: float, sigma_m_s: float) -> dict:
    """For the observed leptons, count how many variants land within
    1 sigma of the strange-quark mass."""
    vals = evaluate_all_variants(m1, m2, m3)
    hits = {}
    misses = {}
    for name, v in vals.items():
        if not math.isfinite(v):
            continue
        if abs(v - m_s_target) < sigma_m_s:
            hits[name] = (v, (v - m_s_target) / sigma_m_s)
        else:
            misses[name] = (v, (v - m_s_target) / sigma_m_s)
    return {
        "n_variants": len(vals),
        "n_finite": len([v for v in vals.values() if math.isfinite(v)]),
        "n_hits": len(hits),
        "hits": hits,
        "closest_miss": min(misses.items(),
                            key=lambda kv: abs(kv[1][1]),
                            default=None),
    }


if __name__ == "__main__":
    # Self-test: evaluate all variants on observed leptons.
    m_e, m_mu, m_tau = 0.51099895e-3, 105.6583755e-3, 1.77693
    # m_s(mu_star = 1.883 GeV) from paper
    m_s, sig = 95.07e-3, 0.69e-3
    result = count_hits_on_observed(m_e, m_mu, m_tau, m_s, sig)
    print(f"Algorithm-space check on observed leptons:")
    print(f"  {result['n_variants']} variants, {result['n_finite']} finite")
    print(f"  {result['n_hits']} variants hit within 1 sigma of m_s")
    for name, (v, nsig) in result["hits"].items():
        print(f"    {name}: {v*1000:.3f} MeV ({nsig:+.2f} sigma)")
