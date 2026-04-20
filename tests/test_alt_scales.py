"""
tests.test_alt_scales — §5.3 alternative-scale prescriptions table.

Audit-requested artifact for the paper's §5.3 scale-prescription table.
Computes F^2 residual (in both MeV and sigma units) at each of eight
lepton-derived / reference scale prescriptions. Non-perturbative
prescriptions are labeled rather than computed.

Uses PDG 2024 lepton pole masses (same as rest of audit) and FLAG 2024
N_f = 2+1+1 strange-quark input:
    m_s(2 GeV) = 93.44 +/- 0.68 MeV
    alpha_s(M_Z) = 0.1180

Usage:
    python -m tests.test_alt_scales
"""
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from engine import algorithm as alg
from engine import running


# PDG 2024 charged-lepton pole masses (GeV)
M_E = 0.51099895e-3
M_MU = 105.6583755e-3
M_TAU = 1.77693

# FLAG 2024 N_f = 2+1+1
M_S_2GEV_GEV = 0.09344
SIGMA_M_S_2GEV_GEV = 0.00068

# Perturbative threshold used in paper §5.3 (below this, running is
# non-perturbative; we label rather than compute).
MU_PERTURBATIVE_GEV = 1.0


def prescription_list():
    """Eight prescriptions from paper Table §5.3 in paper order."""
    return [
        ("mu_star = m_e + m_mu + m_tau (baseline)",
         M_E + M_MU + M_TAU),
        ("m_mu + m_tau",
         M_MU + M_TAU),
        ("m_tau",
         M_TAU),
        ("2 GeV (conventional reference, not lepton-derived)",
         2.0),
        ("(m_e + m_mu + m_tau) / 3",
         (M_E + M_MU + M_TAU) / 3.0),
        ("sqrt(m_mu * m_tau)",
         math.sqrt(M_MU * M_TAU)),
        ("(m_e * m_mu * m_tau)^(1/3)",
         (M_E * M_MU * M_TAU) ** (1.0 / 3.0)),
        ("3 / (1/m_e + 1/m_mu + 1/m_tau) (harmonic)",
         3.0 / (1.0 / M_E + 1.0 / M_MU + 1.0 / M_TAU)),
    ]


def main():
    F2_MeV = alg.outer_soddy_sq(M_E, M_MU, M_TAU) * 1000.0
    print("=" * 70)
    print("Alternative scale prescriptions (paper §5.3)")
    print("=" * 70)
    print(f"F^2 (engine, PDG 2024 leptons):  {F2_MeV:.4f} MeV")
    print(f"m_s(2 GeV) input:                {M_S_2GEV_GEV*1000:.2f} +/- "
          f"{SIGMA_M_S_2GEV_GEV*1000:.2f} MeV")
    print(f"Perturbative threshold:          {MU_PERTURBATIVE_GEV:.1f} GeV")
    print()
    print(f"{'Prescription':<55}{'mu (MeV)':>12}{'Hit?':>18}")
    print("-" * 85)

    results = []
    for label, mu_gev in prescription_list():
        mu_mev = mu_gev * 1000.0
        if mu_gev < MU_PERTURBATIVE_GEV:
            results.append({
                "prescription": label,
                "mu_MeV": mu_mev,
                "perturbative": False,
                "hit_within_1sigma": None,
                "n_sigma": None,
                "comment": "non-perturbative",
            })
            print(f"{label:<55}{mu_mev:>12.2f}{'non-perturbative':>18}")
            continue
        try:
            m_s, sig_m_s = running.m_s_at(
                mu_gev, M_S_2GEV_GEV, SIGMA_M_S_2GEV_GEV
            )
        except Exception as exc:
            results.append({
                "prescription": label,
                "mu_MeV": mu_mev,
                "perturbative": True,
                "hit_within_1sigma": None,
                "n_sigma": None,
                "comment": f"running failed: {exc}",
            })
            print(f"{label:<55}{mu_mev:>12.2f}{'running failed':>18}")
            continue
        m_s_mev = m_s * 1000.0
        sig_mev = sig_m_s * 1000.0
        residual_mev = F2_MeV - m_s_mev
        n_sigma = residual_mev / sig_mev
        hit = abs(residual_mev) < sig_mev
        results.append({
            "prescription": label,
            "mu_MeV": mu_mev,
            "perturbative": True,
            "m_s_mu_MeV": m_s_mev,
            "sigma_m_s_MeV": sig_mev,
            "residual_MeV": residual_mev,
            "n_sigma": n_sigma,
            "hit_within_1sigma": hit,
            "comment": "",
        })
        tag = f"yes ({n_sigma:+.2f}sigma)" if hit else f"no ({n_sigma:+.2f}sigma)"
        print(f"{label:<55}{mu_mev:>12.2f}{tag:>18}")

    out_path = Path(__file__).resolve().parent.parent / "results" / "test_alt_scales.json"
    with open(out_path, "w") as f:
        json.dump({
            "F2_MeV": F2_MeV,
            "m_s_2GeV_MeV": M_S_2GEV_GEV * 1000.0,
            "sigma_m_s_2GeV_MeV": SIGMA_M_S_2GEV_GEV * 1000.0,
            "perturbative_threshold_GeV": MU_PERTURBATIVE_GEV,
            "prescriptions": results,
        }, f, indent=2, default=str)
    print(f"\nResults written to {out_path}")


if __name__ == "__main__":
    main()
