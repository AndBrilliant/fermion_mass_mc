"""
mass_as_curvature_mc_v4.engine

Stronger Monte Carlo for the Koide-Soddy observation.  Replaces the MC
in §5.4-5.6 of the Brilliant 2026 Soddy paper.  Architecture:

    Null A (random-spectrum): four distinct prior SHAPES
        A1 Donoghue log-uniform (primary)
        A2 Hall-Salem-Watari log-normal
        A3 Yukawa-anarchy SVD
        A4 Linear-uniform stress test

    Null B (measurement-noise bootstrap): four variants
        B1 Uncorrelated PDG+FLAG Gaussians (baseline)
        B2 Correlated (rho=0.5) among lattice inputs
        B3 Tighter (0.5x sigma)
        B4 Looser (2x sigma)

    Algorithm-space check (paper Appendix): ~25 similar-complexity
        functions of (m_1, m_2, m_3), hit-counted on observed leptons.

All p-values are empirical tail probabilities under specified nulls.
They are NOT severity assessments in the Mayo sense (Cousins 2018,
arXiv:1807.05996).  See PRE_REGISTRATION.md.
"""
from . import algorithm
from . import algorithm_variants
from . import koide_manifold
from . import priors
from . import running
from . import statistics

__all__ = [
    "algorithm",
    "algorithm_variants",
    "koide_manifold",
    "priors",
    "running",
    "statistics",
]
