# PRE_REGISTRATION.md — v4 MC for Brilliant 2026 Soddy paper

**Pre-registered on:** 2026-04-19
**Author:** A. M. Brilliant
**Paper under test:** "A Monte Carlo null-model test of an outer-Soddy completion of the Koide lepton triple" (Brilliant 2026, bibkey `Brilliant2026strange`). *Title amended post pre-registration; see amendment note below. The pre-registration-frozen state at commit `5ab2b88` preserves the original title "A Descartes–Soddy completion of the Koide lepton triple".*
**Replaces:** MC reported in §5.4–5.6 of the paper
**Repository:** https://github.com/AndBrilliant/fermion_mass_mc

This document specifies all analysis choices BEFORE running the MC, so
that a reviewer can distinguish pre-specified from post-hoc decisions.

---

## 1. Scope of the test

A single numerical observation is being tested:

> $\mathcal{F}^2 = 95.1134$ MeV (outer Soddy curvature squared of the
> Koide lepton triple) vs $m_s^{\overline{\text{MS}}}(\mu_\star)$ at
> $\mu_\star = m_e + m_\mu + m_\tau \approx 1.883$ GeV. Observed residual
> $+0.038$ MeV against FLAG 2024 uncertainty $\pm 0.692$ MeV.

(Residual value per the 2026-04-20 amendment below. The pre-registration
value was $+0.042$ MeV; the updated value reflects corrected σ
propagation identified in the adversarial audit. The observation itself
is unchanged; only its displayed σ-normalization is corrected.)

The MC does NOT test heavy-Koide (c,b,t), the golden quark-lepton mass
relation, or any other downstream claim. Those are either out of scope
for the Soddy paper or handled by the companion v3 LEE paper.

---

## 2. Four test strands

### Null A — random-spectrum, Koide-conditioned

**Question:** Would a random Koide-compatible universe, drawn from a
literature-motivated prior, produce F² within 1σ of m_s(μ_⋆)?

**Sampling:** direct sampling on the Koide manifold via the analytic
quadratic (engine/koide_manifold.py). No rejection, no approximation.

**Prior shapes (pre-registered 4-shape grid, NOT a range grid):**

| ID | Shape               | Citation                              | Role                 |
|----|---------------------|---------------------------------------|----------------------|
| A1 | Log-uniform         | Donoghue 1998, hep-ph/9712333         | **Primary**          |
| A2 | Log-normal          | Hall-Salem-Watari 2007, 0707.3446     | Shape robustness     |
| A3 | Yukawa-anarchy SVD  | Hall-Murayama-Weiner 2000 + dGM 2003  | Shape robustness     |
| A4 | Linear-uniform      | paper §5.5 "D" variant                | Stress test          |

**Hit criterion:** `|F² − m_s^MSbar(μ_⋆)| < σ_m_s^MSbar(μ_⋆)`.
**Sample size:** `N_SAMPLES_PER_PRIOR = 10_000` valid Koide triples.
**Cutoff:** μ_⋆ > 1 GeV (matches paper §5.4).
**Seeds:** `2026, 2027, 2028, 2029` for A1..A4.

### Null B — measurement-noise parametric bootstrap

**Question:** Is the observed residual distinguishable from noise
under PDG+FLAG measurement uncertainties on the actual inputs?
*(This null is missing from the paper.)*

**Inputs perturbed:** `m_e, m_mu, m_tau` (PDG 2024), `m_s(2 GeV)`
(FLAG 2024), `α_s(M_Z)` (PDG 2024). Truncated to positive values.

**Variants (pre-registered):**

| ID | Configuration                    | Role                    |
|----|----------------------------------|-------------------------|
| B1 | Uncorrelated, 1σ                 | **Baseline**            |
| B2 | Lattice-correlated ρ=0.5, 1σ     | Correlation robustness  |
| B3 | Uncorrelated, 0.5σ (tighter)     | σ-scale sensitivity     |
| B4 | Uncorrelated, 2σ (looser)        | σ-scale sensitivity     |

**Test statistics (both reported):**
- `P(perturbed residual ≤ observed +0.038 MeV)` [per 2026-04-20 amendment]
- Fraction producing a paper-style hit.

**Sample size:** `N_BOOTSTRAP = 2_000` per variant.
**Seeds:** `2126, 2127, 2128, 2129`.

### Test C — algorithm-space enumeration (paper Appendix)

**Question:** Is F² the only similar-complexity function of the Koide
triple that lands near m_s, or are there several?

**Procedure:** enumerate ~25 algebraic functions (engine/algorithm_variants.py);
evaluate each on observed leptons, count hits within 1σ of m_s(μ_⋆);
for each variant also compute Null A hit rate for comparison.

**Pre-registered catalog:** fixed at initial commit. Any additions
post-run must be reported as such.
**Seed:** 2226.

### Test D — temporal convergence

**Question:** *Has the residual behaved over time the way a real
regularity's would?*

As experimental precision has tightened across FLAG releases (2013 σ≈2.4
MeV → 2024 σ=0.68 MeV), a **real regularity** should produce a residual
that shrinks in MeV while staying roughly constant in units of σ (i.e.,
F² tracks the narrowing m_s band). A **coincidence at 2024 central
values** would produce a residual that stays stable in MeV as σ shrinks
around it, so residual/σ grows toward the value it has in the narrowest
band.

This is the third discriminator in the v3 LEE paper §6.3 and the
substance of the Soddy paper §1 closing paragraph ("the discriminatory
power to notice it has only recently arrived"). Test D makes it
quantitative.

**Inputs:** historical FLAG m_s(2 GeV) central values + σ for
FLAG 2013, 2016, 2019, 2021, 2024. Loaded from
`results/historical_masses.json` (compiled 2026-04-16, see v3
`results/temporal_history.json` for provenance).

**Test statistics (descriptive):**
1. Trajectory table: residual (MeV), σ (MeV), residual/σ, within-1σ? per epoch.
2. Chi-square against m_s^true = F² hypothesis. Reported for both
   naive N_epochs and conservative N_independent (treating 2019/2021/2024
   as one, which shares lattice inputs).
3. Inverse-variance weighted mean residual: under real-regularity,
   should be consistent with zero.
4. MC coincidence null (N=50,000, seed 2326): draw m_s^true uniformly
   from [60, 130] MeV (historical range of all determinations); for
   each epoch, sample measured values around m_s^true; count fraction
   of MC draws where F² lands within 1σ at ALL epochs. Observed: 5/5.
   Small fraction under null = observed trajectory is unusual for a
   coincidence.

**What Test D does NOT do:**
- It does not forecast future narrowing. The historical data is fixed.
- It does not separate "real-regularity" from "coincidence where
  m_s^true happens to equal F²." Those two are observationally
  indistinguishable from the trajectory alone. Test D characterizes
  whether the trajectory is consistent with the former.

---

## 3. Inputs (all locked at time of pre-registration)

### Lepton pole masses (PDG 2024)
- $m_e = 0.51099895$ MeV, σ = 0.00000015 MeV
- $m_\mu = 105.6583755$ MeV, σ = 0.0000023 MeV
- $m_\tau = 1776.93$ MeV, σ = 0.09 MeV

### Strange quark (FLAG 2024, N_f = 2+1+1)
- $m_s^{\overline{\text{MS}}}(2 \text{ GeV}) = 93.44$ MeV, σ = 0.68 MeV

### QCD coupling (PDG 2024)
- $α_s(M_Z) = 0.1180$, σ = 0.0009

### Running
- Four-loop MSbar mass RG with flavor-threshold matching via
  `rundec` (engine/running.py). The σ used for hit criteria
  throughout is the FLAG lattice input propagated through running
  (±0.692 MeV at μ_⋆); α_s and truncation contributions are
  separately enumerated in the main paper's error budget (§5.1).
- Default thresholds: $m_c(m_c) = 1.273$ GeV, $m_b(m_b) = 4.183$ GeV,
  $m_t = 173.0$ GeV.

### Historical m_s trajectory (Test D only)
See `results/historical_masses.json`. Five FLAG releases (2013, 2016,
2019, 2021, 2024). All m_s(2 GeV) values are MSbar Nf=2+1+1 where
available.
