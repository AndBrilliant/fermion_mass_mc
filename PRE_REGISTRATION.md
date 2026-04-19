# PRE_REGISTRATION.md â v4 MC for Brilliant 2026 Soddy paper

**Pre-registered on:** 2026-04-19
**Author:** A. M. Brilliant
**Paper under test:** "A Descartes–Soddy completion of the Koide lepton triple" (Brilliant 2026, bibkey `Brilliant2026strange`)
**Replaces:** MC reported in Â§5.4â5.6 of the paper
**Repository:** https://github.com/AndBrilliant/fermion_mass_mc

This document specifies all analysis choices BEFORE running the MC, so
that a reviewer can distinguish pre-specified from post-hoc decisions.

---

## 1. Scope of the test

A single numerical observation is being tested:

> $\mathcal{F}^2 = 95.1134$ MeV (outer Soddy curvature squared of the
> Koide lepton triple) vs $m_s^{\overline{\text{MS}}}(\mu_\star)$ at
> $\mu_\star = m_e + m_\mu + m_\tau \approx 1.883$ GeV. Observed residual
> $+0.042$ MeV against FLAG 2024 uncertainty $\pm 0.69$ MeV.

The MC does NOT test heavy-Koide (c,b,t), the golden quark-lepton mass
relation, or any other downstream claim. Those are either out of scope
for the Soddy paper or handled by the companion v3 LEE paper.

---

## 2. Four test strands

### Null A â random-spectrum, Koide-conditioned

**Question:** Would a random Koide-compatible universe, drawn from a
literature-motivated prior, produce FÂ² within 1Ï of m_s(Î¼_â)?

**Sampling:** direct sampling on the Koide manifold via the analytic
quadratic (engine/koide_manifold.py). No rejection, no approximation.

**Prior shapes (pre-registered 4-shape grid, NOT a range grid):**

| ID | Shape               | Citation                              | Role                 |
|----|---------------------|---------------------------------------|----------------------|
| A1 | Log-uniform         | Donoghue 1998, hep-ph/9712333         | **Primary**          |
| A2 | Log-normal          | Hall-Salem-Watari 2007, 0707.3446     | Shape robustness     |
| A3 | Yukawa-anarchy SVD  | Hall-Murayama-Weiner 2000 + dGM 2003  | Shape robustness     |
| A4 | Linear-uniform      | paper Â§5.5 "D" variant                | Stress test          |

**Hit criterion:** `|FÂ² â m_s^MSbar(Î¼_â)| < Ï_m_s^MSbar(Î¼_â)`.
**Sample size:** `N_SAMPLES_PER_PRIOR = 10_000` valid Koide triples.
**Cutoff:** Î¼_â > 1 GeV (matches paper Â§5.4).
**Seeds:** `2026, 2027, 2028, 2029` for A1..A4.

### Null B â measurement-noise parametric bootstrap

**Question:** Is the observed residual distinguishable from noise
under PDG+FLAG measurement uncertainties on the actual inputs?
*(This null is missing from the paper.)*

**Inputs perturbed:** `m_e, m_mu, m_tau` (PDG 2024), `m_s(2 GeV)`
(FLAG 2024), `Î±_s(M_Z)` (PDG 2024). Truncated to positive values.

**Variants (pre-registered):**

| ID | Configuration                    | Role                    |
|----|----------------------------------|-------------------------|
| B1 | Uncorrelated, 1Ï                 | **Baseline**            |
| B2 | Lattice-correlated Ï=0.5, 1Ï     | Correlation robustness  |
| B3 | Uncorrelated, 0.5Ï (tighter)     | Ï-scale sensitivity     |
| B4 | Uncorrelated, 2Ï (looser)        | Ï-scale sensitivity     |

**Test statistics (both reported):**
- `P(perturbed residual â¤ observed +0.042 MeV)`
- Fraction producing a paper-style hit.

**Sample size:** `N_BOOTSTRAP = 2_000` per variant.
**Seeds:** `2126, 2127, 2128, 2129`.

### Test C â algorithm-space enumeration (paper Appendix)

**Question:** Is FÂ² the only similar-complexity function of the Koide
triple that lands near m_s, or are there several?

**Procedure:** enumerate ~25 algebraic functions (engine/algorithm_variants.py);
evaluate each on observed leptons, count hits within 1Ï of m_s(Î¼_â);
for each variant also compute Null A hit rate for comparison.

**Pre-registered catalog:** fixed at initial commit. Any additions
post-run must be reported as such.
**Seed:** 2226.

### Test D â temporal convergence

**Question:** *Has the residual behaved over time the way a real
regularity's would?*

As experimental precision has tightened across FLAG releases (2013 Ïâ2.4
MeV â 2024 Ï=0.68 MeV), a **real regularity** should produce a residual
that shrinks in MeV while staying roughly constant in units of Ï (i.e.,
FÂ² tracks the narrowing m_s band). A **coincidence at 2024 central
values** would produce a residual that stays stable in MeV as Ï shrinks
around it, so residual/Ï grows toward the value it has in the narrowest
band.

This is the third discriminator in the v3 LEE paper Â§6.3 and the
substance of the Soddy paper Â§1 closing paragraph ("the discriminatory
power to notice it has only recently arrived"). Test D makes it
quantitative.

**Inputs:** historical FLAG m_s(2 GeV) central values + Ï for
FLAG 2013, 2016, 2019, 2021, 2024. Loaded from
`results/historical_masses.json` (compiled 2026-04-16, see v3
`results/temporal_history.json` for provenance).

**Test statistics (descriptive):**
1. Trajectory table: residual (MeV), Ï (MeV), residual/Ï, within-1Ï? per epoch.
2. Chi-square against m_s^true = FÂ² hypothesis. Reported for both
   naive N_epochs and conservative N_independent (treating 2019/2021/2024
   as one, which shares lattice inputs).
3. Inverse-variance weighted mean residual: under real-regularity,
   should be consistent with zero.
4. MC coincidence null (N=50,000, seed 2326): draw m_s^true uniformly
   from [60, 130] MeV (historical range of all determinations); for
   each epoch, sample measured values around m_s^true; count fraction
   of MC draws where FÂ² lands within 1Ï at ALL epochs. Observed: 5/5.
   Small fraction under null = observed trajectory is unusual for a
   coincidence.

**What Test D does NOT do:**
- It does not forecast future narrowing. The historical data is fixed.
- It does not separate "real-regularity" from "coincidence where
  m_s^true happens to equal FÂ²." Those two are observationally
  indistinguishable from the trajectory alone. Test D characterizes
  whether the trajectory is consistent with the former.

---

## 3. Inputs (all locked at time of pre-registration)

### Lepton pole masses (PDG 2024)
- $m_e = 0.51099895$ MeV, Ï = 0.00000015 MeV
- $m_\mu = 105.6583755$ MeV, Ï = 0.0000023 MeV
- $m_\tau = 1776.93$ MeV, Ï = 0.09 MeV

### Strange quark (FLAG 2024, N_f = 2+1+1)
- $m_s^{\overline{\text{MS}}}(2 \text{ GeV}) = 93.44$ MeV, Ï = 0.68 MeV

### QCD coupling (PDG 2024)
- $Î±_s(M_Z) = 0.1180$, Ï = 0.0009

### Running
- One-loop MSbar mass RG with flavor-threshold matching
  (engine/running.py). Accurate to ~0.5% RMS vs four-loop CRunDec over
  [1, 100] GeV (v2 Appendix A).
- Default thresholds: $m_c(m_c) = 1.273$ GeV, $m_b(m_b) = 4.183$ GeV,
  $m_t = 173.0$ GeV.

### Historical m_s trajectory (Test D only)
See `results/historical_masses.json`. Five FLAG releases (2013, 2016,
2019, 2021, 2024). All m_s(2 GeV) values are MSbar Nf=2+1+1 where
available.

---

## 4. Terminology (Cousins-compliant)

All probabilities reported are **empirical tail probabilities under the
specified null**, equivalently **Koide-conditioned hit fractions**.
They are NOT severity assessments in the Mayo-Spanos sense (Cousins
2018, arXiv:1807.05996). We say "p-value" and "tail probability"
interchangeably; we do NOT say "severity" or "1âp."

---

## 5. What counts as a pre-registered "success"

**No hypothesis is being tested against an alternative here.** The MC's
output is descriptive: tail probabilities and histograms that
characterize the observation's behavior under each null.

Specific descriptive claims the paper will make on the basis of the MC:

1. **Null A robustness.** If hit fraction under A1 is in [0.1%, 1%]
   AND is stable within a factor of 5 across A1âA4, the paper reports
   "stable across four distinct prior shapes."
2. **Null B distinguishability.** P(perturbed residual â¤ observed)
   under B1 characterizes how unusually small the residual is
   relative to measurement noise. Descriptive quantity.
3. **Algorithm-space character.** If outer_soddy_sq is one of â¤3
   variants that hit on observed leptons AND has a low Null A hit
   rate, the paper reports "outer Soddy is one of [N] functions in
   the enumerated catalog."
4. **Temporal convergence.** If Test D MC yields P(all-5-epochs-hit
   under coincidence null) â¤ 5%, the paper reports "the observed
   trajectory is unusual for a coincidence at 2024 central values."

---

## 6. Amendments after pre-registration

Any amendments will be logged in a new "amendments" section of this
file with date and justification, NOT by editing existing entries.

**Amendments log:**
- 2026-04-19: Added Test D (temporal convergence) at user request.
  This pre-registration was written the same day no MC has run yet,
  so the addition is still pre-specified. Historical data file was
  compiled 2026-04-16 and is itself pre-specified.

---

## 7. Code and reproducibility

- GitHub: https://github.com/AndBrilliant/fermion_mass_mc
- Pre-registration commit: [to be recorded at initial push]
- Production run commit: [to be recorded after running tests]
- All seeds listed in Â§2 above.
- Reproduction:
  ```
  python -m tests.validate_reproducibility   # smoke test
  python -m tests.test_a_random_spectrum     # ~minutes
  python -m tests.test_b_bootstrap           # ~minute
  python -m tests.test_c_algorithm_space     # ~minute
  python -m tests.test_d_temporal_convergence  # ~seconds
  ```

---

## Appendix A â what's different from the paper's Â§5.4â5.6

| Aspect                    | Paper Â§5.4â5.6        | v4 MC                          |
|---------------------------|------------------------|--------------------------------|
| Null A priors             | 4 (range variations)  | 4 distinct SHAPES              |
| Null B                    | **not present**        | 4-variant bootstrap            |
| Algorithm-space check     | Not formalized         | 25-variant enumeration         |
| Temporal convergence MC   | Qualitative in Â§1      | Quantitative Test D            |
| Prior citations           | Not given              | Donoghue, HSW, HMW, dGM        |
| Test statistic histogram  | Not saved              | Saved to results/              |
| Reproducibility           | commit bb8b1ae         | pre-registered seed list       |

The v4 MC is a replacement. `results/` contains a JSON file per test
with all numbers needed to regenerate any table or figure in the
paper revision.
