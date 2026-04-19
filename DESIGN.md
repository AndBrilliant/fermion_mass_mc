# DESIGN.md — v4 Monte Carlo

*Companion to manuscript `main.tex` (Soddy formula / Descartes circle theorem
applied to the SM mass spectrum).*

## 1. The statistical problem

The manuscript makes a pipeline claim: from two measured inputs (m_e, m_μ) and
one mathematical condition (Koide's Q = 2/3), a fixed algorithm — Koide
quadratic, Soddy outer root with curvature identification k = √m, reflection
through the strange-quark mass — produces three predictions (m_τ, m_s, m_d)
that agree with PDG/FLAG measurements to ≤ 0.44σ each. A separate
observational claim extends this to (m_u, m_c, m_b, m_t) via the heavy-quark
Koide relation at μ = m_τ and the m_u ↔ m_b reflection.

The null hypothesis under test is:

> **H₀**: the observed agreement is a coincidence arising from the joint
> distribution of random fermion masses, and the algorithm itself has no
> predictive content beyond what scale-invariance would deliver.

Under H₀, if we replace the SM masses with random draws from a distribution
matching the statistical properties of observed Yukawa couplings, we should
recover the observed level of agreement with non-negligible frequency.

The alternative hypothesis — that the framework reflects structure — predicts
the observed level of joint agreement should be rare under H₀.

## 2. Test statistic

For a set of predictions {pred_i} with measured values {meas_i} and
measurement uncertainties {σ_i}, the test statistic is

        T  =  max_i |pred_i − meas_i| / σ_i.

We choose max-residual over a sum-squared χ² statistic for three reasons:

1. **It matches the paper's falsifiability criterion.** The manuscript
   (§ Falsifiability) asserts the framework is falsified if any single
   measurement moves outside its 3σ window. The worst prediction is the one
   that controls the claim.

2. **It is interpretable as a single number.** "Fraction of pseudo-universes
   in which all predictions simultaneously achieve this level of agreement"
   is a defensible English-language statement.

3. **It is conservative.** A pattern in which many predictions are marginally
   close but one is far would be caught by max, but might be averaged away
   by χ². We prefer to report the more stringent number.

## 3. Prior (Null model)

**Primary prior — N1 (Donoghue scale-invariant):**
log-uniform on [1 MeV, 200 GeV], masses drawn independently then sorted
ascending. Motivated by Donoghue's argument (Phys. Rev. D 57, 5499, 1998;
extended in Donoghue-Dutta-Ross 2008, Phys. Rev. D 78, 033003) that observed
Yukawa couplings are statistically consistent with a scale-invariant random
distribution.

**Sensitivity prior — N2 (Donoghue wider):**
log-uniform on [0.1 MeV, 2 TeV]. Tests insensitivity to range. One log-decade
buffer on each side of N1.

**Sensitivity prior — N3 (Anthropic-restricted):**
log-uniform on [1 MeV, 10 GeV]. Range motivated by the anthropic bounds of
Agrawal-Barr-Donoghue-Seckel (Phys. Rev. Lett. 80, 1822, 1998): masses
compatible with atomic/nuclear structure admitting chemistry.

We report the primary p-value under N1. If the result changes by more than a
factor of 2 between N1, N2, N3, we report "prior-sensitive"; otherwise we
report the least-favorable value as the conservative number.

## 4. The two tests

### Test 1 — Soddy pipeline alone

**Null procedure:**

1. Draw (m_1, m_2, m_3, m_4, m_5) from the prior, sorted ascending.
   Interpret the three smallest as lepton-like, the next two as strange-like
   and down-like.
2. Compute the Koide condition Q(m_1, m_2, m_3); condition on |Q − 2/3| <
   tolerance (i.e., accept only pseudo-universes where the lepton triple lies
   on the Koide manifold). This ensures the null is comparable to the
   observed case, where the leptons *do* satisfy Koide.
3. Apply the algorithm:
   - α_K = √(3/2) − 1 ≈ 0.22474
   - μ_★ = m_1 + m_2 + m_3
   - m_s^pred = α_K² · μ_★
   - m_d^pred = α_K⁴ · μ_★
4. Compute residuals against the drawn m_4 and m_5:
   - r_s = |m_s^pred − m_4| / (0.01 · m_4)    [FLAG-scale σ ≈ 1%]
   - r_d = |m_d^pred − m_5| / (0.02 · m_5)    [FLAG-scale σ ≈ 2%]
5. Test statistic T = max(r_s, r_d).
6. Observed T from the SM: max(0.08, 0.32) = 0.32.
7. Empirical p = fraction of retained pseudo-universes with T ≤ 0.32.

**Note on the Koide conditioning.** The algorithm's predictions depend on
the lepton triple being Koide-compatible. Conditioning on Koide in the null
is essential: it compares the structural claim (given Koide, Soddy predicts
strange) against a null that respects the same constraint. Without this
conditioning, we would be conflating "is the Soddy formula special" with
"is Koide itself special" — and Koide's own significance is a separate,
well-studied question.

### Test 2 — Full joint stack

**Null procedure:** as Test 1, but extended to all nine fermion masses.

1. Draw (m_1, ..., m_9) from prior, sorted.
2. Condition on Koide-compatibility of (m_1, m_2, m_3).
3. Condition on heavy-Koide-compatibility: some triple among (m_4..m_9) and
   some scale μ satisfies |Q − 2/3| < tolerance. (This is a conditioning
   on the *existence* of the heavy-Koide observation, not a discovery of
   it — the discovery question is addressed by v3.)
4. Apply the full algorithm: predict (m_τ, m_s, m_d, m_u, plus the heavy
   Koide best-residual).
5. Test statistic T = max over all predictions of |pred − drawn|/σ.
6. Observed T ≈ 0.63σ (driven by m_u before QED corrections) or 0.44σ
   (excluding m_u). Report both.

Test 2 is more demanding than Test 1 and requires CRunDec for the heavy-Koide
step. It reuses v3's engine where possible.

## 5. Sample size and standard errors

N = 10,000 per prior per test. This gives a binomial standard error of
√(p(1−p)/N) ≈ 0.001 at p = 0.01, so the p-values are reliable to two
significant figures at the 99%-significant level. We match v3's sample size
for continuity.

## 6. Decision thresholds (pre-registered)

- **Strong claim** — "the joint pattern is improbable under standard nulls":
  requires p < 0.01 under N1 for Test 2, and p < 0.05 for Test 1.

- **Weak claim** — "the pattern is unusual but not decisively so":
  requires p < 0.05 under N1 for Test 2.

- **Null result** — "pattern is consistent with coincidence":
  p > 0.1 under N1 for Test 2.

We commit to reporting all outcomes, including if the paper's claim is
weakened. This is the core pre-registration discipline; see
`PRE_REGISTRATION.md` for the full decision tree.

## 7. What we do NOT test in v4 (and why)

- **Algorithm-space LEE** (the set of alternative-algorithms-of-similar-complexity
  test): deferred. Valuable but requires separate pre-registration of the
  algorithm family, and the user chose to keep v4 focused.

- **Gross-Vitells scale-scan LEE** on μ = m_τ: v3 covered this empirically.
  A formal trial-factor correction via Gross-Vitells (arXiv:1005.1891) would
  strengthen the claim but is a separate analysis.

- **Severity** (in the Mayo sense): our tail probabilities are p-values
  computed under the null, not severity (which is attained power under an
  alternative). Per Cousins's April 2026 correspondence, we use the language
  "empirical tail probability" and "look-elsewhere-corrected p-value"
  exclusively and do not claim severity.

## 8. What qualifies as success for this MC

The MC is methodologically successful if:

1. The prior is defensibly motivated (Donoghue citation present).
2. All choices are pre-registered and seeds are frozen.
3. The result is reported honestly, whatever it is.
4. The full analysis is reproducible from the JSON outputs.

This is independent of whether the paper's substantive claim is supported.
A "null result" here is still a successful MC.
