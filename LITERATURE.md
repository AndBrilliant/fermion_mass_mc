# LITERATURE.md — Citations and motivations for v4

Every methodological choice in `DESIGN.md` is backed by one or more
verified references. This document maps each choice to its citation(s),
records the quote or argument being used, and explains relevance.

**Status:** citations verified against primary sources 2026-04-19.
See v4 research report `../LITERATURE_RESEARCH_REPORT.md` for the
underlying survey.

---

## A. Prior distribution: log-uniform / scale-invariant

The primary prior is log-uniform on masses. The supporting literature
establishes this as the **phenomenologically best-fitting** null, not as
a first-principles derivation — we frame it honestly on that basis.

### A1. Donoghue 1998 — the original "weight" paper

**Citation:** J. F. Donoghue, *"The weight for random quark masses,"*
Phys. Rev. D **57**, 5499–5508 (1998). arXiv:hep-ph/9712333.
DOI: 10.1103/PhysRevD.57.5499.

**What the paper actually establishes.** Donoghue argues that if SM
Yukawa couplings are drawn from an underlying multiverse distribution,
the observed pattern of quark masses is best fit by a weight
ρ(m) ∝ 1/m — equivalently uniform in log m, equivalently scale-invariant.
The abstract uses the hedged phrasing "roughly scale invariant." This is
a phenomenological fit, not a derivation from first principles.

**How v4 uses it.** Primary motivation for the log-uniform prior N1.
Framed honestly: "Donoghue (1998) found the observed quark mass spectrum
is best fit by a scale-invariant weight ρ(m) ∝ 1/m; we adopt this form
as our primary null prior."

### A2. Donoghue–Dutta–Ross 2006 — statistical goodness-of-fit

**Citation:** J. F. Donoghue, K. Dutta, A. Ross, *"Quark and lepton
masses and mixing in the landscape,"* Phys. Rev. D **73**, 113002 (2006).
arXiv:**hep-ph/0511219**. DOI: 10.1103/PhysRevD.73.113002.

*(Correction: an earlier draft of this file cited arXiv:0807.2330 and
PRD 78, 033003 (2008). Both are wrong. 0807.2330 is a computer-science
paper; PRD 78, 033003 is unrelated. The correct citation is above.)*

**What the paper actually establishes.** Monte Carlo sampling from
ρ(m) ∝ 1/m and comparison against the full observed SM mass pattern
(quarks, charged leptons, CKM elements, Jarlskog invariant). The
abstract states: "We use statistical techniques to show that the
observed masses appear to be representative of a scale invariant
distribution." Charged leptons and quark masses are consistent. The
neutrino sector is the weak point. This is the **direct methodological
precedent** for a MC tail-probability test under a log-uniform Yukawa
prior.

**How v4 uses it.** Cited as the methodological template for what we
are doing. We are extending their framework (which tested the overall
SM pattern) to test a specific predicted-mass claim (Koide + Soddy).

### A3. Donoghue–Dutta–Ross–Tegmark 2010 — log-uniform as explicit prior

**Citation:** J. F. Donoghue, K. Dutta, A. Ross, M. Tegmark, *"Likely
values of the Higgs vacuum expectation value,"* Phys. Rev. D **81**,
073003 (2010). arXiv:0903.1024.

**What the paper establishes.** Uses ρ(y) ∝ 1/y explicitly as a Bayesian
prior for Yukawa couplings. Figure caption states: "A scale invariant
weight corresponds to a uniform distribution on this [log] scale."

**How v4 uses it.** Direct precedent for treating log-uniform as a
prior distribution in a frequentist/Bayesian hybrid analysis.

### A4. Donoghue 2016 — recent synthesis

**Citation:** J. F. Donoghue, *"The multiverse and particle physics,"*
Annu. Rev. Nucl. Part. Sci. **66**, 1–21 (2016). arXiv:1601.05136.
DOI: 10.1146/annurev-nucl-102115-044644.

**Relevance.** Donoghue's own review, ~20 years after the original,
reaffirms the 1/m weight as the preferred phenomenological form. Not
superseded. Good one-stop citation for reviewer context.

---

## B. Alternative prior shapes (for sensitivity analysis)

The log-uniform prior is not the only defensible shape. Reporting across
multiple shapes is the strongest defense against a "circular prior"
objection.

### B1. Gaussian landscape — Hall, Salem, Watari 2007/2008

**Citations:**
- L. J. Hall, M. P. Salem, T. Watari, *"Statistical understanding of
  quark and lepton masses in Gaussian landscapes,"* Phys. Rev. D **76**,
  093001 (2007). arXiv:0707.3446.
- L. J. Hall, M. P. Salem, T. Watari, *"Quark and lepton masses from
  Gaussian landscapes,"* Phys. Rev. Lett. **100**, 141801 (2008).
  arXiv:0710.4262.

**What these papers establish.** Yukawa couplings arise from overlap
integrals of Gaussian wavefunctions in extra dimensions. The resulting
mass distribution is **log-normal-like** (not log-uniform) and reproduces
the observed flavor hierarchy. This is the principal alternative to
Donoghue's 1/m in the landscape literature.

**How v4 uses it.** Motivation for the sensitivity prior N2 (log-normal).
If our tail probability is stable between N1 (log-uniform) and N2
(log-normal), the result is robust to the prior shape.

### B2. Anarchy — Hall, Murayama, Weiner 2000; de Gouvêa, Murayama

**Citations:**
- L. J. Hall, H. Murayama, N. Weiner, *"Neutrino mass anarchy,"*
  Phys. Rev. Lett. **84**, 2572 (2000). arXiv:hep-ph/9911341.
- N. Haba, H. Murayama, *"Anarchy and hierarchy,"* Phys. Rev. D **63**,
  053010 (2001). arXiv:hep-ph/0009174.
- A. de Gouvêa, H. Murayama, *"A statistical test of anarchy,"*
  Phys. Lett. B **573**, 94 (2003). arXiv:hep-ph/0301050.
- A. de Gouvêa, H. Murayama, *"Neutrino mixing anarchy: alive and
  kicking,"* Phys. Lett. B **747**, 479 (2015). arXiv:1204.1249.

**What these papers establish.** Yukawa matrix entries are drawn
uniformly as O(1) numbers under the Haar measure; mass eigenvalues are
then derived from SVD. The de Gouvêa–Murayama papers are the **direct
methodological ancestors** of our test: they use a random-prior MC and
a KS-style statistic to compute an empirical tail probability
(44%–64% for anarchy's consistency with observed mixing).

**How v4 uses it.** Cited as the direct methodological precedent. Also
motivates an optional sensitivity prior N4 (Yukawa-matrix entries drawn
from a Gaussian-O(1) distribution, eigenvalues interpreted as masses).

### B3. Froggatt–Nielsen — the symmetry-based alternative

**Citation:** C. D. Froggatt, H. B. Nielsen, Nucl. Phys. B **147**, 277
(1979). DOI: 10.1016/0550-3213(79)90316-X.

**Relevance.** Power-law-suppressed Yukawas under a U(1) flavor
symmetry. Different prior philosophy from both landscape approaches.
Mentioned as a known alternative; not implemented as a separate
sampler in v4 for simplicity.

---

## C. Anthropic constraints (restricted scope)

**Critical correction from our literature review.** Anthropic bounds in
the published literature exist only for m_u, m_d, and m_e. No anthropic
bound exists on m_s, m_c, m_b, m_t, m_μ, or m_τ. This is stated
explicitly in the primary source (ABDS 1998). The v4 "anthropic" prior
must therefore be restricted to the first generation or clearly
relabeled.

### C1. Agrawal–Barr–Donoghue–Seckel 1998 — the short paper

**Citation:** V. Agrawal, S. M. Barr, J. F. Donoghue, D. Seckel,
*"Anthropic considerations in multiple-domain theories and the scale
of electroweak symmetry breaking,"* Phys. Rev. Lett. **80**, 1822
(1998). arXiv:**hep-ph/9801253**. DOI: 10.1103/PhysRevLett.80.1822.

### C2. Agrawal–Barr–Donoghue–Seckel 1998 — the long paper

**Citation:** V. Agrawal, S. M. Barr, J. F. Donoghue, D. Seckel, *"The
anthropic principle and the mass scale of the Standard Model,"*
Phys. Rev. D **57**, 5480 (1998). arXiv:**hep-ph/9707380**.
DOI: 10.1103/PhysRevD.57.5480.

*(Correction: an earlier draft of this file may have transposed these
two arXiv IDs between the short PRL and the long PRD.)*

**What these papers establish.** Direct quote from the long paper:
"The masses of the quarks and leptons of the second and third
generations have very little impact on the particles and reactions
which occur naturally in our universe. Therefore anthropic arguments
would not place constraints on their masses. In an ensemble of
anthropically allowed universes, these masses could be randomly
distributed." The first generation bounds establish (m_u + m_d)_max
roughly 2× observed, m_d − m_u bounded above by the octet-decuplet
splitting, and m_e ≲ few MeV from BBN.

**How v4 uses it.** Anthropic restriction applied **only** to
(m_u, m_d, m_e). For the other six fermions, the anthropic-labeled
prior uses the same log-uniform range as N1 with explicit acknowledgment
in the paper that no bound exists.

### C3. Damour–Donoghue 2008 — tightest modern first-generation bound

**Citation:** T. Damour, J. F. Donoghue, *"Constraints on the
variability of quark masses from nuclear binding,"* Phys. Rev. D **78**,
014014 (2008). arXiv:0712.2968. DOI: 10.1103/PhysRevD.78.014014.

**What the paper establishes.** Nuclear-EFT-tightened bounds:
(m_u + m_d)/(m_u + m_d)_phys < 1.64 at 95% CL (so m_u + m_d ≲ 18 MeV);
hydrogen stability requires m_d − m_u − 1.67 m_e ≥ 0.83 MeV; combined
Higgs-vev bound 0.39 ≤ v/v_phys ≤ 1.64. Explicitly: "m_u and m_e have
no lower anthropic bounds, while m_d is constrained to be non-zero."

**How v4 uses it.** Quantitative source for the (u, d, e) tight window
in sensitivity prior N3. Light quark masses drawn from a narrow
log-uniform window consistent with m_u + m_d < 18 MeV.

### C4. Other anthropic references (context)

- M. Tegmark, A. Aguirre, M. J. Rees, F. Wilczek, Phys. Rev. D **73**,
  023505 (2006). arXiv:astro-ph/0511774. *Concurs with ABDS on limited
  scope of anthropic constraints.*
- R. Harnik, G. D. Kribs, G. Perez, Phys. Rev. D **74**, 035006 (2006).
  arXiv:hep-ph/0604027. *Counter-example: habitable "weakless" universe
  with Yukawas rescaled.*
- R. L. Jaffe, A. Jenkins, I. Kimchi, Phys. Rev. D **79**, 065014 (2009).
  arXiv:0809.1647. *Light-quark environmental impact; parallel analysis
  to Damour–Donoghue.*

---

## D. Test statistic and methodology

### D1. Test statistic rationale

Max-sigma-residual is the Kolmogorov-style aggregation. Supplementary:
report the **full histogram** of Q (or the test statistic) under each
prior, not only a tail probability. Koide's Q is algebraically bounded
[1/3, 1), so the null distribution is narrow and non-Gaussian; summary
p-values alone hide the shape.

### D2. Cousins — frequentist framing

**Citation:** R. D. Cousins, *"Lectures on statistics in theory: prelude
to statistics in practice,"* arXiv:1807.05996 [physics.data-an], v4
(2024). *(Not peer-reviewed; widely cited.)*

**Relevance.** Foundational framing for frequentist p-values, coverage,
and the choice between frequentist and Bayesian reporting. Our tail
probabilities are frequentist p-values under the null; we do not compute
Bayes factors (though Beaujean et al. would be the reference if we did).

### D3. Cousins–Mayo correspondence 2026 on severity

**Cousins, personal communication to A. Brilliant, April 10, 2026.**

Per Cousins: severity is post-data power, computed under an alternative
pdf. It is **not** the same as 1 − p. Our quantities are empirical tail
probabilities under H₀, NOT severity. We use the language "empirical
tail probability" and "look-elsewhere-corrected p-value" exclusively.

---

## E. Look-elsewhere effect (for future scale-scan work; not Test 1)

### E1. Gross–Vitells 2010 — canonical HEP LEE

**Citation:** E. Gross, O. Vitells, *"Trial factors for the look
elsewhere effect in high energy physics,"* Eur. Phys. J. C **70**, 525
(2010). arXiv:1005.1891. DOI: 10.1140/epjc/s10052-010-1470-8.

**Relevance.** The canonical reference for 1D scan trial factors in HEP.
Applies to our μ-scan for heavy-Koide with caveats (Koide's Q is not a
χ² random field; Wilks' theorem may not hold; validate by MC).

### E2. Davies 1987 — foundational statistical result

**Citation:** R. B. Davies, Biometrika **74**, 33 (1987).
DOI: 10.1093/biomet/74.1.33.

**Relevance.** The bound Gross–Vitells operationalize. Cite as
methodological foundation.

### E3. Algeri et al. 2016 — modern refinement

**Citations:**
- S. Algeri, J. Conrad, D. A. van Dyk, B. Anderson, arXiv:1602.03765
  (2016).
- S. Algeri, D. A. van Dyk, Statistica Sinica **30** (2020),
  arXiv:1701.06820.

**Relevance.** Modern comparison of Gross–Vitells against alternatives;
identifies regimes where Gross–Vitells is biased. Cite when using
parametric bootstrap as the primary method.

### E4. Lyons 2008, 2013 — narrative

- L. Lyons, Ann. Appl. Stat. **2**, 887 (2008). arXiv:0811.1663.
- L. Lyons, arXiv:1310.1284 (2013). *"Discovering the significance of 5σ."*

**Relevance.** Narrative justification for why trial factors are needed
and for HEP's 5σ threshold convention.

---

## F. Koide formula — context (no prior MC null test exists)

**Key finding:** to the best of a thorough literature survey, no prior
HEP paper has implemented a MC null test for the Koide Q = 2/3
observation with a random-mass prior. We are extending the methodology
of the anarchy program (de Gouvêa–Murayama) and the landscape program
(Donoghue–Dutta–Ross) to this observable.

### F1. Original Koide papers

- Y. Koide, Lett. Nuovo Cim. **34**, 201 (1982).
- Y. Koide, Phys. Lett. B **120**, 161 (1983).
- Y. Koide, Phys. Rev. D **28**, 252 (1983).

### F2. Geometric interpretation

- R. Foot, arXiv:hep-ph/9402242. *Q as cos² of the angle between
  (√m_e, √m_μ, √m_τ) and (1,1,1).*
- A. Rivero, A. Gsponer, arXiv:hep-ph/0505220. *"The strange formula
  of Dr. Koide" — historical review; notes Q = 2/3 is the midpoint of
  [1/3, 1]. This is a deterministic fact about the range, not a
  probability, and must not be confused with a null tail probability.*

### F3. Quark extensions and running analyses

- J.-M. Gérard, F. Goffinet, M. Herquet, arXiv:hep-ph/0510289.
- A. Kartavtsev, arXiv:1111.0480.
- N. Li, B.-Q. Ma et al., Eur. Phys. J. C **76**, 167 (2016).
  arXiv:1601.00805.
- Y. Sumino, JHEP **0905**, 075 (2009); Phys. Lett. B **671**, 477
  (2009). *Gauged family symmetry explanation.*

None of these is a statistical null test under a random-mass prior.

---

## G. Internal references

- `../mass_as_curvature_mc_v3/` — predecessor. Tested the heavy-Koide
  claim (p ≈ 0.002 for (c,b,t) at m_τ). Inherited/cited by v4.
- `../cousins_simulation/` — corpus of Cousins's papers summarized for
  context; April 10 2026 email on severity semantics.
- `../main.tex` — the manuscript being tested.
- `LITERATURE_RESEARCH_REPORT.md` — the underlying verified survey
  from which these citations were drawn.
