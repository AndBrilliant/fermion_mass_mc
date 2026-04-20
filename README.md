# fermion_mass_mc

Pre-registered Monte Carlo and reproducibility artifacts for

> Brilliant, A. M., "A Descartes–Soddy completion of the Koide lepton triple,"
> Preprints.org (2026). [`Brilliant2026strange`]

This repository contains the code, pre-registration, and committed result
artifacts that back Sections 3–5 of the paper and the numerical values
reported in Section 4.

## Status

**Pre-registered:** 2026-04-19 (Zenodo; see `PRE_REGISTRATION.md`)
**Production run:** complete. Artifacts in `results/`.
**Paper version:** v2 (post-adversarial audit, 2026-04-20)
**Pre-registration commit:** `5ab2b88` (2026-04-19 JST).
**Production / audit-remediation commit:** `3f09519` (2026-04-20 JST).

## What this code tests

A single observation from the paper: the outer-Soddy curvature of the
Koide lepton triple, squared, equals the Standard-Model strange-quark
MS-bar mass at the lepton-sum scale within current lattice precision.

$$\mathcal{F}^2 = 95.1134 \text{ MeV} \quad \text{vs} \quad m_s^{\overline{\text{MS}}}(\mu_\star) = 95.075 \pm 0.692 \text{ MeV}$$

at $\mu_\star = m_e + m_\mu + m_\tau \approx 1883.10$ MeV. Observed
residual **+0.038 MeV / +0.055σ** against the lattice uncertainty.


## Six tests

| Test | Paper section | Question |
|------|---------------|----------|
| A    | §5.4          | Would a random Koide-compatible universe produce this residual? |
| B    | §5.6          | Is this residual distinguishable from today's measurement noise? |
| C    | §5.6          | Is $\mathcal{F}^2$ one of many similar functions that could have hit? |
| D    | §5.6          | Has the residual's historical trajectory behaved as a real relation's would? |
| filter sensitivity | §5.5 | How does the lower cutoff $\mu_\star > \mu_{\text{min}}$ affect the hit fraction? |
| alt scales | §5.3    | How does varying the renormalization-scale prescription affect the comparison? |

See `PRE_REGISTRATION.md` for exact pre-specifications of A–D and their
pre-committed seeds 2026–2029 (Test A) and 2126–2129 (Test B).

## Dependencies

- Python 3.9+
- numpy
- scipy
- [rundec](https://pypi.org/project/rundec/) (Python wrapper for CRunDec;
  Herren & Steinhauser 2018)

Pinned versions in `requirements.txt`:

```
rundec==0.7
numpy==2.4.4
scipy==1.17.1
```

The engine uses four-loop MS-bar mass running with threshold matching at
$m_b$ through `rundec`; see `engine/running.py`.


## Running it

```bash
# From repository root
pip install -r requirements.txt

python -m tests.validate_reproducibility        # smoke test, seconds
python -m tests.test_a_random_spectrum          # Null A, 4 priors, ~3 min
python -m tests.test_b_bootstrap                # Null B, 4 variants, ~1 min
python -m tests.test_c_algorithm_space          # 24-function enumeration, ~30 sec
python -m tests.test_d_temporal_convergence     # temporal trajectory, ~seconds
python -m tests.test_filter_sensitivity         # §5.5 cutoff scan, ~2 min
python -m tests.test_alt_scales                 # §5.3 scale prescriptions, ~5 sec
```

Each test writes a JSON artifact to `results/`. Core four tests (A/B/C/D):
under 5 minutes. Full robustness suite (core + filter sensitivity +
alternative scales): under 10 minutes.

## Repository layout

```
fermion_mass_mc/
├── README.md                     ← you are here
├── LICENSE                       ← MIT
├── requirements.txt              ← pinned dependencies
├── PRE_REGISTRATION.md           ← pre-registered MC specification
├── DESIGN.md                     ← architecture / design narrative
├── LITERATURE.md                 ← citation stack for prior shapes
├── engine/
│   ├── algorithm.py              ← Soddy pipeline (outer/inner curvature, F²)
│   ├── algorithm_variants.py     ← 24-function enumeration for Test C
│   ├── koide_manifold.py         ← analytic direct sampler on the Koide manifold
│   ├── priors.py                 ← four prior shapes (A1–A4)
│   ├── running.py                ← MS-bar running via rundec 4-loop
│   └── statistics.py             ← Clopper–Pearson, tail probs, hit counters
├── tests/
│   ├── test_a_random_spectrum.py       ← Null A (§5.4)
│   ├── test_b_bootstrap.py             ← Null B measurement-noise bootstrap (§5.6)
│   ├── test_c_algorithm_space.py       ← 24-variant enumeration (§5.6)
│   ├── test_d_temporal_convergence.py  ← temporal trajectory (§5.6)
│   ├── test_filter_sensitivity.py      ← §5.5 cutoff sensitivity
│   ├── test_alt_scales.py              ← §5.3 alternative scale prescriptions
│   ├── test_1_soddy_pipeline.py        ← earlier-generation pipeline validator
│   ├── test_2_joint_stack.py           ← earlier-generation pipeline validator
│   └── validate_reproducibility.py     ← smoke test
├── results/
│   ├── test_a_random_spectrum.json         ← Table 1 hit fractions + CPs
│   ├── test_b_bootstrap.json               ← §5.6 bootstrap variants
│   ├── test_c_algorithm_space.json         ← §5.6 Table 2 (24 variants)
│   ├── test_d_temporal_convergence.json    ← five-epoch FLAG trajectory
│   ├── test_filter_sensitivity.json        ← §5.5 four cutoffs
│   ├── test_alt_scales.json                ← §5.3 eight prescriptions
│   ├── historical_masses.json              ← FLAG + PDG trajectory (Test D input)
│   └── ms_history.csv                      ← CSV form of historical masses
├── figures/
│   └── fig_strange_convergence.png  ← Figure 2 of the paper
└── REPO/
    ├── README.md                            ← pointer notes
    ├── koide_papers_master_list.txt         ← ~60 Koide literature entries (cited as `KoideSoddyCode`)
    ├── empirical_mass_formulas_tradition.txt
    └── draft_paragraph_for_main_tex.tex
```


## Reproducibility

Every Monte Carlo result in the paper's §5.3, §5.4, §5.5, and §5.6 is
backed by a committed JSON artifact in `results/` produced by the
corresponding test driver. Re-running each test with the pre-committed
seeds reproduces the artifact byte-for-byte (MD5-verified 2026-04-20).
The §5.1 error-budget decomposition is an analytic propagation from
the stated input uncertainties and is verifiable by inspection, not an
MC output.

Four-loop CRunDec running is the one environment-sensitive piece of the
pipeline. The pinned `rundec==0.7` in `requirements.txt` matches the
build used to produce the committed artifacts.

## What the Monte Carlo does not provide

All probabilities reported are **empirical tail probabilities under the
specified null**. They are NOT $p$-values, are NOT severity assessments
in the Mayo–Spanos sense, and should not be composed across tests. The
Koide relation itself is treated as an empirical input rather than a
hypothesis under test; the MC measures the conditional frequency of
the strange-sector match *given* Koide. See paper §5.7 and Cousins
2018 (arXiv:1807.05996) for the framing.

## Terminology

- $\mathcal{F} = k_4^{-} = e_1 - 2\sqrt{e_2}$ — the exact outer Soddy
  curvature of the Koide lepton triple (dimensions MeV$^{1/2}$)
- $\mathcal{F}^2$ — its square, compared against $m_s^{\overline{\text{MS}}}(\mu_\star)$
- $e_k$ — elementary symmetric polynomial of degree $k$ in
  $\sqrt{m_\ell}$; $p_k$ — power sum. Koide exact is $p_2 = \tfrac{2}{3}e_1^2$.
- $\alpha_K = \sqrt{3/2} - 1$; $\alpha_K^2 = 5/2 - \sqrt{6}$
- $\mu_\star = m_e + m_\mu + m_\tau \approx 1883.10$ MeV
- Surrogate: $e_1 - \sqrt{p_2}$, equal to $k_4^{-}$ only under exact Koide;
  used as an algebraic simplification in some paper equations but the
  committed numerical comparisons use $k_4^{-}$ throughout.

## License

MIT (see `LICENSE`). Cited works retain their own copyrights.

## Citation

If you use this code please cite the paper:

```bibtex
@article{Brilliant2026strange,
  author  = {Brilliant, A. M.},
  title   = {A Descartes--Soddy completion of the Koide lepton triple},
  journal = {Preprints.org},
  year    = {2026},
  doi     = {TBD-Preprints.org pending}
}

@misc{Brilliant2026methods,
  author = {Brilliant, A. M.},
  title  = {Pre-registered Monte Carlo protocols for the outer Soddy--Koide numerical observation},
  year   = {2026},
  doi    = {10.5281/zenodo.PENDING},
  url    = {https://doi.org/10.5281/zenodo.PENDING},
  note   = {Companion methods note; Zenodo DOI pending upload}
}
```

and link back to this repository.
