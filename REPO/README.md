# Koide Literature Repository

A curated, citation-verified bibliography of the Koide formula literature
and the broader tradition of empirical mass formulas in particle physics.

## Contents

- [`koide_papers_master_list.txt`](koide_papers_master_list.txt) — the
  master bibliography of Koide-specific papers, organized by topic
  (original papers, geometric interpretations, quark/neutrino extensions,
  dynamical explanations, RG-running, reviews, critical analyses).
- [`empirical_mass_formulas_tradition.txt`](empirical_mass_formulas_tradition.txt) —
  context on the broader tradition of empirical mass formulas in particle
  physics (Balmer, Gell-Mann-Okubo, Barut, Harari-Haut-Weyers, Regge
  trajectories, Fritzsch texture ansatze, anarchy hypothesis).
- This `README.md`.

## Purpose

The Koide formula (Y. Koide, 1981-82) is a strikingly precise empirical
relation among the charged lepton masses:

> Q = (m_e + m_mu + m_tau) / (sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau))^2 = 2/3

that has held to within experimental uncertainty across 40+ years and a
factor of 35 improvement in tau-mass measurement precision. The formula
has a substantial literature (~60-100 papers) spanning geometric
interpretations, extensions to quarks and neutrinos, dynamical models
(most notably Sumino's family gauge symmetry), RG-running analyses,
and reviews.

This repository exists to:

1. Provide a maintained, citation-verified list of the Koide literature.
2. Contextualize the Koide formula within the broader tradition of
   empirical mass formulas in particle physics, several of which
   (Balmer, Gell-Mann-Okubo, Regge) received dynamical explanations
   decades after the empirical discovery.
3. Support future work applying geometric, statistical, or landscape
   methods to fermion mass structure.

## Verification status

Entries are marked:

- `[V]` Verified against primary source (arXiv/DOI/journal) within past week.
- `[S]` From secondary source (Wikipedia, review paper). Verify before citing.
- `[R]` From a review paper's citation list.
- `[?]` Uncertain; verify before using.

## How to contribute

When adding an entry:

1. Provide full bibliographic citation (authors, title, journal, vol,
   page, year) AND arXiv number AND (ideally) DOI.
2. Mark verification status.
3. Write a one-paragraph `Comment:` explaining what the paper contributes.
4. Place under the appropriate section.

When two sources give inconsistent citation details for the same paper,
prefer in order: (1) DOI resolution, (2) arXiv abstract page,
(3) INSPIRE-HEP record, (4) journal publisher's page. Wikipedia and
HandWiki are tertiary.

## License

Content in this repository is made available under the MIT License for
the bibliography compilation and curatorial choices. The individual
citations point to third-party works governed by their own copyrights.

## Related

This repository accompanies the paper *"The Soddy formula applied to the
Standard Model mass spectrum"* (A. Brilliant, 2026), which extends the
Koide lepton relation via the Descartes circle theorem and Soddy's
formula to predict light quark masses. The Monte Carlo null-hypothesis
test of that paper's claims is in the sibling directory
`mass_as_curvature_mc_v4/`.

## Maintainer

A. Brilliant (Anthropic Research affiliation as of 2026).

Questions, corrections, and pull requests welcome.
