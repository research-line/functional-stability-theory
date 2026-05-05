# The Spectral Zookeeper

**The Riemann Hypothesis via Spectral Microcluster Closure in the CCM Framework**

Lukas Geiger (2026)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19673126.svg)](https://doi.org/10.5281/zenodo.19673126)

## Structure

```text
zookeeper/
|-- paper/
|   |-- RH_Zookeeper_v1_en.tex      Main paper source
|   |-- RH_Zookeeper_v1_en.pdf      Main paper PDF
|   |-- RH_Zookeeper_v1_de.tex      German paper source
|   |-- RH_Zookeeper_v1_de.pdf      German paper PDF
|   |-- RH_Zookeeper_v1_kombi.pdf   Combined EN+DE PDF
|   `-- fst-rh-references.bib       Bibliography
|-- scripts/                        Numerical verification
|   |-- c2bt_spectral_mass.py       Three-lemma mass-based gap computation
|   |-- c2_quasimode_scaling.py     Quasimode scaling verification
|   |-- c2_cancellation_*.py        Secular cancellation scans
|   `-- b10_analytical_derivation.py Resolvent identity derivation
`-- results/                        Compact numerical logs and JSON/MD outputs
```

## Key Results

The paper studies the Riemann Hypothesis conditionally on external CCM inputs via a three-lemma microcluster closure:

1. **Scalar Secular Cancellation** - |mu| <= 3e-7 (12-digit verified)
2. **Projected Poisson Quasimode** - ||h|| <= 3e-6 (12-digit verified)
3. **Coercive Complement** - g* >= 5 (mass-based effective gap)

Combined: ||(I - P_V)k|| <= 6e-7, implying MS2, hence RH via CCM + Hurwitz.

## Public / Private Boundary

This directory contains the public preprint, reproducibility scripts, and compact numerical results for the published Zookeeper package. Internal proof logs, review chains, planning files, credentials, and `BEWEISNOTIZ.md` remain excluded by repository ignore rules until a later journal-submission or completion gate.
