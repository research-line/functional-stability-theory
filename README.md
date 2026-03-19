# Functional Stability Theory (FST)

[![DOI Framework](https://zenodo.org/badge/DOI/10.5281/zenodo.19036190.svg)](https://doi.org/10.5281/zenodo.19036190)
[![DOI Dark Energy](https://zenodo.org/badge/DOI/10.5281/zenodo.19036235.svg)](https://doi.org/10.5281/zenodo.19036235)
[![DOI Yang-Mills](https://zenodo.org/badge/DOI/10.5281/zenodo.19087433.svg)](https://doi.org/10.5281/zenodo.19087433)
[![DOI Hodge](https://zenodo.org/badge/DOI/10.5281/zenodo.19087439.svg)](https://doi.org/10.5281/zenodo.19087439)
[![DOI BSD](https://zenodo.org/badge/DOI/10.5281/zenodo.19087443.svg)](https://doi.org/10.5281/zenodo.19087443)

**Functional Stability Theory** is a unified mathematical framework that connects thermodynamic stability, game-theoretic equilibrium selection, and renormalized energy functionals to derive rigorous results across multiple domains of mathematical physics and complexity theory.

The central object is a **renormalized free-energy functional** whose critical-point structure, combined with selection principles from evolutionary game theory, yields domain-specific theorems when instantiated with the appropriate state spaces and energy densities. The framework is governed by the **Dissipative Selection Principle** (DSP): among all critical points of the renormalized functional, physical evolution selects the unique dissipation-minimising attractor -- a meta-theorem that unifies the domain-specific results below.

## Repository Structure

- **framework/** -- The Unified Renormalized Energy Framework (meta-level paper)
- **domain-proofs/** -- Domain-specific instantiations and proofs
- **applications/** -- Empirical applications (FST-I Particles, FST-II Chemistry, FST-III Biology)
- **scripts/** -- Numerical validation scripts
- **fst_references.bib** -- Shared bibliography

## Domain Proof Papers -- Current Status

| Paper | Version | Status | Open Problem | Next Step |
|-------|---------|--------|--------------|-----------|
| [**Turbulence**](https://doi.org/10.5281/zenodo.19056813) | v1.3 | Journal-ready | DFC1 empirical (only input) | Figures + falsification protocol |
| [**Dark Energy**](https://doi.org/10.5281/zenodo.19036235) | v1.6 | Framework Note | Hu-Sawicki parameters quantitatively open | MCMC fit against DESI+Planck+Cassini |
| [**Yang-Mills**](https://doi.org/10.5281/zenodo.19087433) | v2.1 | Conditional | Analytical proof of lambda < 0 | Doeblin minorisation programme |
| [**Navier-Stokes**](https://doi.org/10.5281/zenodo.19087449) | v2.1 | Conditional | Assumption G2 (projection regularity) | Independent geometric verification |
| [**NS Log-Distance**](https://doi.org/10.5281/zenodo.19056807) | v1.3 | Proof of Life verified | TLL for 3D NS analytically open | Kuramoto-Sivashinsky / reaction-diffusion |
| [**BSD**](https://doi.org/10.5281/zenodo.19087443) | v1.1 | Reformulation | Higher Gross-Zagier (rank >= 2) | Rank 2-4 numerics, conditional labelling |
| [**Hodge**](https://doi.org/10.5281/zenodo.19087439) | v1.1 | No-Go Theorem | = Deligne's question (1982) | Prismatic cohomology / AP4 |
| [**P vs NP**](https://doi.org/10.5281/zenodo.19056809) | v1.2 | Reformulation | Uniformity Bridge | Instance compression formalisation |
| [**Framework (RFEP)**](https://doi.org/10.5281/zenodo.19036190) | v1.6 | Meta-Theorem | Pattern A falsifiability clarified | -- |

## Proof Architecture

```
                    Framework (Pattern A)
                    +-- DS1-DS3 (Dissipative Selection)
                    +-- Second-Order Resolvent Dominance
                         |
          +--------------+--------------+
          |              |              |
     PROVEN        CONDITIONAL        OPEN
     (rigorous)   (reduced to       (= open
                  threshold axiom)   research)
          |              |              |
         TU:            YM:           BSD:
      DFC => NL'    Doeblin a > 0   Rank >= 2
      (journal)     (Kingman l < 0)  Gross-Zagier
                        |
         DE:           NS:          Hodge:
      Screening    G2 (projection)  Deligne's
      (validated)  G3 (Gronwall)    question
                        |
        NS-LDI:       PvNP:
      TLL + LDI     Uniformity
      (Lorenz OK)    Bridge
```

## Numerical Validation Scripts

The `scripts/` directory contains computational validation scripts:

| Script | Paper | Description |
|--------|-------|-------------|
| `compute_F_spectrum.py` | Turbulence | Verifies K41 as unique minimiser of F[E]; strict convexity test |
| `compute_ds3_lorenz.py` | Navier-Stokes | DS3 stress test on Lorenz attractor; TV saturation |
| `compute_tll_ldi_lorenz.py` | NS-LDI | **Proof of Life**: TLL+LDI on Lorenz attractor (5/5 tests passed) |
| `compute_w_vs_desi.py` | Dark Energy | w_eff(z) comparison with DESI constraints |
| `compute_w_mapping.py` | Dark Energy | Correct w_eff -> w_DE mapping + DESI grid scan |
| `compute_height_saturation.py` | BSD | Height saturation test for quadratic twists |
| `compute_ghr_spectrum.py` | Hodge | GHR spectrum numerical verification |

## Related Repositories

- [rh-even-dominance](https://github.com/research-line/rh-even-dominance) -- Riemann Hypothesis: Even-dominance proof (foundation)
- [crm-cosmology](https://github.com/research-line/crm-cosmology) -- Curvature Relaxation Model (foundation)
- [rfep-framework](https://github.com/research-line/rfep-framework) -- archived, now integrated here under `framework/`

## Author

Lukas Geiger -- ORCID: [0009-0005-7296-1534](https://orcid.org/0009-0005-7296-1534)

## License

This work is licensed under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).
