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

## Published Papers

| Paper | Domain | DOI | Version |
|-------|--------|-----|---------|
| Unified Framework (RFEP) | Meta-Framework | [10.5281/zenodo.19036190](https://doi.org/10.5281/zenodo.19036190) | v1.5 |
| Yang-Mills Mass Gap | Mathematical Physics | [10.5281/zenodo.19087433](https://doi.org/10.5281/zenodo.19087433) | v1.0 |
| Navier-Stokes Regularity | Mathematical Physics | [10.5281/zenodo.19087449](https://doi.org/10.5281/zenodo.19087449) | v1.1 |
| NS Log-Distance Integrability | Navier-Stokes | [10.5281/zenodo.19056807](https://doi.org/10.5281/zenodo.19056807) | v1.2 |
| Turbulence (K41 Spectrum) | Turbulence | [10.5281/zenodo.19056813](https://doi.org/10.5281/zenodo.19056813) | v1.2 |
| Dark Energy (Cosmological Constant) | Cosmology | [10.5281/zenodo.19036235](https://doi.org/10.5281/zenodo.19036235) | v1.5 |
| Hodge Conjecture (Arithmetic Positivity) | Algebraic Geometry | [10.5281/zenodo.19087439](https://doi.org/10.5281/zenodo.19087439) | v1.0 |
| BSD Conjecture (Normal Form) | Number Theory | [10.5281/zenodo.19087443](https://doi.org/10.5281/zenodo.19087443) | v1.0 |
| P vs NP (Entropic No-Go) | Complexity Theory | [10.5281/zenodo.19056809](https://doi.org/10.5281/zenodo.19056809) | v1.1 |
| FST-RH Trilogy | Riemann Hypothesis | [10.5281/zenodo.19035845](https://doi.org/10.5281/zenodo.19035845) | Published |
| CRM I-IV | Cosmology | [10.5281/zenodo.18728936](https://doi.org/10.5281/zenodo.18728936) | Published |

## Numerical Validation Scripts

The `scripts/` directory contains computational validation scripts:

| Script | Paper | Description |
|--------|-------|-------------|
| `compute_F_spectrum.py` | Turbulence | Verifies K41 as unique minimiser of F[E]; strict convexity test |
| `compute_ds3_lorenz.py` | Navier-Stokes | DS3 stress test on Lorenz attractor; TV saturation |
| `compute_w_vs_desi.py` | Dark Energy | w_eff(z) comparison with DESI constraints |
| `compute_w_mapping.py` | Dark Energy | Correct w_eff -> w_DE mapping + DESI grid scan |
| `compute_height_saturation.py` | BSD | Height saturation test for quadratic twists |
| `compute_ghr_spectrum.py` | Hodge | GHR spectrum numerical verification |

## Related Repositories

- [rfep-framework](https://github.com/research-line/rfep-framework) -- Renormalized Free-Energy Principle (the mathematical bridge)
- [rh-even-dominance](https://github.com/research-line/rh-even-dominance) -- Riemann Hypothesis: Even-dominance proof (foundation)
- [crm-cosmology](https://github.com/research-line/crm-cosmology) -- Cosmic Recursion Model (foundation)

## Author

Lukas Geiger -- ORCID: [0009-0005-7296-1534](https://orcid.org/0009-0005-7296-1534)

## License

This work is licensed under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).
