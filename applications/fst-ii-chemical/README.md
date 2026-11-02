# FST-II — Chemical Stability and Evolutionary Game Dynamics

Companion repository and reproducibility artifacts for:
**Functional Stability Theory — Part II: Chemical Reaction Networks, Evolutionary Game Dynamics, and Non-Normal Resolvent Bounds**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21844410.svg)](https://doi.org/10.5281/zenodo.21844410)
[![Concept DOI](https://img.shields.io/badge/Concept%20DOI-10.5281%2Fzenodo.20130563-blue.svg)](https://doi.org/10.5281/zenodo.20130563)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-blue.svg)](https://creativecommons.org/licenses/by/4.0/)

## Overview

This application module provides mathematical proofs, numerical simulations, and bilingual LaTeX manuscripts for the application of Functional Stability Theory (FST) to chemical reaction networks, replicator equations, and non-normal reaction dynamics.

### Core Theoretical Contributions

1. **Permanence & Boundary Transversality (Theorem 1):** Strict proof of interior simplex permanence and inward transversal vector flows for open autocatalytic networks.
2. **Non-Normal Resolvent & Pseudospectral Stability (Theorem 2):** Exact calculation of the non-normal transient amplification (peak growth $\approx 1.2056$ at $t \approx 1.843$), numerical abscissa $\mu_2 > 0$ under negative spectral abscissa $\alpha_T < 0$, exact $H_\infty$ resolvent norm bound ($M_{\mathrm{res}} \approx 6.0353$), and stability radius contrast ($> 1.88\times$ over naive eigenvalues).
3. **Replicator-Nash & Potential-Lyapunov Demarcation (Theorem 3):** Strict demarcation establishing that interior limit points of replicator dynamics are not generically Nash equilibria, and identifying tangent space negative definiteness for Kullback--Leibler divergence monotonicity.

## Repository Contents

- `papers/`: Bilingual LaTeX sources and compiled PDFs:
  - `FST_II_Chemical_Stability_v3_en.pdf` (English manuscript)
  - `FST_II_Chemical_Stability_v3_de.pdf` (German manuscript)
  - `FST_II_Chemical_Stability_v3_kombi.pdf` (Combined bilingual publication)
- `scripts/`: Reproducibility, gate verification, and ledger scripts:
  - `verify_fst2_permanence_resolvent_and_nash_gate.py`: Verifies Theorems 1–3 and resolvent bounds.
  - `verify_fst2_transient_peak_2026_10_18.py`: High-precision verified matrix exponential and transient peak derivation.
  - `autocatalysis_3species.py`: Numerical ODE simulation of 3-species autocatalytic replicator system.
  - `homochirality_mechanism_taxonomy_gate.py`: Homochirality mechanism gate audit.
  - `homochirality_open_drive_gate.py`: Open drive enantiomeric excess drift controls.
  - `dissipation_topology_waterline_ledger.py`: Waterline ledger for dissipation topology.
- `results/`: Numerical benchmark outputs and JSON data.

## Verification

Run automated contract tests across the FST-II domain:

```bash
pytest tests/test_compute_fst2.py -v
```

Author: Lukas Geiger | Functional Stability Theory Program
