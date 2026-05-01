# The Zeta Zoo

**Branch-Local Hilbert-Polya via the Trinity of UCU, SGE, and Weil-Form Transversality**

Lukas Geiger (2026)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19673227.svg)](https://doi.org/10.5281/zenodo.19673227)

## Summary

The mathematical classification paper of the FST programme. Introduces:

- **UCU** (Universal Convexity Uniqueness): strict convexity implies unique saturation
- **SGE** (Semigroup-Group Equivalence): HP-BL as binary algebraic invariant
- **Weil-form transversality**: gap positivity across both sides of the SGE dichotomy

## Zoo Population (v1.1)

| Family | HP-BL Class | Mechanism |
|--------|-------------|-----------|
| Riemann zeta | NO (classical) | NE-A + NE-B obstructions |
| Selberg zeta | YES | Casimir / Laplace-Beltrami |
| Prime-Hub | OPEN | Three obstructions documented |
| CRM flow-zeta | YES | Ruelle-type, Lie group flow |
| Dedekind Q(sqrt(-5)) | NO | Prime-ideal semigroup |
| Ihara-Petersen | YES | Bass-Hashimoto matrix |

## Companion Papers

- [The Spectral Zookeeper](../zookeeper/) -- RH proof via microcluster closure (DOI: 10.5281/zenodo.19673127)
- [FST Spectrum Duality / RFEP](../spectrum-duality/) -- Physical instantiation (DOI: 10.5281/zenodo.19162705)

## Paper Files

- [`paper/NE_B_BOUNDARY_v1_en.tex`](paper/NE_B_BOUNDARY_v1_en.tex)
- [`paper/NE_B_BOUNDARY_v1_en.pdf`](paper/NE_B_BOUNDARY_v1_en.pdf)
- [`paper/NE_B_BOUNDARY_v1_de.tex`](paper/NE_B_BOUNDARY_v1_de.tex)
- [`paper/NE_B_BOUNDARY_v1_de.pdf`](paper/NE_B_BOUNDARY_v1_de.pdf)
- [`paper/NE_B_BOUNDARY_v1_kombi.pdf`](paper/NE_B_BOUNDARY_v1_kombi.pdf)

## Reproducibility

Validation scripts are kept in [`../../scripts/zeta-zoo/`](../../scripts/zeta-zoo/):

- `dedekind_ne_b_test.py` -- Dedekind Q(sqrt(-5)) NE-B analog probe
- `ihara_petersen_sge_test.py` -- Ihara/Petersen SGE YES-side test
- `sge_control_experiment.py` -- SGE YES/NO discriminating control experiment

Compact public outputs are stored in [`results/`](results/) as JSON plus Markdown summaries. Runtime logs and internal proof/planning notes are intentionally not part of the public repository.
