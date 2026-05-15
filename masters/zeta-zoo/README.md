# The Zeta Zoo

**The mathematical side of Functional Stability Theory**

Lukas Geiger (2026)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19673226.svg)](https://doi.org/10.5281/zenodo.19673226)

## Summary

The Zeta Zoo is the mathematical classification paper of the FST programme. It develops the branch-local Hilbert-Polya taxonomy around three organizing principles:

- **UCU**: Universal Convexity Uniqueness.
- **SGE**: Semigroup-Group Equivalence as a HP-BL classification axis.
- **Weil-form transversality**: gap positivity across both sides of the SGE dichotomy.

Version 2.0 is a source-audit release. It updates internal Zenodo references to stable Concept-DOIs, removes a plaintext email from the paper, and publishes English, German, and combined PDFs as the current public file set.

## Zoo Population

| Family | HP-BL class | Mechanism |
|---|---|---|
| Riemann zeta | NO (classical) | NE-A + NE-B obstructions |
| Selberg zeta | YES | Casimir / Laplace-Beltrami |
| Prime-Hub | OPEN | OP5 construction under documented obstructions |
| CRM flow-zeta | YES | Ruelle-type / Lie group flow |
| Dedekind Q(sqrt(-5)) | NO | Prime-ideal semigroup |
| Ihara-Petersen | YES | Bass-Hashimoto matrix |

## Companion Papers

- [The Spectral Zookeeper](../zookeeper/) -- CCM microcluster route, Concept DOI: `10.5281/zenodo.19673126`.
- [FST Spectrum Duality / RFEP](../spectrum-duality/) -- physical instantiation, Concept DOI: `10.5281/zenodo.19036190`.
- [Dirichlet Character Atlas](../atlas/) -- negative method validation, Concept DOI: `10.5281/zenodo.19960809`.
- [Selberg](../selberg/) -- SGE-YES validation, Concept DOI: `10.5281/zenodo.19962588`.

## Paper Files

- [`paper/NE_B_BOUNDARY_v2_en.tex`](paper/NE_B_BOUNDARY_v2_en.tex)
- [`paper/NE_B_BOUNDARY_v2_en.pdf`](paper/NE_B_BOUNDARY_v2_en.pdf)
- [`paper/NE_B_BOUNDARY_v2_de.tex`](paper/NE_B_BOUNDARY_v2_de.tex)
- [`paper/NE_B_BOUNDARY_v2_de.pdf`](paper/NE_B_BOUNDARY_v2_de.pdf)
- [`paper/NE_B_BOUNDARY_v2_kombi.pdf`](paper/NE_B_BOUNDARY_v2_kombi.pdf)

## Reproducibility

Validation scripts are kept in [`../../scripts/zeta-zoo/`](../../scripts/zeta-zoo/):

- `dedekind_ne_b_test.py` -- Dedekind Q(sqrt(-5)) NE-B analogue probe.
- `ihara_petersen_sge_test.py` -- Ihara/Petersen SGE YES-side test.
- `sge_control_experiment.py` -- SGE YES/NO discriminating control experiment.
- `mobius_graded_primehub_toy.py` -- finite-prime toy model for the Mobius-graded Prime-Hub determinant.
- `prime_weil_bridge_defect.py` -- diagnostic and operator-Fourier Prime-Weil bridge defect probes.

Compact public outputs are stored in [`results/`](results/) as JSON plus Markdown summaries. Runtime logs, internal derivation notes, handoffs, and planning files are intentionally not part of the public repository.
