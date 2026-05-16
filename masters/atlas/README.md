# Dirichlet Character Atlas

**The Dirichlet Character Atlas — A Weil-Kernel Numerical Atlas for Dirichlet L-Function Sector Gaps (Galerkin Diagnostics and the Boundaries of Leading-Order Theory)**

The Atlas is the *micro-cartography* companion of the Zeta Zoo: where the Zeta Zoo provides macro-classification (SGE taxonomy), the Atlas resolves a single Dirichlet family at three Galerkin resolutions and exposes the structural limits of leading-order numerics.

## DOI

- **Concept-DOI (always latest):** [10.5281/zenodo.19960809](https://doi.org/10.5281/zenodo.19960809)
- **v2 (current):** [10.5281/zenodo.20241612](https://doi.org/10.5281/zenodo.20241612)
- **v1:** [10.5281/zenodo.19960810](https://doi.org/10.5281/zenodo.19960810)

## Status

| | |
|---|---|
| Status | DRAFT preprint (negative-result atlas + structural diagnosis) |
| Review | 7-phase review chain completed 2026-04-29 (score 9.5/10) |
| Version | v2 |
| Pages (EN / DE / Kombi) | 13 / 14 / 27 |
| Resolution | N ∈ {200, 400, 600} |
| Real Dirichlet characters covered | D ∈ {5, 8, 12, 13, 17, 21, 24, 29, 33, 60} |

## Main results

1. **Falsification** of the leading-order approximation φ⁺ ≈ φ⁻ (empirical, from Session-6 data).
2. **Tautology diagnosis** of the full Galerkin formula (Session-7 T5: matrix-decomposition tautology).
3. **N-truncation artifact** identification for χ₂₁ (oscillation vanishes in higher resolutions).
4. **Roadmap** for the Dirichlet-twisted CCM Route D (Zookeeper follow-up).

## Files

| File | Purpose |
|---|---|
| `paper/DIRICHLET_CHARACTER_ATLAS_v2_en.pdf` | English version (13 p.) |
| `paper/DIRICHLET_CHARACTER_ATLAS_v2_de.pdf` | German version (14 p.) |
| `paper/DIRICHLET_CHARACTER_ATLAS_v2_kombi.pdf` | Combined EN+DE (27 p.) |
| `scripts/` | Galerkin computation pipeline (Python; basis, kappa grid, asymptotic scans, χ-specific tests) |
| `results/` | Numerical results (JSON) and analysis reports (Markdown) |

## v2 change note

Version v2 is a source-correction release. The LMFDB passage was narrowed from a broad zero-certificate claim to computed low-lying zeros as finite consistency data; EN/DE/Kombi PDFs were rebuilt and PDF metadata was normalized before upload.

## Reproducibility

Main computation entry points in `scripts/`:

- `analytic_gap_formula_test.py` / `_v2.py` — analytic gap formula tests
- `basis_kappa_full_grid.py` — basis and κ scan
- `asymptotic_scan_server.py` — N-asymptotic scan (server-side)
- `all10_high_N_server.py` — full ten-character high-N study
- `arch_term_analysis.py` — archimedean term diagnostics

Numerical outputs (JSON + per-experiment MD reports) live in `results/`.

## Position in the FST programme

The Atlas is the **fourth CoreCore Master** of FST — together with the Zookeeper (RH proof), Zeta Zoo (classification), and Spectrum Duality / RFEP (physics). Atlas + Selberg form the *method-validation pair*: the Atlas is the **negative** test (where leading-order Galerkin diagnostics fall short for Dirichlet characters), Selberg is the **positive** test (where v2.0 reproduces a classical operator-based result).

## License

CC-BY-4.0
