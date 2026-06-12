# BSD Positivity Normal Form

Public reproducibility package for the BSD domain supplement in the Functional Stability Theory program.

Zenodo concept DOI: <https://doi.org/10.5281/zenodo.19087443>

## Status

This folder contains the public v1.4 paper files released on Zenodo as record <https://doi.org/10.5281/zenodo.20671962>. Version 1.4 is a maintenance release: it publishes the June 2026 citation, source-metadata, and English-style corrections after v1.3, without upgrading the theorem status or adding a new BSD proof claim. The paper remains a structural reformulation: rank <= 1 is verified through the Gross-Zagier/Kolyvagin regime, while the rank >= 2 Higher Gross-Zagier bridge remains open.

## Files

| File | Purpose |
|------|---------|
| `BSD_Positivity_EN.tex` / `BSD_Positivity_EN.pdf` | English paper source and PDF |
| `BSD_Positivity_DE.tex` / `BSD_Positivity_DE.pdf` | German paper source and PDF |
| `BSD_Positivity_kombi.pdf` | Combined EN+DE PDF released with v1.4 |
| `../../scripts/bsd/compute_bsd_verification.py` | BSD formula sanity checks for selected LMFDB curves |
| `../../scripts/bsd/compute_height_saturation.py` | Rank-1 identity and quadratic-twist heuristic plot |
| `../../scripts/bsd/compute_rank2_lmfdb.py` | Rank-2 regulator positivity sample and plot |
| `../../scripts/bsd/compute_height_saturation.png` | Generated plot from the height-saturation script |
| `../../scripts/bsd/compute_rank2_lmfdb.png` | Generated plot from the rank-2 regulator script |

## Reproduce

From the repository root:

```bash
PYTHONIOENCODING=utf-8 python scripts/bsd/compute_bsd_verification.py
PYTHONIOENCODING=utf-8 python scripts/bsd/compute_height_saturation.py
PYTHONIOENCODING=utf-8 python scripts/bsd/compute_rank2_lmfdb.py
```

The plot scripts write their PNG outputs next to the scripts in `scripts/bsd/`.

## Publication Gate

Internal proof notes, review chains, planning files, and Zenodo credentials are intentionally not part of this public package. They remain local-only until the project reaches the required completion gate.
