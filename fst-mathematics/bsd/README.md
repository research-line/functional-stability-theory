# BSD Positivity Normal Form

Public reproducibility package for the BSD domain supplement in the Functional Stability Theory program.

Zenodo concept DOI: <https://doi.org/10.5281/zenodo.19087443>

## Status

This folder contains the public v1.5 paper files released on Zenodo (Record <https://doi.org/10.5281/zenodo.21610295> and <https://doi.org/10.5281/zenodo.21610292>). Version 1.5 incorporates the 3-page frontmatter architecture, Chicago/APA title casing, author-pair en-dash typography, comprehensive RevTeX 4-2 table hardening (Tables 1–8), and the 4-tier TikZ vector architecture schema (`fig:bsd_architecture`) visualizing the positivity normal-form structure. In strict compliance with research governance (§ 79 and § 80), all mathematical statements, proofs, theorems, and epistemic boundaries remain unchanged: rank <= 1 is verified through the Gross-Zagier/Kolyvagin regime, while the rank >= 2 Higher Gross-Zagier bridge remains open.

## Files

| File | Purpose | Pages |
|------|---------|-------|
| `BSD_Positivity_EN.tex` / `BSD_Positivity_EN.pdf` | English paper source and PDF | 29 |
| `BSD_Positivity_DE.tex` / `BSD_Positivity_DE.pdf` | German paper source and PDF | 31 |
| `BSD_Positivity_kombi.pdf` | Combined bilingual EN+DE PDF with interactive bookmarks | 60 |
| `../../scripts/bsd/compute_bsd_verification.py` | BSD formula sanity checks for selected LMFDB curves | - |
| `../../scripts/bsd/compute_height_saturation.py` | Rank-1 identity and quadratic-twist heuristic plot | - |
| `../../scripts/bsd/compute_rank2_lmfdb.py` | Rank-2 regulator positivity sample and plot | - |
| `../../scripts/bsd/compute_height_saturation.png` | Generated plot from the height-saturation script | - |
| `../../scripts/bsd/compute_rank2_lmfdb.png` | Generated plot from the rank-2 regulator script | - |

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
