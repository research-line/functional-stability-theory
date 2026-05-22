# Yang-Mills Mass Gap

Public reproducibility package for the Yang-Mills domain supplement in the
Functional Stability Theory program.

Zenodo concept DOI: <https://doi.org/10.5281/zenodo.19087433>

## Status

This folder contains the public post-v2.5 guardrail candidate synchronized
from the local strict-review, design-check, and citation-checked working
state. The latest released Zenodo record remains v2.5:
<https://doi.org/10.5281/zenodo.20315811>. The GitHub tree is therefore
slightly ahead of Zenodo in claim-boundary wording and supporting public
artifacts, while the continuum mass-gap step remains conditional on the
analytical renormalization-group contraction input.

## Files

| File | Purpose |
|------|---------|
| `FST-YM_YangMills_MassGap_en.tex` / `FST-YM_YangMills_MassGap_en.pdf` | English paper source and PDF |
| `FST-YM_YangMills_MassGap_de.tex` / `FST-YM_YangMills_MassGap_de.pdf` | German paper source and PDF |
| `FST-YM_YangMills_MassGap_kombi.pdf` | Combined bilingual PDF |
| `../../scripts/yang-mills/compute_dobrushin_su2.py` | SU(2) lattice Dobrushin influence scan |
| `../../scripts/yang-mills/compute_birkhoff_rg.py` | Birkhoff contraction scan for hierarchical RG steps |
| `../../scripts/yang-mills/compute_dobrushin_su2.png` | Generated Dobrushin result plot |
| `../../scripts/yang-mills/compute_birkhoff_rg.png` | Generated Birkhoff/RG result plot |

## Reproduce

From the repository root:

```bash
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_dobrushin_su2.py
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_birkhoff_rg.py
```

The scripts write their PNG outputs next to the scripts in
`scripts/yang-mills/`.

## Publication Gate

Internal proof notes, `_results/`, design-check directories, extracted text
artefacts, review chains, planning files, Zenodo credentials, and private
comparison notes are intentionally not part of this public package. They remain
local-only until the project reaches the required completion gate.
