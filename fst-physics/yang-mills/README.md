# Yang-Mills Mass Gap

Public reproducibility package for the Yang-Mills domain supplement in the
Functional Stability Theory program.

Zenodo concept DOI: <https://doi.org/10.5281/zenodo.19087433>
Latest Zenodo v2.6 record DOI: <https://doi.org/10.5281/zenodo.20716608>

## Status

This folder contains the public v2.6 release package synchronized from the
local strict-review, design-check, and citation-checked working state. The
Zenodo v2.6 record is live under the record DOI above. The continuum mass-gap
step remains conditional on the analytical renormalization-group contraction
input.

## Files

| File | Purpose |
|------|---------|
| `FST-YM_YangMills_MassGap_en.tex` / `FST-YM_YangMills_MassGap_en.pdf` | English paper source and PDF |
| `FST-YM_YangMills_MassGap_de.tex` / `FST-YM_YangMills_MassGap_de.pdf` | German paper source and PDF |
| `FST-YM_YangMills_MassGap_kombi.pdf` | Combined bilingual PDF |
| `../../scripts/yang-mills/compute_dobrushin_su2.py` | SU(2) lattice Dobrushin influence scan |
| `../../scripts/yang-mills/compute_birkhoff_rg.py` | Birkhoff contraction scan for hierarchical RG steps |
| `../../scripts/yang-mills/compute_os_capacity_ledger.py` | OS-danger capacity ledger and negative-control diagnostic |
| `../../scripts/yang-mills/compute_rp_os_transfer_ledger.py` | RP/OS transfer matrix positivity ledger |
| `../../scripts/yang-mills/compute_rp_os_rfep_transfer_ledger.py` | RFEP transfer matrix diagnostic ledger |
| `../../scripts/yang-mills/compute_u1_2d_strong_coupling_positive_control.py` | 2D U(1) strong-coupling positive control ledger |
| `../../scripts/yang-mills/compute_ym_waisen_transfer_ledger.py` | Yang-Mills waisen transfer verification ledger |
| `../../scripts/yang-mills/u1_strong_coupling_positive_control.py` | Direct U(1) character expansion positive control script |
| `../../scripts/yang-mills/compute_dobrushin_su2.png` | Generated Dobrushin result plot |
| `../../scripts/yang-mills/compute_birkhoff_rg.png` | Generated Birkhoff/RG result plot |

## Reproduce

From the repository root:

```bash
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_dobrushin_su2.py
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_birkhoff_rg.py
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_os_capacity_ledger.py
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_rp_os_transfer_ledger.py
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_rp_os_rfep_transfer_ledger.py
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_u1_2d_strong_coupling_positive_control.py
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_ym_waisen_transfer_ledger.py
PYTHONIOENCODING=utf-8 python scripts/yang-mills/u1_strong_coupling_positive_control.py
```

The plotting scripts write their PNG outputs next to the scripts in
`scripts/yang-mills/`. The ledger script writes local `_data/` and `_results/`
folders under `scripts/yang-mills/`; those generated working outputs are not
versioned here.

## Publication Gate

Internal proof notes, `_results/`, design-check directories, extracted text
artefacts, review chains, planning files, Zenodo credentials, and private
comparison notes are intentionally not part of this public package. They remain
local-only until the project reaches the required completion gate.
