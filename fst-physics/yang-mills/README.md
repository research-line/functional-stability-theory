# Yang-Mills Mass Gap

Public reproducibility package for the Yang-Mills domain supplement in the
Functional Stability Theory program.

Zenodo concept DOI: <https://doi.org/10.5281/zenodo.19087433>

## Status

This folder contains the public v2.2 candidate paper files synchronized from
the local Post-Zookeeper/RFEP working state. The latest Zenodo release is still
v2.1. The continuum mass-gap step remains conditional on the analytical
renormalization-group contraction input.

## Files

| File | Purpose |
|------|---------|
| `FST-YM_YangMills_Skeleton_v1_en.tex` / `FST-YM_YangMills_Skeleton_v1_en.pdf` | English paper source and PDF |
| `FST-YM_YangMills_Skeleton_v1_de.tex` / `FST-YM_YangMills_Skeleton_v1_de.pdf` | German paper source and PDF |
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

Internal proof notes, review chains, planning files, Zenodo credentials, and
private comparison notes are intentionally not part of this public package.
They remain local-only until the project reaches the required completion gate.
