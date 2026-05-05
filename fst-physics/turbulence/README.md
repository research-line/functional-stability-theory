# Turbulence Cascade

Public reproducibility package for the turbulence domain supplement in the
Functional Stability Theory program.

Zenodo concept DOI: <https://doi.org/10.5281/zenodo.19056813>

## Status

This folder contains the public Skeleton/anomalous-dissipation paper files
synchronized from the local v1.4-candidate working state. The latest Zenodo
release is still v1.3. The DFC hierarchy is the empirical input; the separate
K41 variational-minimizer Paper-A files remain local until their release
strategy and Zenodo record are decided.

## Files

| File | Purpose |
|------|---------|
| `FST-TU_Turbulence_Skeleton_v1_en.tex` / `FST-TU_Turbulence_Skeleton_v1_en.pdf` | English paper source and PDF |
| `FST-TU_Turbulence_Skeleton_v1_de.tex` / `FST-TU_Turbulence_Skeleton_v1_de.pdf` | German paper source and PDF |
| `../../scripts/turbulence/compute_F_spectrum.py` | K41 free-energy minimizer and strict convexity test |
| `../../scripts/turbulence/compute_goy_shell_dfc.py` | Sabra/GOY shell-model DFC1/DFC2 verification |
| `../../scripts/turbulence/compute_F_spectrum.png` | Generated K41 free-energy result plot |
| `../../scripts/turbulence/compute_goy_shell_dfc.png` | Generated shell-model DFC result plot |

## Reproduce

From the repository root:

```bash
PYTHONIOENCODING=utf-8 python scripts/turbulence/compute_F_spectrum.py
PYTHONIOENCODING=utf-8 python scripts/turbulence/compute_goy_shell_dfc.py
```

The scripts write their PNG outputs next to the scripts in
`scripts/turbulence/`.

## Publication Gate

Internal proof notes, review chains, planning files, Zenodo credentials,
revision notes, and private split-strategy notes are intentionally not part of
this public package. They remain local-only until the project reaches the
required completion or submission gate.
