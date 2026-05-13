# Turbulence Cascade

Public reproducibility package for the turbulence domain supplement in the
Functional Stability Theory program.

Zenodo concept DOI: <https://doi.org/10.5281/zenodo.19056813>

Companion K41 Paper-A concept DOI: <https://doi.org/10.5281/zenodo.20131305>

## Status

This folder contains the public Skeleton/anomalous-dissipation companion paper
files synchronized from the Zenodo v1.5 release. The latest Zenodo record is
<https://doi.org/10.5281/zenodo.20173281>. The DFC hierarchy is the conditional
input.

The unconditional K41 variational-minimizer theorem is now separated into
`../k41-variational-minimiser/` and published under concept DOI
<https://doi.org/10.5281/zenodo.20131305>.

## Files

| File | Purpose |
|------|---------|
| `FST-TU_Turbulence_Skeleton_v1_en.tex` / `FST-TU_Turbulence_Skeleton_v1_en.pdf` | English paper source and PDF |
| `FST-TU_Turbulence_Skeleton_v1_de.tex` / `FST-TU_Turbulence_Skeleton_v1_de.pdf` | German paper source and PDF |
| `FST-TU_Turbulence_Skeleton_v1_kombi.pdf` | Combined bilingual English/German PDF |
| `../../scripts/turbulence/compute_goy_shell_dfc.py` | Sabra/GOY shell-model DFC1/DFC2 verification |
| `../../scripts/turbulence/compute_goy_shell_dfc.png` | Generated shell-model DFC result plot |

## Reproduce

From the repository root:

```bash
PYTHONIOENCODING=utf-8 python scripts/turbulence/compute_goy_shell_dfc.py
```

The script writes its PNG output next to the script in `scripts/turbulence/`.

## Publication Gate

Internal proof notes, review chains, planning files, Zenodo credentials,
revision notes, and private split-strategy notes are intentionally not part of
this public package. They remain local-only until the project reaches the
required completion or submission gate.
