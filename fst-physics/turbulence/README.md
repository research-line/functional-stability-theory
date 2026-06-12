# Turbulence Cascade

Public reproducibility package for the turbulence domain supplement in the
Functional Stability Theory program.

Zenodo concept DOI: <https://doi.org/10.5281/zenodo.19056813>

Companion K41 Paper-A concept DOI: <https://doi.org/10.5281/zenodo.20131305>

## Status

This folder contains the public Skeleton/anomalous-dissipation companion paper
files synchronized for the published v1.7 maintenance release:
<https://doi.org/10.5281/zenodo.20670703>. The DFC hierarchy is the
conditional input. Version 1.7 adds the synchronized bibliography/disclosure
maintenance, the initial Dual-DFC1 ledger, and the Shell-DFC waterline guardrail
while keeping DFC1^vee as an empirical/projection bridge rather than a proved
Navier-Stokes consequence.

The unconditional K41 variational-minimizer theorem is now separated into
`../k41-variational-minimiser/` and published under concept DOI
<https://doi.org/10.5281/zenodo.20131305> (latest v1.3:
<https://doi.org/10.5281/zenodo.20562341>).

## Files

| File | Purpose |
|------|---------|
| `FST-TU_Turbulence_Skeleton_v1_en.tex` / `FST-TU_Turbulence_Skeleton_v1_en.pdf` | English paper source and PDF |
| `FST-TU_Turbulence_Skeleton_v1_de.tex` / `FST-TU_Turbulence_Skeleton_v1_de.pdf` | German paper source and PDF |
| `FST-TU_Turbulence_Skeleton_v1_kombi.pdf` | Combined bilingual English/German PDF |
| `../../scripts/turbulence/compute_goy_shell_dfc.py` | Sabra/GOY shell-model DFC1/DFC2 verification |
| `../../scripts/turbulence/compute_goy_shell_dfc.png` | Generated shell-model DFC result plot |
| `../../scripts/turbulence/compute_dual_dfc.py` | Toy/ledger evaluator for the dual Downhill Flux Condition DFC1^vee |
| `../../scripts/turbulence/compute_shell_dfc_waterline_ledger.py` | Sabra-shell waterline smoke test with matched controls |
| `../../scripts/turbulence/_results/DUAL_DFC_LEDGER_2026-05-27.*` | Published snapshot of the initial Dual-DFC1 ledger |
| `../../scripts/turbulence/_results/SHELL_DFC_WATERLINE_LEDGER_2026-06-05.*` | Published snapshot of the Shell-DFC waterline smoke test |

## Reproduce

From the repository root:

```bash
PYTHONIOENCODING=utf-8 python scripts/turbulence/compute_goy_shell_dfc.py
PYTHONIOENCODING=utf-8 python scripts/turbulence/compute_dual_dfc.py
PYTHONIOENCODING=utf-8 python scripts/turbulence/compute_shell_dfc_waterline_ledger.py
```

The scripts write their generated outputs next to the scripts under
`scripts/turbulence/_results/`. The paper-cited local `_results/...` snapshots
are mirrored there for the public repository package.

## Publication Gate

Internal proof notes, review chains, planning files, Zenodo credentials,
revision notes, and private split-strategy notes are intentionally not part of
this public package. They remain local-only until the project reaches the
required completion or submission gate.
