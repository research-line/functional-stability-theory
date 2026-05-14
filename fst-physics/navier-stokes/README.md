# Navier-Stokes Regularity and NS-LDI

Public reproducibility package for the Navier-Stokes domain supplements in the
Functional Stability Theory program.

Navier-Stokes concept DOI: <https://doi.org/10.5281/zenodo.19087449>
NS-LDI concept DOI: <https://doi.org/10.5281/zenodo.19056807>

## Status

This folder contains the public Navier-Stokes v2.2 paper files and
the NS-LDI v1.4 paper files. The latest Navier-Stokes Zenodo release is
v2.2 (record DOI <https://doi.org/10.5281/zenodo.20078143>), and the
latest NS-LDI Zenodo release is v1.4 (record DOI
<https://doi.org/10.5281/zenodo.20178115>). The main regularity result remains conditional on Assumption G2
(projection regularity); the NS-LDI package remains a proof-of-life and
diagnostic bridge, not an unconditional 3D Navier-Stokes proof.
The local research workspace already contains a stricter v2.3 candidate, but
this public repository intentionally stays on the latest DOI-backed Zenodo
release until the next upload is complete.

## Files

| File | Purpose |
|------|---------|
| `FST-NS_NavierStokes_Skeleton_v1_en.tex` / `FST-NS_NavierStokes_Skeleton_v1_en.pdf` | English Navier-Stokes paper source and PDF |
| `FST-NS_NavierStokes_Skeleton_v1_de.tex` / `FST-NS_NavierStokes_Skeleton_v1_de.pdf` | German Navier-Stokes paper source and PDF |
| `FST-NS_NavierStokes_Skeleton_v1_kombi.pdf` | Combined bilingual Navier-Stokes PDF |
| `FST-NS_LogDistanceIntegrability_v1_en.tex` / `FST-NS_LogDistanceIntegrability_v1_en.pdf` | English NS-LDI paper source and PDF |
| `FST-NS_LogDistanceIntegrability_v1_de.tex` / `FST-NS_LogDistanceIntegrability_v1_de.pdf` | German NS-LDI paper source and PDF |
| `FST-NS_LogDistanceIntegrability_v1_kombi.pdf` | Combined bilingual NS-LDI PDF generated from the public EN/DE pair |
| `../../scripts/navier-stokes/compute_ds3_lorenz.py` | DS3 stress test on the Lorenz attractor |
| `../../scripts/navier-stokes/compute_bv_selection.py` | Balanced-viscosity selection test on the Lorenz attractor |
| `../../scripts/navier-stokes/compute_bv_multi_attractor.py` | BV-selection stress test on Lorenz, Roessler, and Chen attractors |
| `../../scripts/navier-stokes/compute_mu_reach.py` | Measure-theoretic reach scan on Lorenz and KS attractors |
| `../../scripts/navier-stokes/compute_tll_ldi_lorenz.py` | TLL+LDI proof-of-life test on the Lorenz attractor |
| `../../scripts/navier-stokes/compute_tll_ldi_ks.py` | TLL+LDI diagnostics and grid refinement on the KS attractor |
| `../../scripts/navier-stokes/_results/` | Compact JSON and PNG result artifacts |

## Reproduce

From the repository root:

```bash
PYTHONIOENCODING=utf-8 python scripts/navier-stokes/compute_ds3_lorenz.py
PYTHONIOENCODING=utf-8 python scripts/navier-stokes/compute_bv_selection.py
PYTHONIOENCODING=utf-8 python scripts/navier-stokes/compute_bv_multi_attractor.py
PYTHONIOENCODING=utf-8 python scripts/navier-stokes/compute_mu_reach.py
PYTHONIOENCODING=utf-8 python scripts/navier-stokes/compute_tll_ldi_lorenz.py
PYTHONIOENCODING=utf-8 python scripts/navier-stokes/compute_tll_ldi_ks.py
```

The scripts write PNG and JSON outputs next to the scripts or under
`scripts/navier-stokes/_results/`.

## Publication Gate

Internal proof notes, review chains, planning files, Zenodo credentials,
revision notes, and private comparison notes are intentionally not part of this
public package. The repository contains the current public reproducibility
scripts and result artifacts for the Zenodo-backed Navier-Stokes v2.2 release
and the DOI-backed NS-LDI v1.4 package only. Local working notes and the
stricter Navier-Stokes v2.3 candidate remain local-only until the next release
gate is cleared.
