# K41 Variational Minimiser

Public reproducibility package for Paper A in the turbulence split of the
Functional Stability Theory program.

Zenodo concept DOI: <https://doi.org/10.5281/zenodo.20131305>  
Latest v1.3 DOI: <https://doi.org/10.5281/zenodo.20562341>

Companion Turbulence/DFC Paper-B concept DOI:
<https://doi.org/10.5281/zenodo.19056813>

## Status

This folder contains the current public K41 variational-minimizer paper package:
the Kolmogorov spectrum is isolated as the unique global minimizer of a
scale-resolved free-energy functional under the stated joint minimization
problem. The latest live maintenance release is Zenodo v1.3. The conditional
anomalous-dissipation and Downhill Flux route lives in the companion folder
`../turbulence/`.

## Files

| File | Purpose |
|------|---------|
| `K41_Variational_Minimiser_v1_en.tex` / `K41_Variational_Minimiser_v1_en.pdf` | English paper source and PDF |
| `K41_Variational_Minimiser_v1_de.tex` / `K41_Variational_Minimiser_v1_de.pdf` | German paper source and PDF |
| `K41_Variational_Minimiser_v1_kombi.pdf` | Combined English/German PDF |
| `../../scripts/k41/compute_F_spectrum.py` | K41 free-energy minimizer and strict convexity test |
| `../../scripts/k41/compute_F_spectrum.png` | Generated K41 free-energy result plot |

## Reproduce

From the repository root:

```bash
PYTHONIOENCODING=utf-8 python scripts/k41/compute_F_spectrum.py
```

The script writes its PNG output next to the script in `scripts/k41/`.

## Publication Gate

Internal proof notes, planning files, Zenodo credentials, and review chains are
not part of this public package. They remain local-only until a publication or
submission gate explicitly opens them.
