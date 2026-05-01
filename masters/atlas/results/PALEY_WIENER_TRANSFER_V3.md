# Paley-Wiener Transfer Validation v3

**Date:** 2026-04-17 (Session 11)
**Script:** `_scripts/paley_wiener_transfer_v3.py`
**Approach:** Full character-specific multiplier resolvent kernel.

## Motivation

v1 (Gauss kernel): 4/10, R²=0.03. v2 (ad-hoc pole-cos): 6/10, R²=0.31. v3 uses the **physically correct** Weil multiplier $\widehat K_\chi(\xi) = h_\mathrm{arch}(\xi) + 2\sum w_{p,m}\cos(m\log p\,\xi)$ and takes the character-specific resolvent $R_\chi(\xi) = 1/(\widehat K_\chi(\xi) - \mu_\chi + i\epsilon)$ whose Fourier inverse is the position kernel.

## Parameters

- lambda=20000, L=9.9035
- log_pm_max=9.0
- xi-grid: [-30.0, 30.0], n_xi=8192
- u-grid: [-18.0, 18.0], n_u=2048
- epsilon sweep: [0.1, 0.3, 0.5, 1.0, 2.0]

## Epsilon sweep results

| epsilon | sign(A_off) | R²(A_off) | sign(A_total) | R²(A_total) | sign_n21(A_off) | R²_n21(A_off) |
|---:|:---:|---:|:---:|---:|:---:|---:|
| 0.10 | 7/10 | 0.032 | 6/10 | 0.042 | 7/9 | 0.026 |
| 0.30 | 6/10 | 0.058 | 6/10 | 0.113 | 5/9 | 0.052 |
| 0.50 | 7/10 | 0.059 | 6/10 | 0.151 | 6/9 | 0.055 |
| 1.00 | 5/10 | 0.045 | 6/10 | 0.169 | 4/9 | 0.044 |
| 2.00 | 4/10 | 0.023 | 6/10 | 0.093 | 3/9 | 0.025 |

**Best epsilon**: 0.5: off 7/10, R²=0.059

## Per-character details (epsilon=0.5)

| chi | D | mu_chi | A_diag | A_off | A_total | gap_emp | sign |
|---|---|---:|---:|---:|---:|---:|:---:|
| chi_5 | 5 | -14.9205 | +32.0532 | -3.0932 | +28.9600 | +0.01603 | **FAIL** |
| chi_8 | 8 | -14.2620 | +31.8746 | -4.4814 | +27.3932 | -0.04902 | OK |
| chi_12 | 12 | -14.8678 | +30.0890 | +0.7685 | +30.8575 | +0.01154 | OK |
| chi_13 | 13 | -14.0012 | +33.2419 | -0.2272 | +33.0147 | -0.12504 | OK |
| chi_17 | 17 | -13.8631 | +31.0777 | -1.9585 | +29.1193 | +0.00886 | **FAIL** |
| chi_21 | 21 | -14.6210 | +27.6691 | -0.7815 | +26.8876 | -0.00423 | OK |
| chi_24 | 24 | -13.5967 | +33.1801 | +2.7364 | +35.9165 | +0.00866 | OK |
| chi_29 | 29 | -15.6311 | +25.5727 | +0.1705 | +25.7432 | +0.01439 | OK |
| chi_33 | 33 | -14.5332 | +28.6933 | +4.0729 | +32.7662 | -0.14221 | **FAIL** |
| chi_60 | 60 | -14.1791 | +26.6247 | +2.5144 | +29.1391 | +0.00521 | OK |
