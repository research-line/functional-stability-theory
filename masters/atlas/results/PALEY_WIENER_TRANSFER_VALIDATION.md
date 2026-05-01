# Paley-Wiener Transfer Validation (Session 8)

**Date:** 2026-04-16
**Script:** `_scripts/paley_wiener_transfer_validation.py`
**Source:** `PALEY_WIENER_DIRICHLET_TRANSFER.md` §6.3

## Parameters

- lambda = 20000, L = 9.9035
- sigma_Phi^2 = 1.5337
- eps_L = L^(-1/2) = 0.3178
- max_m = 8, phi_width_cutoff = 6.0
- total walltime: 0.03 s

## Results per character

| chi | D | n_terms | A_diag | A_off | A_total | gap_emp | sign(A_off) | sign(A_total) |
|---|---|---|---:|---:|---:|---:|:---:|:---:|
| chi_5 | 5 | 286 | +8.2324 | -7.5140 | +0.7184 | +0.01603 | **FAIL** | OK |
| chi_8 | 8 | 282 | +8.2865 | -7.6950 | +0.5916 | -0.04902 | OK | **FAIL** |
| chi_12 | 12 | 276 | +8.0924 | -7.6204 | +0.4720 | +0.01154 | **FAIL** | OK |
| chi_13 | 13 | 288 | +8.2651 | -7.8290 | +0.4361 | -0.12504 | OK | **FAIL** |
| chi_17 | 17 | 288 | +8.2796 | -7.8116 | +0.4681 | +0.00886 | **FAIL** | OK |
| chi_21 | 21 | 281 | +8.0439 | -7.5250 | +0.5189 | -0.00423 | OK | **FAIL** |
| chi_24 | 24 | 276 | +8.0924 | -7.0463 | +1.0461 | +0.00866 | **FAIL** | OK |
| chi_29 | 29 | 288 | +8.3104 | -7.2316 | +1.0788 | +0.01439 | **FAIL** | OK |
| chi_33 | 33 | 281 | +8.0615 | -7.7127 | +0.3488 | -0.14221 | OK | **FAIL** |
| chi_60 | 60 | 272 | +7.8841 | -7.6036 | +0.2806 | +0.00521 | **FAIL** | OK |

## Aggregate

### All 10 characters
- sign(A_off) ok: **4/10**
- sign(A_total) ok: **6/10**
- Pearson R(A_off, gap): **+0.5113**, R^2 = 0.2614
- Pearson R(A_total, gap): **+0.4340**, R^2 = 0.1883

### Without chi_21 (known N-oscillator)
- sign(A_off) ok: **3/9**
- sign(A_total) ok: **6/9**
- Pearson R(A_off, gap): **+0.5099**, R^2 = 0.2600
- Pearson R(A_total, gap): **+0.4521**, R^2 = 0.2044

## Critical test

**chi_33 prediction:** A_off = -7.71268, gap_emp = -0.14221 -> **PASSED**

## Success criterion

Conjecture C2 (PALEY_WIENER_DIRICHLET_TRANSFER §6.2):
  sign(A_off) agreement >= 9/10 AND chi_33 negative.

**Result: NOT CONFIRMED** (sign_off = 4/10, chi_33 sign = OK)
