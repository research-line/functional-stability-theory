# Paley-Wiener Transfer Validation v2 (Session 8 Iteration 2)

**Date:** 2026-04-16
**Script:** `_scripts/paley_wiener_transfer_v2.py`
**Approach:** Scan over kernel variants (Gaussian, pole, mixed)

## Parameters
- lambda=20000, L=9.9035, eps_L=0.3178
- log_pm_max=9.0 (primes cover up to exp(9.0))
- Total walltime: 25.65 s, 24 variants

## Variant scan results

| Variant | sign(A_off) | R2(A_off) | sign(A_total) | R2(A_total) |
|---|:---:|---:|:---:|---:|
| `gauss_sigma1.2` | 4/10 | 0.441 | 6/10 | 0.311 |
| `gauss_sigma2.0` | 4/10 | 0.145 | 6/10 | 0.084 |
| `gauss_sigma3.0` | 4/10 | 0.028 | 6/10 | 0.017 |
| `gauss_sigma4.0` | 4/10 | 0.007 | 6/10 | 0.004 |
| `gauss_sigma5.0` | 4/10 | 0.003 | 6/10 | 0.001 |
| `pole_hadamard` | 6/10 | 0.001 | 6/10 | 0.000 |
| `pole_hadamard2` | 6/10 | 0.047 | 6/10 | 0.032 |
| `pole_equal` | 6/10 | 0.185 | 6/10 | 0.118 |
| `mix_sigma2.0_alpha0.1_hadamard` | 4/10 | 0.005 | 6/10 | 0.000 |
| `mix_sigma2.0_alpha0.1_hadamard2` | 4/10 | 0.117 | 6/10 | 0.055 |
| `mix_sigma2.0_alpha0.5_hadamard` | 6/10 | 0.001 | 6/10 | 0.000 |
| `mix_sigma2.0_alpha0.5_hadamard2` | 5/10 | 0.058 | 6/10 | 0.036 |
| `mix_sigma2.0_alpha1.0_hadamard` | 6/10 | 0.001 | 6/10 | 0.000 |
| `mix_sigma2.0_alpha1.0_hadamard2` | 6/10 | 0.052 | 6/10 | 0.034 |
| `mix_sigma2.0_alpha3.0_hadamard` | 6/10 | 0.001 | 6/10 | 0.000 |
| `mix_sigma2.0_alpha3.0_hadamard2` | 5/10 | 0.049 | 6/10 | 0.033 |
| `mix_sigma3.0_alpha0.1_hadamard` | 5/10 | 0.002 | 6/10 | 0.000 |
| `mix_sigma3.0_alpha0.1_hadamard2` | 4/10 | 0.086 | 6/10 | 0.045 |
| `mix_sigma3.0_alpha0.5_hadamard` | 6/10 | 0.001 | 6/10 | 0.000 |
| `mix_sigma3.0_alpha0.5_hadamard2` | 6/10 | 0.054 | 6/10 | 0.035 |
| `mix_sigma3.0_alpha1.0_hadamard` | 6/10 | 0.001 | 6/10 | 0.000 |
| `mix_sigma3.0_alpha1.0_hadamard2` | 6/10 | 0.050 | 6/10 | 0.033 |
| `mix_sigma3.0_alpha3.0_hadamard` | 6/10 | 0.001 | 6/10 | 0.000 |
| `mix_sigma3.0_alpha3.0_hadamard2` | 6/10 | 0.048 | 6/10 | 0.033 |

## Best variants

- **Best A_off:** `pole_equal` - 6/10, R²=0.185
- **Best A_total:** `gauss_sigma1.2` - 6/10, R²=0.311

## Per-character results for best A_off: `pole_equal`

| chi | D | A_diag | A_off | A_total | gap_emp | sign_off_ok |
|---|---|---:|---:|---:|---:|:---:|
| chi_5 | 5 | +230.0572 | +131.7284 | +361.7856 | +0.01603 | OK |
| chi_8 | 8 | +269.5816 | +142.2936 | +411.8752 | -0.04902 | **FAIL** |
| chi_12 | 12 | +303.2662 | +222.4164 | +525.6826 | +0.01154 | OK |
| chi_13 | 13 | +345.9801 | +173.4798 | +519.4600 | -0.12504 | **FAIL** |
| chi_17 | 17 | +346.3978 | +93.2890 | +439.6868 | +0.00886 | OK |
| chi_21 | 21 | +377.5602 | +73.8939 | +451.4540 | -0.00423 | **FAIL** |
| chi_24 | 24 | +379.0827 | +158.4091 | +537.4918 | +0.00866 | OK |
| chi_29 | 29 | +424.4434 | +199.3004 | +623.7437 | +0.01439 | OK |
| chi_33 | 33 | +453.7475 | +255.2330 | +708.9805 | -0.14221 | **FAIL** |
| chi_60 | 60 | +484.3918 | +177.5348 | +661.9266 | +0.00521 | OK |
