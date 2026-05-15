# MOBIUS_GRADED_PRIMEHUB_TOY

**Purpose.** Finite check of the identity

`det(I-A_s) = Str(Lambda A_s) = sum mu(n)n^-s = prod_p(1-p^-s)`.

**Prime set:** `[2, 3, 5, 7, 11]`
**Squarefree terms:** `32`

| s | product determinant | exterior supertrace | squarefree Mobius sum | max abs error | reciprocal orientation |
|---|---:|---:|---:|---:|---:|
| s=2 | `0.621757463316` | `0.621757463316` | `0.621757463316` | `2.220e-16` | `1.60834418403` |
| s=3 | `0.832788808987` | `0.832788808987` | `0.832788808987` | `2.220e-16` | `1.20078462776` |
| s=0.75+5i | `1.25957583912-0.585044357749i` | `1.25957583912-0.585044357749i` | `1.25957583912-0.585044357749i` | `8.951e-16` | `0.653033226615+0.303319095832i` |

## Interpretation

The three left columns agree to floating-point precision. This is the finite
exterior-algebra core of the Möbius-graded Prime-Hub route.

The final column is the reciprocal orientation. It corresponds to the
ordinary positive Euler-product zeta direction. This is why the OP5
orientation must be audited before looking for eigenvalue-one events
at zeros of `zeta`: `1/zeta` has poles where `zeta` has zeros.

## First squarefree terms

| n | mu | support |
|---:|---:|---|
| 1 | 1 | `[]` |
| 2 | -1 | `[2]` |
| 3 | -1 | `[3]` |
| 5 | -1 | `[5]` |
| 6 | 1 | `[2, 3]` |
| 7 | -1 | `[7]` |
| 10 | 1 | `[2, 5]` |
| 11 | -1 | `[11]` |
| 14 | 1 | `[2, 7]` |
| 15 | 1 | `[3, 5]` |
| 21 | 1 | `[3, 7]` |
| 22 | 1 | `[2, 11]` |
| 30 | -1 | `[2, 3, 5]` |
| 33 | 1 | `[3, 11]` |
| 35 | 1 | `[5, 7]` |
| 42 | -1 | `[2, 3, 7]` |
