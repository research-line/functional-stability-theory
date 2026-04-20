# FST Spectrum Duality — Three-Lemma Endgame for Route B

**Status:** Draft v0.8 (conditional closure of CCM missing step MS2)
**Paper:** `paper/FST_SPECTRUM_DUALITY_v0-8_en.tex`

## Overview

This directory contains the FST axiom package for spectrum-zero identification
and the **Three-Lemma Endgame** that conditionally resolves the CCM (Connes--Consani--Moscovici)
missing step (MS2) via mass-based spectral cluster analysis.

### Key Result (Theorem, conditional)

Under a Galerkin-faithfulness hypothesis (numerically verified to 12 digits):

```
||(I - P_{V_lambda}) k_lambda|| <= r_lambda / g_* -> 0   as N -> infinity
```

Three independent lemmas + one standard inequality:
- **C2ca.1** (Scalar Secular Cancellation): |mu| ~ 3e-7
- **C2ca.2** (Projected Poisson Quasimode): ||h|| ~ 3e-6
- **C2ca.5** (Coercive Complement): g_* >= 5 (order one, lambda-independent)
- **C2ca.4** (Mass Concentration): COROLLARY of the above three

## Structure

```
paper/                  LaTeX paper (v0.8, Three-Lemma Endgame for Route B)
zookeeper/              Companion paper: numerical verification ("The Spectral Zookeeper")
scripts/                Numerical validation (Python, numpy/scipy)
results/                Server computation outputs (JSON)
proof-notes/            Formal proof documents (C2CA*.md)
```

## Scripts

| Script | Purpose |
|--------|---------|
| `c2bt_spectral_mass.py` | **Main diagnostic:** tolerance-free spectral mass analysis |
| `c2_quasimode_test.py` | Quasimode residual R0 computation |
| `c2_cancellation_lambda_scan.py` | Secular cancellation T1+T2 across lambda |
| `c2_poisson_decomposition.py` | Poisson kernel bulk/boundary decomposition |
| `c2_n_scaling.py` | Galerkin dimension N scaling of residuals |
| `b10_analytical_derivation.py` | Resolvent identity B10 verification |

## Dependencies

- Python 3.10+
- numpy, scipy, matplotlib (standard scientific stack)
- For server runs: SSH access to compute server (CCX13)

## Relation to FST Hierarchy

```
Ebene 0:  RH v2.1 (published, DOI 10.5281/zenodo.19546593)
Ebene 1:  THIS PAPER (conditional Route B endgame)    <-- you are here
          Math-Master "The Zeta Zoo" (classification)
          Physics-Master RFEP (physical instantiation)
Ebene 3:  Domain Supplements (Selberg, Dirichlet, BSD, ...)
```
