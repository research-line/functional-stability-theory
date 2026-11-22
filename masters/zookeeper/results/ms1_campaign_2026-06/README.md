# MS1 determinant / zero-counting campaign (June 2026)

This folder holds the measurement campaign that *From Landscape to Atlas*,
Part III (Zenodo v3.2, [10.5281/zenodo.21952620](https://doi.org/10.5281/zenodo.21952620)),
cites as the “full audit trail in the Zookeeper companion”. Scripts:

| Script | Purpose | Outputs in this folder |
|---|---|---|
| [`scripts/ms1_det_zero_counting_harness.py`](../../scripts/ms1_det_zero_counting_harness.py) | Argument-principle counting of the degree-2N secular polynomial against the completed Ξ on contours around zero packets | `MS1_DET_ZERO_COUNTING_*.json` (c = 9, 14) |
| [`scripts/ms1_newton_power_sum_balance.py`](../../scripts/ms1_newton_power_sum_balance.py) | Newton power-sum balance of the degree-paired normalisation | `MS1_NEWTON_BALANCE_*.json` (c = 9, 14, 25; N ≤ 40) |
| [`scripts/ms1_cartwright_ledger.py`](../../scripts/ms1_cartwright_ledger.py) | Cartwright/log-spread ledger over the Newton outputs | `MS1_CARTWRIGHT_LEDGER_2026-06-11.json`, `MS1_CARTWRIGHT_LEDGER_SCAN_2026-06-11.json` |

`gamma_zeros_cache.json` holds the zeta ordinates that the ledger reads.

All statements from this campaign are numerical. No MS1 or MS2 closure is
claimed (see Atlas Part III). The Zookeeper paper and its claim are
unchanged.

## Dependencies

Python 3 with `mpmath`. The harness and the Newton balance also import
`build_galerkin_matrix` from the external package **connes-cvs**
(Akiva Groskin, MIT licence, <https://github.com/akivag613/connes-cvs->, arXiv:2605.20224).
Point the environment variable `CCV_SRC` at a clone of that repository.

## Curation (2026-11-22)

The scripts are the versions used for the campaign. The only edits for
publication: paths point at this folder instead of the private working
directory, a host-specific search path was removed, and internal agent and
file references in comments and in the `purpose` metadata field were
neutralised. Numerics and parameters are unchanged.

Check on 2026-11-22: `ms1_cartwright_ledger.py` with the three N = 40, T = 800
Newton outputs reproduces the published ledger values exactly (margins and
log-spreads per block). The current script additionally writes the fields
`window` and `n_S`. The harness runs were not repeated (high-precision
compute).
