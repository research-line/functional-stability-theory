# Beal Height Dominance — published corrective v1.1

*A Conditional Height-Dominance Architecture for Beal’s Conjecture: A Zookeeper-Style Reduction with a Power-Residue Sieve Companion*

This curated companion contains manuscript sources and historical numerical diagnostics. The published predecessor is [Zenodo v1.0](https://doi.org/10.5281/zenodo.21916518); the [concept DOI](https://doi.org/10.5281/zenodo.21916517) follows subsequent versions. The corrective [v1.1 record](https://zenodo.org/records/23091637) is published (DOI [10.5281/zenodo.23091637](https://doi.org/10.5281/zenodo.23091637)). Both public PDF downloads were checked against the manuscript hashes; the description, latest-version link and DataCite DOI identity were verified.

The Beal conjecture and the HD3′ arithmetic bridge remain open. Elementary normalizations, formal asymptotic sieve inputs, historical finite diagnostics, conditional reductions, and source-reported external results have distinct scopes. Local residue passes do not establish a global root, a common divisor, statistical significance, or a numerical height threshold. The proposed exponential bound on cumulative positive exact counts is of Beal strength and remains unproved.

The corrected manuscripts restrict the asymptotic sieve to fixed K ≥ 4, sufficiently large H, and the stated polylogarithmic-to-height Q range. They correct zero-residue divisibility, genus assumptions, resonant character sums, scalar-family monodromy, and the scope of the recent Chocian, Sahoo, and Pasten references. External proof/replay validation of those preprints is not supplied by this package.

## Contents and provenance

- `paper/`: EN/DE TeX and PDFs, with the two canonical disclosure inputs.
- `evidence/`: eight original historical diagnostic producers and 17 original CSV files.
- `manifest.json`: SHA-256, byte lengths, CSV schemas and row counts for the curated artifacts.
- `tests/test_results.py`: offline integrity, row-alignment and finite-diagnostic checks using the Python standard library.

The CSVs and producer scripts are byte-for-byte copies of the local research archives. Dates embedded in their filenames identify historical runs. Some producer comments and generated text contain historical interpretations; the corrected manuscripts and the claim boundaries above govern their current use. Neither copying an artifact nor passing integrity tests certifies a mathematical proof. The near-power producer uses probable-prime factorization checks; its factor columns are historical diagnostic output, not independently certified primality evidence. Original private proof/review notes, operating files, downloaded third-party PDFs, logs, and credentials are excluded.

The near-miss corpus has 5000 rows, ordered by relative gap. Each variable-Q table uses its first 1000 rows. The three May root ledgers and four June survivor ledgers each have 750 rows: 250 actual near-misses, 250 matched random-root controls, and 250 exact-power controls. The rank-8 table contains seven nontrivial prime channels. The above-threshold and S-unit-shadow tables contain five selected rows each; the synthetic table contains 20 control rows. These nested or selected outputs are not independent random samples.

The S-unit and synthetic-control CSVs are curated historical exports whose filenames identify their June runs. Their original JSON envelopes and private reports remain internal. Historical source-fit classifications do not validate the proof of any cited paper or transfer its theorem to Beal.

## Offline validation

From this directory:

```powershell
python -m unittest discover -s tests -v
```

The tests verify manifest identities, arithmetic row consistency, input/output alignment, the rank-8 q=37 channel, and the stated control labels. They do not rerun the full sieve sweep or upgrade historical claims.

## Replaying diagnostics

Run replays in a copy of this directory: several historical producers default to overwriting their output tables. The following commands give the recorded setup, rather than evidence that a new full computation was performed for v1.1:

```powershell
cd evidence
python scan_near_power_fractures.py --max-base 250 --min-exp 3 --max-exp 7 --keep 5000 --max-rel-gap 1e-4
python power_residue_sieve.py --prime-bound 251 --limit 5000
python run_variable_Q_sieve.py --input near_power_b250_e3-7_z3-9_k5000.csv --limit 1000 --primes 251 503 1009 2003
python root_branch_drift_ledger.py --limit 250 --primes 251 503 1009
python root_branch_drift_ledger.py --limit 250 --primes 7 13 31 61 --out-prefix ROOT_BRANCH_SURVIVOR_SCAN_2026-06-01
python root_branch_rank8_micro_audit.py --rank 8 --max-q 61
python above_threshold_survivor_search.py --limit 5000
python s_unit_shadow_audit.py --help
python q97_plus_synthetic_control.py --help
```

The final two scripts expose output-path options via `--help`; supply paths within your copied directory when replaying. The original full-sweep timing is historical and varies with hardware.

## Building the papers

From `paper/`, run `pdflatex -interaction=nonstopmode -halt-on-error BEAL_Height_Dominance_v1_1_en.tex` and the corresponding `_de.tex` command three times each. Both disclosure inputs are included. A current TeX distribution with the declared standard packages is required.

Author: Lukas Geiger. License: [CC BY 4.0](../../LICENSE), consistent with the existing research repository and predecessor Zenodo record.
