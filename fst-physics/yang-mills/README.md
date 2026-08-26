# Yang-Mills Mass Gap

Public reproducibility package for the Yang-Mills domain supplement in the
Functional Stability Theory program.

Zenodo concept DOI: <https://doi.org/10.5281/zenodo.19087433>
Latest Zenodo v2.6 record DOI: <https://doi.org/10.5281/zenodo.20716608>

## Status

This folder contains the public v2.6 release package synchronized from the
local strict-review, design-check, and citation-checked working state. The
Zenodo v2.6 record is live under the record DOI above. The continuum mass-gap
step remains conditional on the analytical renormalization-group contraction
input.

## Files

| File | Purpose |
|------|---------|
| `FST-YM_YangMills_MassGap_en.tex` / `FST-YM_YangMills_MassGap_en.pdf` | English paper source and PDF |
| `FST-YM_YangMills_MassGap_de.tex` / `FST-YM_YangMills_MassGap_de.pdf` | German paper source and PDF |
| `FST-YM_YangMills_MassGap_kombi.pdf` | Combined bilingual PDF |
| `../../scripts/yang-mills/compute_dobrushin_su2.py` | SU(2) lattice Dobrushin influence scan |
| `../../scripts/yang-mills/compute_birkhoff_rg.py` | Birkhoff scan with finite good-scale density, bad-defect sums, and saturation-zone diagnostics |
| `../../scripts/yang-mills/compute_casimir_corridor_ledger.py` | Predefined gauge-sector Casimir-window audit with commutator and matched post-hoc controls |
| `../../scripts/yang-mills/compute_coercive_complement_ledger.py` | Finite self-adjoint outer-gap and residual/gap gate with bad-scale and Gribov matched controls |
| `../../scripts/yang-mills/compute_orbit_rigidity_ledger.py` | Hodge-inspired bounded-component, gauge-origin, residue-gap, and limit-closure audit |
| `../../scripts/yang-mills/compute_os_capacity_ledger.py` | OS-danger capacity ledger with pre-registered RG-window Tide-Clock fields and negative-control diagnostics |
| `../../scripts/yang-mills/compute_rp_os_transfer_ledger.py` | RP/OS transfer matrix positivity ledger |
| `../../scripts/yang-mills/compute_rp_os_rfep_transfer_ledger.py` | RFEP transfer matrix diagnostic ledger |
| `../../scripts/yang-mills/compute_u1_2d_strong_coupling_positive_control.py` | 2D U(1) strong-coupling positive control ledger |
| `../../scripts/yang-mills/compute_ym_waisen_transfer_ledger.py` | Yang-Mills waisen transfer verification ledger |
| `../../scripts/yang-mills/u1_strong_coupling_positive_control.py` | Direct U(1) character expansion positive control script |
| `../../scripts/yang-mills/compute_dobrushin_su2.png` | Generated Dobrushin result plot |
| `../../scripts/yang-mills/compute_birkhoff_rg.png` | Generated Birkhoff/RG result plot |

## Reproduce

From the repository root:

```bash
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_dobrushin_su2.py
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_birkhoff_rg.py
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_casimir_corridor_ledger.py
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_coercive_complement_ledger.py
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_orbit_rigidity_ledger.py
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_os_capacity_ledger.py
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_rp_os_transfer_ledger.py
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_rp_os_rfep_transfer_ledger.py
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_u1_2d_strong_coupling_positive_control.py
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_ym_waisen_transfer_ledger.py
PYTHONIOENCODING=utf-8 python scripts/yang-mills/u1_strong_coupling_positive_control.py
```

For a bounded diagnostic run that does not overwrite the tracked plot:

```bash
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_birkhoff_rg.py --no-plot --epsilon 0.01 --n-levels 6
```

The OS-capacity ledger can likewise write an isolated, disposable control run:

```bash
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_os_capacity_ledger.py --rows 20 --date-tag smoke --data-dir /tmp/ym-os-data --output-dir /tmp/ym-os-results
```

The coercive-complement successor can likewise write an isolated ledger:

```bash
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_coercive_complement_ledger.py --output-dir /tmp/ym-coercive-complement
```

It computes `s_lambda`, `p_lambda`, an explicit outer `g_*`, and
`(s_lambda+p_lambda)/g_*` for one finite positive control and matched bad-scale
and Gribov-complement negatives. A numerically identical post-hoc-cluster
control is rejected separately as circular. Self-adjointness, orthogonal
projection, and the reducing-subspace identity are fail-closed gates.
`transfer_decision` becomes review-eligible only when target cluster and
complement were fixed before leakage measurement and an independently verified
physical Yang--Mills predefinition certificate is supplied. All bundled
fixtures are synthetic, no bundled row is transfer-eligible, and `claim_pass`
remains zero.

The Selberg-inspired Casimir-corridor audit can write a separate disposable
ledger:

```bash
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_casimir_corridor_ledger.py --output-dir /tmp/ym-casimir-corridor
```

It derives a Casimir-window projector inside a declared physical/Gauss-law
sector and checks separately that the interval was fixed before spectral or
commutator inspection, that both Hamiltonian and Casimir preserve the gauge
sector, and that the Hamiltonian reduces the derived corridor. The matched
post-hoc control uses the same matrices and interval as the positive fixture
but is rejected on provenance. A second matched control is gauge invariant yet
mixes the corridor with its physical exterior, showing why gauge invariance
alone is insufficient. The finite compact-U(1) strong-coupling control uses the
gauge-reduced loop-flux basis with electric Casimir `n^2`; its `n=±1` corridor
has a unit spectral separation from the retained exterior. A reducing matched
negative with a collapsed Hamiltonian separation still fails. The corridor
provenance is explicitly gauge/Casimir-only: neither RG coercivity nor numerical
residual smallness is used to define it. The U(1) row passes only as an analytic
nonphysical calibration. None of these finite controls supplies an independently
verified physical non-Abelian Yang--Mills corridor, and `claim_pass=0`.

The Hodge-inspired orbit-rigidity successor can write an isolated ledger:

```bash
PYTHONIOENCODING=utf-8 python scripts/yang-mills/compute_orbit_rigidity_ledger.py --output-dir /tmp/ym-orbit-rigidity
```

It records `bounded_scale_component`, `gauge_orbit_certificate`,
`continuum_residue_gap`, and `bad_component_escape` for one predeclared finite
component family and matched negatives. Perfect good-scale density cannot pay
for unbounded component complexity; a gauge-orbit certificate derived from the
successful gap signal is rejected as circular; residual-gap collapse and new
remainder modes fail separately. A correct finite Dirichlet-form test is kept
as diagnostic when independent physical origin or Mosco/compact limit closure
is missing. All bundled rows remain finite controls and `claim_pass=0`.

The OS-capacity ledger's window CSV and JSON expose `rg_window_id`,
`window_predefined`, scale occupancy, safe-signal and OS-capacity shares,
nonlocal defect concentration,
bad-run switches, an independently supplied alternate blocking-control ratio,
nonlocal tail cost, and a fail-closed `transfer_status`. The generated windows
depend only on scale index and run length; external CSV inputs must declare the
window and pre-registration flag explicitly. These finite-window diagnostics
do not establish OS compactness or a continuum mass gap.

The generated controls include a dedicated Bad-Channel false positive: its
mean log contraction is negative and its aggregate is labelled summable, so
the base ledger accepts it, while a pre-registered one-percent half-scale
window carries almost all OS-danger and nonlocal defect mass. The separate
`transfer_decision` must therefore reject it as
`rejected_bad_channel_false_positive`; the strong-coupling and summable
positive controls remain `control_pass_windows_clear_no_claim`.

For a chosen margin `epsilon`, the script reports the observed fraction of
scales satisfying `tau_B <= 1 - epsilon`, the cumulative excess
`sum_k max(0, tau_B(R_k) - (1 - epsilon))`, and maximal contiguous saturation
zones. These are finite-sequence diagnostics only; they do not establish
Kingman hypotheses, scale-uniform coercivity, OS compactness, or a continuum
Yang--Mills mass gap.

The plotting scripts write their PNG outputs next to the scripts in
`scripts/yang-mills/`. The ledger script writes local `_data/` and `_results/`
folders under `scripts/yang-mills/`; those generated working outputs are not
versioned here.

## Publication Gate

Internal proof notes, `_results/`, design-check directories, extracted text
artefacts, review chains, planning files, Zenodo credentials, and private
comparison notes are intentionally not part of this public package. They remain
local-only until the project reaches the required completion gate.
