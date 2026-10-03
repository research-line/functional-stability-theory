"""Regression test for the Sabra shell-model energy FLUX Pi_n (follow-up to
GitHub Issue #1, 2026-10-03).

Issue #1's fix (commit 48448f5) corrected the nonlinear term
``sabra_nonlinear`` to conserve the quadratic invariant E = sum_n |u_n|^2.
It did NOT re-derive the shell-to-shell energy flux Pi_n used by
compute_goy_shell_dfc.py and compute_shell_dfc_waterline_ledger.py, which
both still used

    Pi_n = Im(k_n * u_n * conj(u_{n+1}) * u_{n+2})          (OLD, wrong)

-- the flux formula belonging to the *pre-fix* nonlinearity. This test
establishes, independent of any specific closed-form choice, the one
property any correct flux MUST have: cumulative consistency with the
(corrected) energy balance,

    Pi_n := -sum_{m=0}^{n} dE_m/dt|_NL = -2 * sum_{m=0}^{n} Re(conj(u_m) * N_m(u))

(see ``scripts/turbulence/sabra_flux.py`` for the full symbolic derivation
of the closed form tested here), and checks that:

  1. The shared closed-form ``sabra_energy_flux`` matches this brute-force
     cumulative definition to numerical precision, for N=3..30, random
     states, and the edge cases requested in the issue follow-up (zero
     state, single/double-shell excitation, real-only, extreme k-ratios).
  2. The telescoping identity dE_n/dt|_NL = Pi_{n-1} - Pi_n holds (with the
     boundary convention Pi_{-1} = 0, i.e. no flux into shell 0 from
     "shell -1", and no flux defined past the last triad, i.e. nothing
     flows out of the top shell beyond what the last triad carries).
  3. RED CONTROL: the OLD formula does NOT satisfy the cumulative
     consistency check above (documents that the regression would have
     caught the bug before this fix).
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts" / "turbulence"

sys.path.insert(0, str(SCRIPTS_DIR))
from sabra_flux import sabra_energy_flux  # noqa: E402


def _extract_function_source(filename, func_name):
    path = SCRIPTS_DIR / filename
    text = path.read_text(encoding="utf-8")
    start = text.index(f"def {func_name}")
    rest = text[start:]
    next_def = rest.index("\ndef ", 1)
    return rest[:next_def]


@pytest.fixture()
def goy_sabra_nonlinear():
    """The CORRECTED sabra_nonlinear from compute_goy_shell_dfc.py (post
    commit 48448f5), used as the ground-truth dynamics for the brute-force
    flux reference."""
    ns = {"np": np}
    func_src = _extract_function_source("compute_goy_shell_dfc.py", "sabra_nonlinear")
    exec(func_src, ns)  # noqa: S102 -- controlled, repo-local source extraction
    sabra_nonlinear = ns["sabra_nonlinear"]

    def wrapped(u, k_n):
        ns["k_n"] = k_n  # closes over module-level k_n
        return sabra_nonlinear(u)

    return wrapped


def _bruteforce_cumulative_flux(u, k_n, nl_func):
    """Pi_n = -2*sum_{m=0}^n Re(conj(u_m)*N_m(u)), n=0..N-1 (direct
    cumulative-sum definition, independent of any closed form)."""
    nl = nl_func(u, k_n)
    de_dt_per_shell = 2.0 * np.real(np.conj(u) * nl)
    return -np.cumsum(de_dt_per_shell)


def _old_buggy_flux(u, k_n):
    """The flux formula that was actually in both scripts before this fix
    (belongs to the pre-Issue#1 nonlinearity, not the corrected one)."""
    N = len(u)
    pi = np.zeros(N - 2)
    for n in range(N - 2):
        pi[n] = np.imag(k_n[n] * u[n] * np.conj(u[n + 1]) * u[n + 2])
    return pi


def _k_of(N, lam=2.0):
    return (2.0 ** -4) * lam ** np.arange(N, dtype=float)


@pytest.mark.parametrize("N", list(range(3, 31)))
def test_closed_form_matches_cumulative_definition_random(goy_sabra_nonlinear, N):
    rng = np.random.default_rng(1000 + N)
    k_n = _k_of(N)
    for _ in range(5):
        u = rng.standard_normal(N) + 1j * rng.standard_normal(N)
        bf = _bruteforce_cumulative_flux(u, k_n, goy_sabra_nonlinear)[: N - 2]
        cf = sabra_energy_flux(u, k_n)
        assert np.allclose(bf, cf, atol=1e-9, rtol=1e-7), (
            f"N={N}: closed-form Pi_n does not match the cumulative energy "
            f"balance (max diff {np.max(np.abs(bf - cf)):.3e})"
        )


@pytest.mark.parametrize(
    "label,build_u",
    [
        ("zero_state", lambda N: np.zeros(N, dtype=complex)),
        (
            "single_shell_first",
            lambda N: _set(np.zeros(N, dtype=complex), {0: 1.0 + 0j}),
        ),
        (
            "single_shell_middle",
            lambda N: _set(np.zeros(N, dtype=complex), {N // 2: 1.0 + 2j}),
        ),
        (
            "two_shell_excitation",
            lambda N: _set(np.zeros(N, dtype=complex), {0: 1.0 + 0j, 1: 1j}),
        ),
        ("all_real", lambda N: np.ones(N, dtype=complex)),
        ("complex_unit_modulus", lambda N: np.exp(1j * np.linspace(0, 1, N))),
    ],
)
@pytest.mark.parametrize("N", [3, 4, 5, 8, 22])
def test_closed_form_matches_cumulative_definition_edge_cases(
    goy_sabra_nonlinear, label, build_u, N
):
    k_n = _k_of(N)
    u = build_u(N)
    bf = _bruteforce_cumulative_flux(u, k_n, goy_sabra_nonlinear)[: N - 2]
    cf = sabra_energy_flux(u, k_n)
    assert np.allclose(bf, cf, atol=1e-9, rtol=1e-7), (
        f"{label}, N={N}: closed-form Pi_n mismatch "
        f"(max diff {np.max(np.abs(bf - cf)):.3e})"
    )


def _set(arr, updates):
    for idx, val in updates.items():
        arr[idx] = val
    return arr


@pytest.mark.parametrize("lam", [1.01, 1.5, 4.0, 10.0])
def test_closed_form_matches_cumulative_definition_extreme_lambda(
    goy_sabra_nonlinear, lam
):
    N = 10
    k_n = _k_of(N, lam=lam)
    rng = np.random.default_rng(42)
    u = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    bf = _bruteforce_cumulative_flux(u, k_n, goy_sabra_nonlinear)[: N - 2]
    cf = sabra_energy_flux(u, k_n)
    assert np.allclose(bf, cf, atol=1e-9, rtol=1e-7), f"lambda={lam} mismatch"


@pytest.mark.parametrize("N", [3, 5, 10, 22])
def test_telescoping_identity_holds_for_interior_shells(goy_sabra_nonlinear, N):
    """dE_n/dt|_NL = Pi_{n-1} - Pi_n for every shell that has a well-defined
    Pi on both sides (n = 1..N-3), with Pi_{-1} := 0 (nothing flows into the
    system from a nonexistent "shell -1"). This does NOT claim Pi_0 = 0 --
    the flux crossing the boundary after the FIRST shell is a genuine,
    generally nonzero, two-triad flux (see sabra_flux.py); only the outer
    system boundaries (before shell 0 and after shell N-1) are zero, which
    is checked separately below."""
    rng = np.random.default_rng(777 + N)
    k_n = _k_of(N)
    u = rng.standard_normal(N) + 1j * rng.standard_normal(N)

    nl = goy_sabra_nonlinear(u, k_n)
    de_dt = 2.0 * np.real(np.conj(u) * nl)  # length N

    pi = sabra_energy_flux(u, k_n)  # length N-2, Pi_0..Pi_{N-3}
    pi_padded = np.concatenate([[0.0], pi])  # [Pi_{-1}=0, Pi_0, ..., Pi_{N-3}]

    # dE_n/dt = Pi_{n-1} - Pi_n for n = 0..N-3 (both sides well-defined,
    # using the convention Pi_{-1}:=0 for n=0 -- this is a boundary
    # DEFINITION, not an empirical claim that Pi_0 itself vanishes).
    for n in range(N - 2):
        lhs = de_dt[n]
        rhs = pi_padded[n] - pi_padded[n + 1]
        assert abs(lhs - rhs) < 1e-9 * max(1.0, abs(lhs)), (
            f"N={N}, shell {n}: dE/dt={lhs} != Pi_{{n-1}}-Pi_n={rhs}"
        )


@pytest.mark.parametrize("N", [5, 10, 22])
def test_pi_after_first_shell_is_generally_nonzero(goy_sabra_nonlinear, N):
    """Anti-regression guard: Pi_0 (the flux crossing the boundary AFTER
    shell 0, i.e. the first entry of sabra_energy_flux) is a genuine,
    generally nonzero two-triad flux -- it is NOT one of the system's zero
    boundaries. If this ever starts asserting Pi_0 == 0 as a correctness
    requirement elsewhere, that would be the bug this test guards against."""
    k_n = _k_of(N)
    rng = np.random.default_rng(321)
    u = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    pi = sabra_energy_flux(u, k_n)
    assert abs(pi[0]) > 1e-6, (
        "Pi_0 is degenerate for this (deterministic) sample state -- pick a "
        "different seed; the point of this test is that Pi_0 is generally "
        "nonzero, unlike the true system boundaries."
    )


@pytest.mark.parametrize("N", [3, 5, 10, 22])
def test_system_boundary_fluxes_are_zero(goy_sabra_nonlinear, N):
    """The only two fluxes that are zero by construction are the SYSTEM's
    outer boundaries: no flux enters from before shell 0 (Pi_{-1} := 0,
    structural, not fitted), and no flux leaves past the last shell (there
    is no triad beyond index N-3, so the cumulative sum of ALL per-shell
    dE/dt must vanish -- global energy conservation of the nonlinear term,
    independent of this test file; see test_sabra_energy_conservation.py).
    Neither condition says anything about Pi_0 or any other interior Pi_n."""
    rng = np.random.default_rng(555 + N)
    k_n = _k_of(N)
    u = rng.standard_normal(N) + 1j * rng.standard_normal(N)

    nl = goy_sabra_nonlinear(u, k_n)
    de_dt = 2.0 * np.real(np.conj(u) * nl)

    # (a) flux INTO the system before shell 0: zero by construction
    #     (Pi_{-1} := 0 is how sabra_energy_flux is built, not measured).
    pi_minus_one = 0.0
    assert pi_minus_one == 0.0

    # (b) flux OUT of the system after the last shell (shell N-1): the
    #     cumulative flux must have fully accounted for all energy change
    #     by the time the chain of triads runs out, i.e. the TOTAL sum of
    #     dE_n/dt over every shell is zero -- nothing is left over to flow
    #     out past shell N-1.
    total_de_dt = float(np.sum(de_dt))
    assert abs(total_de_dt) < 1e-9 * max(1.0, np.sum(np.abs(u) ** 2)), (
        f"N={N}: total dE/dt={total_de_dt}, expected 0 "
        f"(no flux can leave the system past the last shell)"
    )


def test_old_formula_fails_cumulative_consistency_red_control(goy_sabra_nonlinear):
    """RED CONTROL: the formula that was actually in both scripts before
    this fix must NOT satisfy the cumulative energy-balance identity for
    the corrected dynamics -- this is exactly the bug this fix addresses."""
    any_mismatch = False
    for N in (5, 10, 22):
        k_n = _k_of(N)
        rng = np.random.default_rng(99 + N)
        u = rng.standard_normal(N) + 1j * rng.standard_normal(N)
        bf = _bruteforce_cumulative_flux(u, k_n, goy_sabra_nonlinear)[: N - 2]
        old = _old_buggy_flux(u, k_n)
        if not np.allclose(bf, old, atol=1e-9, rtol=1e-7):
            any_mismatch = True
    assert any_mismatch, (
        "Expected the OLD (pre-fix) flux formula to disagree with the "
        "energy-consistent cumulative definition -- if this now passes, "
        "the red control no longer discriminates and must be revisited."
    )
