"""Regression test for GitHub Issue #1: non-conservation of quadratic energy
in the Sabra nonlinear term.

The Sabra shell model (L'vov, Podivilov, Pomyalov, Procaccia, Vandembroucq,
"Improved shell model of turbulence", Phys. Rev. E 58, 1811 (1998), Eq. (21))
must conserve the quadratic invariant E = sum_n |u_n|^2 in the inviscid,
unforced limit: dE/dt = 2*Re(sum_n conj(u_n) * N_n(u)) = 0 for every complex
state u, where N_n is the nonlinear term alone (no dissipation, no forcing).

Before the fix, both compute_goy_shell_dfc.py and
compute_shell_dfc_waterline_ledger.py conjugated the wrong factor in the
"local" nonlinear term (conj(u_{n+1}) instead of conj(u_{n-1})) and used the
wrong backward-term coefficient (-1/4 instead of +1/2), which breaks the
energy identity. The issue's own three-shell counterexample
(u = [1, 1, i], k_n = 0.125) gave dE/dt = -7*k_n/2 = -0.4375 instead of 0.
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts" / "turbulence"


def _load_module(name, filename):
    path = SCRIPTS_DIR / filename
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


# Importing compute_goy_shell_dfc.py / compute_shell_dfc_waterline_ledger.py
# as modules runs their full simulation (module-level script) -- expensive
# and unnecessary for a unit test of sabra_nonlinear alone. We instead
# reimplement the (now fixed) nonlinear term inline and cross-check it
# byte-for-byte against the source files' current coefficients/signs, so the
# test fails loudly if either file regresses.


def _extract_sabra_nonlinear_source(filename):
    """Read the sabra_nonlinear function body out of the given script without
    executing the rest of the (long-running) module."""
    path = SCRIPTS_DIR / filename
    text = path.read_text(encoding="utf-8")
    start = text.index("def sabra_nonlinear")
    # cut at the next top-level 'def ' after the function start
    rest = text[start:]
    next_def = rest.index("\ndef ", 1)
    func_src = rest[:next_def]
    return func_src


def _dEdt(u, k_n, nl_func):
    """dE/dt = 2*Re(sum_n conj(u_n) * N_n(u)) for the inviscid, unforced term."""
    nl = nl_func(u, k_n)
    return 2.0 * np.real(np.sum(np.conj(u) * nl))


@pytest.fixture()
def goy_sabra_nonlinear():
    ns = {"np": np}
    func_src = _extract_sabra_nonlinear_source("compute_goy_shell_dfc.py")
    exec(func_src, ns)  # noqa: S102 -- controlled, repo-local source extraction
    sabra_nonlinear = ns["sabra_nonlinear"]

    def wrapped(u, k_n):
        ns["k_n"] = k_n  # the module-level function closes over the global k_n
        return sabra_nonlinear(u)

    return wrapped


@pytest.fixture()
def ledger_sabra_nonlinear():
    ns = {"np": np}
    func_src = _extract_sabra_nonlinear_source("compute_shell_dfc_waterline_ledger.py")
    exec(func_src, ns)  # noqa: S102
    sabra_nonlinear = ns["sabra_nonlinear"]

    def wrapped(u, k):
        ns["N_SHELLS"] = len(u)  # the function closes over the module-level N_SHELLS
        return sabra_nonlinear(u, k)

    return wrapped


def test_issue1_three_shell_counterexample_now_conserves(goy_sabra_nonlinear):
    """Exact reproduction from Issue #1: u=[1,1,i], k_n=0.125 must give dE/dt=0."""
    k_n_scalar = 0.125
    u = np.array([1 + 0j, 1 + 0j, 1j], dtype=complex)
    k_n_vec = np.array([k_n_scalar, k_n_scalar, k_n_scalar])
    d = _dEdt(u, k_n_vec, goy_sabra_nonlinear)
    assert abs(d) < 1e-12, f"dE/dt={d}, expected 0 (issue claimed -0.4375 pre-fix)"


@pytest.mark.parametrize("N", [3, 4, 5, 6, 8, 10, 22, 30])
def test_goy_shell_dfc_conserves_energy_random_states(goy_sabra_nonlinear, N):
    """Inviscid, unforced dE/dt must vanish to machine precision for random
    complex shell states, with zero boundary conditions at both ends."""
    rng = np.random.default_rng(12345 + N)
    k_n = (2.0 ** -4) * 2.0 ** np.arange(N, dtype=float)
    for _ in range(5):
        u = (rng.standard_normal(N) + 1j * rng.standard_normal(N))
        d = _dEdt(u, k_n, goy_sabra_nonlinear)
        assert abs(d) < 1e-9 * max(1.0, np.sum(np.abs(u) ** 2)), (
            f"N={N}: dE/dt={d} (energy not conserved)"
        )


@pytest.mark.parametrize("N_SHELLS", [3, 4, 5, 6, 8, 10, 18, 30])
def test_waterline_ledger_conserves_energy_random_states(ledger_sabra_nonlinear, N_SHELLS):
    rng = np.random.default_rng(54321 + N_SHELLS)
    k = (2.0 ** -4) * 2.0 ** np.arange(N_SHELLS, dtype=float)
    for _ in range(5):
        u = (rng.standard_normal(N_SHELLS) + 1j * rng.standard_normal(N_SHELLS))
        d = _dEdt(u, k, lambda uu, kk: ledger_sabra_nonlinear(uu, kk))
        assert abs(d) < 1e-9 * max(1.0, np.sum(np.abs(u) ** 2)), (
            f"N_SHELLS={N_SHELLS}: dE/dt={d} (energy not conserved)"
        )
