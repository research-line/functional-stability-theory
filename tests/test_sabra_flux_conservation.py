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
     state, single/double-shell excitation, real-only, extreme k-ratios)
     -- over ALL N-1 inner boundaries it returns, not a truncated prefix.
  2. The telescoping identity dE_n/dt|_NL = Pi_{n-1} - Pi_n holds for
     EVERY shell n = 0, ..., N-1 (with Pi_{-1} = Pi_{N-1} = 0 the only two
     system-boundary conventions).
  3. RED CONTROL: the OLD formula does NOT satisfy the cumulative
     consistency check above (documents that the regression would have
     caught the bug before this fix).
  4. RED CONTROL: an energy-conserving but WRONG coefficient choice
     (a, b, c) = (1, -0.3, -0.7), which still satisfies a+b+c=0 and thus
     passes a pure energy-conservation test, is caught by the SECOND
     invariant H = sum_n (-2)^n |u_n|^2 (Eq. "Hinv" of the primary
     source, valid for this model's standard (a,b,c)=(1,-1/2,-1/2)).
  5. A literal line-by-line transcription of the primary source's Eq. (21)
     (L'vov, Podivilov, Pomyalov, Procaccia, Vandembroucq, Phys. Rev. E 58,
     1811 (1998)) -- using the PAPER's own (k_{n+1}, k_n, k_{n-1}) index
     convention, not the code's (k_n, k_{n-1}, k_{n-2}) -- matches the
     code's ``sabra_nonlinear`` exactly under k -> k/LAM (see sabra_flux.py
     module docstring, "Index convention vs. the primary source").
"""

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


@pytest.fixture()
def ledger_sabra_nonlinear():
    """The CORRECTED sabra_nonlinear from
    compute_shell_dfc_waterline_ledger.py (post commit 48448f5) -- the
    OTHER actual production function, alongside goy_sabra_nonlinear.
    Required by the 2026-10-03 Endabnahme (A3): dynamics/flux invariants
    must be anchored to BOTH shipped production functions, not only to
    the GOY script or to a helper reimplementation such as
    ``_paper_literal_nonlinear``. Unlike the GOY version, this function
    takes ``k`` as an explicit argument but still closes over the
    module-level ``N_SHELLS`` for its slice bounds."""
    ns = {"np": np}
    func_src = _extract_function_source(
        "compute_shell_dfc_waterline_ledger.py", "sabra_nonlinear"
    )
    exec(func_src, ns)  # noqa: S102 -- controlled, repo-local source extraction
    sabra_nonlinear = ns["sabra_nonlinear"]

    def wrapped(u, k_n):
        ns["N_SHELLS"] = len(u)  # closes over module-level N_SHELLS
        return sabra_nonlinear(u, k_n)

    return wrapped


def _paper_literal_nonlinear(u, k_n, a=1.0, b=-0.5, c=-0.5):
    """Literal, line-by-line transcription of Eq. (21) of L'vov, Podivilov,
    Pomyalov, Procaccia, Vandembroucq, "Improved shell model of
    turbulence", Phys. Rev. E 58, 1811 (1998):

        du_n/dt = i*(a*k_{n+1}*u_{n+2}*conj(u_{n+1})
                      + b*k_n*u_{n+1}*conj(u_{n-1})
                      - c*k_{n-1}*u_{n-1}*u_{n-2})

    using the PAPER's own index convention (k_{n+1}, k_n, k_{n-1}), not the
    code's (k_n, k_{n-1}, k_{n-2}). Boundary convention u_{-2}=u_{-1}=
    u_N=u_{N+1}=0."""
    N = len(u)
    nl = np.zeros(N, dtype=complex)
    for n in range(N):
        u_np2 = u[n + 2] if n + 2 < N else 0j
        u_np1 = u[n + 1] if n + 1 < N else 0j
        u_nm1 = u[n - 1] if n - 1 >= 0 else 0j
        u_nm2 = u[n - 2] if n - 2 >= 0 else 0j
        k_np1 = k_n[n + 1] if n + 1 < len(k_n) else k_n[n] * 2.0  # LAM=2 geometric
        term1 = a * k_np1 * u_np2 * np.conj(u_np1)
        term2 = b * k_n[n] * u_np1 * np.conj(u_nm1)
        term3 = -c * (k_n[n - 1] if n - 1 >= 0 else k_n[0] / 2.0) * u_nm1 * u_nm2
        nl[n] = 1j * (term1 + term2 + term3)
    return nl


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


def _energy_conserving_wrong_coeff_nonlinear(u, k_n):
    """Negative control (2026-10-03, independent review follow-up):
    (a, b, c) = (1, -0.3, -0.7) still satisfies a+b+c=0 (so dE/dt=0 to
    machine precision -- a pure energy test alone would PASS this), but it
    is NOT the physical model (a, b, c) = (1, -1/2, -1/2) and must be
    caught by a second, independent invariant (see
    ``test_h_invariant_catches_energy_conserving_wrong_coefficients``)."""
    N = len(u)
    a, b, c = 1.0, -0.3, -0.7
    nl = np.zeros(N, dtype=complex)
    nl[: N - 2] += a * k_n[: N - 2] * np.conj(u[1 : N - 1]) * u[2:N]
    nl[1 : N - 1] += b * k_n[: N - 2] * u[2:N] * np.conj(u[: N - 2])
    nl[2:N] += -c * k_n[: N - 2] * u[1 : N - 1] * u[: N - 2]
    return 1j * nl


def _k_of(N, lam=2.0):
    return (2.0 ** -4) * lam ** np.arange(N, dtype=float)


@pytest.mark.parametrize("N", list(range(3, 31)))
def test_closed_form_matches_cumulative_definition_random(goy_sabra_nonlinear, N):
    rng = np.random.default_rng(1000 + N)
    k_n = _k_of(N)
    for _ in range(5):
        u = rng.standard_normal(N) + 1j * rng.standard_normal(N)
        bf = _bruteforce_cumulative_flux(u, k_n, goy_sabra_nonlinear)[: N - 1]
        cf = sabra_energy_flux(u, k_n)
        assert np.allclose(bf, cf, atol=1e-9, rtol=1e-7), (
            f"N={N}: closed-form Pi_n does not match the cumulative energy "
            f"balance over ALL {N - 1} inner boundaries "
            f"(max diff {np.max(np.abs(bf - cf)):.3e})"
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
    bf = _bruteforce_cumulative_flux(u, k_n, goy_sabra_nonlinear)[: N - 1]
    cf = sabra_energy_flux(u, k_n)
    assert np.allclose(bf, cf, atol=1e-9, rtol=1e-7), (
        f"{label}, N={N}: closed-form Pi_n mismatch over all N-1 boundaries "
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
    bf = _bruteforce_cumulative_flux(u, k_n, goy_sabra_nonlinear)[: N - 1]
    cf = sabra_energy_flux(u, k_n)
    assert np.allclose(bf, cf, atol=1e-9, rtol=1e-7), f"lambda={lam} mismatch"


@pytest.mark.parametrize("N", [3, 5, 10, 22])
def test_telescoping_identity_holds_for_every_inner_boundary(goy_sabra_nonlinear, N):
    """dE_n/dt|_NL = Pi_{n-1} - Pi_n for EVERY shell n = 0, ..., N-1 (not
    just the interior ones), using the two SYSTEM-boundary conventions
    Pi_{-1} := 0 (nothing flows into the system before shell 0) and
    Pi_{N-1} := 0 (nothing flows out past the last shell -- equivalent to
    global energy conservation of the nonlinear term). This does NOT claim
    Pi_0 = 0, nor Pi_{N-2} = 0 -- both are genuine, generally nonzero,
    interior fluxes; only the two system boundaries vanish."""
    rng = np.random.default_rng(777 + N)
    k_n = _k_of(N)
    u = rng.standard_normal(N) + 1j * rng.standard_normal(N)

    nl = goy_sabra_nonlinear(u, k_n)
    de_dt = 2.0 * np.real(np.conj(u) * nl)  # length N

    pi = sabra_energy_flux(u, k_n)  # length N-1, Pi_0..Pi_{N-2}
    # Pad with the two true system-boundary zeros: [Pi_{-1}=0, Pi_0, ...,
    # Pi_{N-2}, Pi_{N-1}=0].
    pi_padded = np.concatenate([[0.0], pi, [0.0]])

    for n in range(N):
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


@pytest.mark.parametrize("N", [5, 10, 22])
def test_pi_before_last_shell_is_generally_nonzero(goy_sabra_nonlinear, N):
    """Anti-regression guard for the Fable-review boundary fix (2026-10-03):
    Pi_{N-2} (the LAST entry of sabra_energy_flux, i.e. the flux INTO the
    final shell) is a genuine, generally nonzero flux -- it is NOT the
    system boundary Pi_{N-1} (which is zero by construction and not part
    of the returned array). Before the fix this entry was silently
    zero-padded by the two calling scripts; this test would have caught
    that regression."""
    k_n = _k_of(N)
    rng = np.random.default_rng(654)
    u = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    pi = sabra_energy_flux(u, k_n)
    assert abs(pi[-1]) > 1e-6, (
        "Pi_{N-2} is degenerate for this (deterministic) sample state -- "
        "pick a different seed; the point of this test is that Pi_{N-2} "
        "(the last returned entry) is generally nonzero."
    )
    # And it must actually equal dE_{N-1}/dt|_NL (the flux INTO the last
    # shell), not the structural zero Pi_{N-1}.
    nl = goy_sabra_nonlinear(u, k_n)
    de_dt_last = 2.0 * np.real(np.conj(u[-1]) * nl[-1])
    assert abs(pi[-1] - de_dt_last) < 1e-9 * max(1.0, abs(de_dt_last)), (
        f"N={N}: Pi_{{N-2}}={pi[-1]} should equal dE_{{N-1}}/dt|_NL={de_dt_last}"
    )


@pytest.mark.parametrize("N", [3, 5, 10, 22])
def test_system_boundary_fluxes_are_zero(goy_sabra_nonlinear, N):
    """The only two fluxes that are zero by construction are the SYSTEM's
    outer boundaries: no flux enters from before shell 0 (Pi_{-1} := 0,
    structural, not fitted), and no flux leaves past the last shell (there
    is no triad beyond index N-3, so the cumulative sum of ALL per-shell
    dE/dt must vanish -- global energy conservation of the nonlinear term,
    independent of this test file; see test_sabra_energy_conservation.py).
    Neither condition says anything about Pi_0, Pi_{N-2}, or any other
    *interior* Pi_n -- those are checked elsewhere. NOTE (Endabnahme
    2026-10-03, A3): check (a) below is a DEFINITIONAL/structural fact
    about an empty prefix sum, not an independent measurement -- it is
    recorded for documentation, not claimed as a correctness test. Check
    (c) is the genuine flux-anchored last-boundary comparison; it is the
    one that actually uses the closed-form ``sabra_energy_flux`` output,
    unlike a bare re-sum of the same per-shell derivatives would."""
    rng = np.random.default_rng(555 + N)
    k_n = _k_of(N)
    u = rng.standard_normal(N) + 1j * rng.standard_normal(N)

    nl = goy_sabra_nonlinear(u, k_n)
    de_dt = 2.0 * np.real(np.conj(u) * nl)

    # (a) flux INTO the system before shell 0: there is nothing to sum
    #     over before shell 0 -- this is definitional (an empty sum is
    #     0.0 by construction), not a physical measurement.
    empty_prefix_sum = float(np.sum(de_dt[:0]))
    assert empty_prefix_sum == 0.0

    # (b) flux OUT of the system after the last shell (shell N-1): the
    #     cumulative flux must have fully accounted for all energy change
    #     by the time the chain of triads runs out, i.e. the TOTAL sum of
    #     dE_n/dt over every shell is zero -- nothing is left over to flow
    #     out past shell N-1. This is an actual measurement (full-sum
    #     reduction), not a hand-written literal.
    total_de_dt = float(np.sum(de_dt))
    assert abs(total_de_dt) < 1e-9 * max(1.0, np.sum(np.abs(u) ** 2)), (
        f"N={N}: total dE/dt={total_de_dt}, expected 0 "
        f"(no flux can leave the system past the last shell)"
    )

    # (c) Real last-boundary cross-check that actually USES the closed-form
    #     flux (fixes the previous tautological re-sum flagged by the
    #     2026-10-03 Endabnahme, A3): with the telescoping identity
    #     T_{N-1} = Pi_{N-2} - Pi_{N-1} and the structural Pi_{N-1}=0, the
    #     quantity -Pi_{N-2} + T_{N-1} must equal the independently
    #     measured sum_n T_n (= total_de_dt from (b)) -- this genuinely
    #     depends on sabra_energy_flux's last returned entry, Pi_{N-2}.
    pi = sabra_energy_flux(u, k_n)
    last_shell_de_dt = de_dt[-1]
    flux_anchored_boundary_check = -pi[-1] + last_shell_de_dt
    assert abs(flux_anchored_boundary_check - total_de_dt) < 1e-9 * max(
        1.0, abs(total_de_dt)
    ), (
        f"N={N}: -Pi_{{N-2}}+T_{{N-1}}={flux_anchored_boundary_check} "
        f"should equal sum_n T_n={total_de_dt}"
    )
    assert pi.shape[0] == N - 1, f"N={N}: expected {N - 1} inner boundaries, got {pi.shape[0]}"


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


@pytest.mark.parametrize("N", [5, 10, 22])
def test_h_invariant_catches_energy_conserving_wrong_coefficients(N):
    """RED CONTROL (requested by independent review, 2026-10-03): a
    coefficient choice that still satisfies a+b+c=0 (so a bare dE/dt=0
    energy test PASSES it) but is NOT the physical model must be caught by
    a SECOND invariant. For the standard model (a,b,c)=(1,-1/2,-1/2), the
    primary source's second invariant is
        H = sum_n (a/c)^n |u_n|^2 = sum_n (-2)^n |u_n|^2,
    conserved (dH/dt=0, no dissipation/forcing) ONLY for a+b+c=0 AND the
    specific ratio a/c=-2 pinned by (a,b,c)=(1,-1/2,-1/2) -- not for every
    energy-conserving triple. (1,-0.3,-0.7) satisfies a+b+c=0 but has
    a/c = 1/-0.7 =/= -2, so H must drift for it."""
    k_n = _k_of(N)
    rng = np.random.default_rng(42 + N)
    u = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    weights = (-2.0) ** np.arange(N)

    nl_correct = _paper_literal_nonlinear(u, k_n, a=1.0, b=-0.5, c=-0.5)
    dH_correct = float(np.sum(weights * 2.0 * np.real(np.conj(u) * nl_correct)))
    assert abs(dH_correct) < 1e-9 * max(1.0, np.sum(np.abs(u) ** 2) * 2.0 ** N), (
        f"N={N}: dH/dt={dH_correct} for the correct model, expected ~0"
    )

    nl_wrong = _energy_conserving_wrong_coeff_nonlinear(u, k_n)
    # The wrong-coefficient model must still conserve ENERGY (a+b+c=0)...
    dE_wrong = float(np.sum(2.0 * np.real(np.conj(u) * nl_wrong)))
    assert abs(dE_wrong) < 1e-9 * max(1.0, np.sum(np.abs(u) ** 2)), (
        f"N={N}: the negative control (1,-0.3,-0.7) should still conserve "
        f"energy (a+b+c=0) -- got dE/dt={dE_wrong}"
    )
    # ...but H must generically drift for it (it is not pinned to a/c=-2).
    dH_wrong = float(np.sum(weights * 2.0 * np.real(np.conj(u) * nl_wrong)))
    assert abs(dH_wrong) > 1e-6, (
        f"N={N}: expected the H-invariant to catch the energy-conserving "
        f"wrong coefficients (1,-0.3,-0.7) -- got dH/dt={dH_wrong} (too "
        f"close to 0, red control no longer discriminates)"
    )


@pytest.mark.parametrize("N", [5, 10, 22])
@pytest.mark.parametrize(
    "fixture_name", ["goy_sabra_nonlinear", "ledger_sabra_nonlinear"]
)
def test_h_invariant_conserved_by_production_functions(request, fixture_name, N):
    """Endabnahme requirement (2026-10-03, A3): anchor the H-invariant
    check to BOTH actual shipped production functions, not only to the
    helper ``_paper_literal_nonlinear`` used by
    ``test_h_invariant_catches_energy_conserving_wrong_coefficients``
    above. H = sum_n (-2)^n |u_n|^2 must be conserved (dH/dt=0) by the
    real ``sabra_nonlinear`` of both compute_goy_shell_dfc.py and
    compute_shell_dfc_waterline_ledger.py for the standard model
    (a,b,c)=(1,-1/2,-1/2). This is the regression anchor: a coefficient
    mutation in either production function (e.g. -0.5->-0.3,
    +0.5->+0.7 -- the same energy-conserving but H-violating triple as
    the negative control above) must turn this test red for the mutated
    fixture while leaving the other one green."""
    nonlinear = request.getfixturevalue(fixture_name)
    k_n = _k_of(N)
    rng = np.random.default_rng(42 + N)
    u = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    weights = (-2.0) ** np.arange(N)

    nl = nonlinear(u, k_n)
    dH = float(np.sum(weights * 2.0 * np.real(np.conj(u) * nl)))
    assert abs(dH) < 1e-9 * max(1.0, np.sum(np.abs(u) ** 2) * 2.0 ** N), (
        f"{fixture_name}, N={N}: dH/dt={dH} for the production "
        f"nonlinearity, expected ~0"
    )


@pytest.mark.parametrize("N", [3, 5, 8, 15])
@pytest.mark.parametrize(
    "fixture_name", ["goy_sabra_nonlinear", "ledger_sabra_nonlinear"]
)
def test_literal_eq21_transcription_matches_code_under_k_over_lambda(
    request, fixture_name, N
):
    """Direct check (requested by independent review, 2026-10-03) of the
    code's sabra_nonlinear against a LITERAL, line-by-line transcription
    of Eq. (21) of the primary source, using the paper's own
    (k_{n+1}, k_n, k_{n-1}) index convention. Per sabra_flux.py's module
    docstring ("Index convention vs. the primary source"), the code
    convention (k_n, k_{n-1}, k_{n-2}) is equivalent to the paper form
    evaluated at k -> k/LAM (LAM=2 here, geometric k_n=K0*LAM^n).
    Endabnahme requirement (2026-10-03, A3): run against BOTH actual
    production functions (``fixture_name``), not only the GOY script."""
    code_nonlinear = request.getfixturevalue(fixture_name)
    lam = 2.0
    k_n = (2.0 ** -4) * lam ** np.arange(N, dtype=float)
    rng = np.random.default_rng(1234 + N)
    u = rng.standard_normal(N) + 1j * rng.standard_normal(N)

    code_nl = code_nonlinear(u, k_n)
    paper_nl_over_lambda = _paper_literal_nonlinear(u, k_n / lam, a=1.0, b=-0.5, c=-0.5)

    assert np.allclose(code_nl, paper_nl_over_lambda, atol=1e-9, rtol=1e-7), (
        f"{fixture_name}, N={N}: code sabra_nonlinear does not match the "
        f"literal Eq. (21) transcription under k -> k/LAM "
        f"(max diff {np.max(np.abs(code_nl - paper_nl_over_lambda)):.3e})"
    )
