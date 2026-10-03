"""
sabra_flux.py
==============
Shared, energy-consistent shell-to-shell energy flux Pi_n for the Sabra
shell model, matching the corrected ``sabra_nonlinear`` (GitHub Issue #1,
commit 48448f5, 2026-10-03).

Why this module exists
-----------------------
``compute_goy_shell_dfc.py`` and ``compute_shell_dfc_waterline_ledger.py``
each computed the energy flux locally as

    Pi_n = Im(k_n * u_n * conj(u_{n+1}) * u_{n+2})                 (OLD, wrong)

This is the flux formula that belonged to the *pre-fix* (buggy) nonlinear
term. It was never re-derived when the nonlinearity's conjugation/
coefficient bug was fixed, so it no longer matches the energy balance of
the corrected dynamics (independent finding, 2026-10-03, follow-up to
Issue #1).

Derivation (cumulative energy balance, not "from memory")
-----------------------------------------------------------
For the quadratic invariant E_n = |u_n|^2 and the corrected nonlinearity

    du_n/dt|_NL = i*(a*k_n*conj(u_{n+1})*u_{n+2}
                      + b*k_{n-1}*u_{n+1}*conj(u_{n-1})
                      + c*k_{n-2}*u_{n-1}*u_{n-2})
    with (a, b, c) = (1, -1/2, +1/2)   [L'vov, Podivilov, Pomyalov, Procaccia,
    Vandembroucq, "Improved shell model of turbulence", Phys. Rev. E 58, 1811
    (1998), Eq. (21), already verified and cited in commit 48448f5]

the nonlinear part of dE_n/dt is 2*Re(conj(u_n) * du_n/dt|_NL), and the
energy flowing OUT of the cumulative block of shells {0, ..., n} per unit
time (no dissipation, no forcing) is, by definition,

    Pi_n := -sum_{m=0}^{n} dE_m/dt|_NL = -2 * sum_{m=0}^{n} Re(conj(u_m) * N_m(u))

This module's ``sabra_energy_flux`` is a *closed-form* (O(N), vectorized)
expression for exactly this cumulative sum, obtained by symbolic
telescoping (sympy, see
``_proof-notes/SABRA_ENERGIE_ISSUE1_2026-10-03.md`` NACHTRAG 2026-10-03 for
the full derivation and verification script) and cross-checked numerically
against the brute-force cumulative-sum definition above for N=3..30,
random and edge-case states, to machine/solver precision. The closed form:

    W_m      := Im(u_m * u_{m+1} * conj(u_{m+2}))           for m = 0..N-3
    Pi_0      = -2*k_0*W_0
    Pi_n      = -2*k_n*W_n - k_{n-1}*W_{n-1}                 for n = 1..N-3

i.e. the flux across the boundary after shell n is carried by exactly the
two triads that straddle that boundary: triad (n-1, n, n+1) (weight
k_{n-1}, since only its "outside" leg n+1 crosses) and triad
(n, n+1, n+2) (weight 2*k_n, since two of its three legs are outside).
This two-triad combination is the standard structure of shell-model energy
flux (consistent with the energy-conserving coefficient choice a+b-c=0
used by ``sabra_nonlinear``); it collapses to a single boundary term only
in special/approximate conventions, which is why the single-term OLD
formula above is not a valid flux for this (correct) dynamics.

By construction (global energy conservation of the corrected nonlinearity,
verified in test_sabra_energy_conservation.py), Pi_n so defined also obeys
the telescoping identity

    dE_n/dt|_NL = Pi_{n-1} - Pi_n      for n = 1, ..., N-2
    dE_0/dt|_NL = -Pi_0                (no flux INTO shell 0 from "shell -1")
    dE_{N-1}/dt|_NL = Pi_{N-3}         (no flux OUT of the last shell; there
                                        is no triad beyond index N-3)

which also answers the "flux before the first shell is zero" / "no flux
leaves the last shell" boundary conditions: Pi is only ever defined for a
genuine triad (n = 0..N-3); there is nothing before shell 0 or after the
last triad by construction, not by a special-cased zero.
"""

from __future__ import annotations

import numpy as np


def sabra_energy_flux(u: np.ndarray, k_n: np.ndarray) -> np.ndarray:
    """Energy-consistent cumulative shell flux Pi_n for the Sabra model.

    Parameters
    ----------
    u : complex ndarray, shape (N,)
        Shell amplitudes.
    k_n : real ndarray, shape (>= N,)
        Shell wavenumbers (only k_n[:N-2] is used).

    Returns
    -------
    ndarray, shape (N-2,), real
        Pi[n] for n = 0, ..., N-3: the net nonlinear energy flow out of the
        cumulative block of shells {0, ..., n} (positive = forward/
        downscale cascade), consistent with
        ``dE_n/dt|_NL = Pi_{n-1} - Pi_n`` (Pi_{-1} := 0).

    Raises
    ------
    ValueError
        If N < 3 (no triad exists, flux is undefined).
    """
    N = u.shape[0]
    if N < 3:
        raise ValueError(f"sabra_energy_flux needs at least 3 shells, got N={N}")

    w = np.imag(u[: N - 2] * u[1 : N - 1] * np.conj(u[2:N]))  # W_m, m=0..N-3
    pi = -2.0 * k_n[: N - 2] * w
    if N > 3:
        pi[1:] += -k_n[: N - 3] * w[:-1]
    return pi
