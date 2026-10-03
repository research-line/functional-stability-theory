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

Index/boundary correction (2026-10-03, independent review follow-up)
----------------------------------------------------------------------
An independent mathematical review (``_proof-notes/SABRA_FABLE_REVIEW_2026-10-03.md``)
found that the first version of this module, while correct on its interior
boundaries, still dropped the LAST valid inner boundary: it returned only
N-2 values (Pi_0..Pi_{N-3}) and both calling scripts zero-padded the
missing Pi_{N-2} = -k_{N-3}*W_{N-3} entry, which is in general nonzero (it
is the flux INTO the last shell, not a system boundary). This module now
returns the full N-1 values Pi_0..Pi_{N-2}, i.e. the flux across every one
of the N-1 *inner* boundaries of the shell chain; see "Boundary conditions"
below for exactly which two fluxes (not N-2!) are zero by construction.

Index convention vs. the primary source (explicit, 2026-10-03)
------------------------------------------------------------------
``sabra_nonlinear`` (and this module) index the triad with (k_n, k_{n-1},
k_{n-2}) instead of the primary source's (k_{n+1}, k_n, k_{n-1})
[L'vov et al. 1998, Eq. (21)]. For geometric k_n = K0*LAM^n this is *not*
merely a relabeling: it is a global 1/LAM rescaling of the nonlinear term,
equivalent to running the LITERAL paper dynamics in rescaled time
t' = t/LAM with an EFFECTIVE viscosity nu_eff = LAM*nu and an EFFECTIVE
forcing f_eff = LAM*f (at fixed nu, f this is not a mere time relabeling,
since dissipation and forcing do not rescale the same way as the
quadratic nonlinear term). The inviscid energy-conservation condition
a+b+c=0 is unaffected by this rescaling (it is a statement about the
nonlinear term alone), but any nu/f value quoted elsewhere (e.g. a
docstring saying "nu=1e-7") is this *code*-convention nu, which equals
1/LAM times the nu one would plug into the literal Eq. (21) form (i.e.
nu_code = nu_eff / LAM, consistent with nu_eff = LAM*nu_code above). This
matches the independent review's finding (section 1.5 of the review
note above); it does not change any already-verified energy-conservation
result in this module.

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
random and edge-case states, to machine/solver precision. The closed form,
for EVERY inner boundary n = 0..N-2 (i.e. all N-1 of them, not just N-2):

    W_m      := Im(u_m * u_{m+1} * conj(u_{m+2}))           for m = 0..N-3
    Pi_0      = -2*k_0*W_0
    Pi_n      = -2*k_n*W_n - k_{n-1}*W_{n-1}                 for n = 1..N-3
    Pi_{N-2}  = -k_{N-3}*W_{N-3}                             (only the
                "incoming" triad (N-3,N-2,N-1) contributes; the triad
                (N-2,N-1,N) does not exist, u_N := 0)

i.e. the flux across the boundary after shell n is carried by (up to) the
two triads that straddle that boundary: triad (n-1, n, n+1) (weight
k_{n-1}, since only its "outside" leg n+1 crosses) and triad
(n, n+1, n+2) (weight 2*k_n, since two of its three legs are outside) --
except at the very last inner boundary (n = N-2), where the second triad
does not exist and only the first contributes. This two-triad structure
is the standard shape of shell-model energy flux (consistent with the
energy-conserving coefficient choice a+b-c=0 used by ``sabra_nonlinear``);
it collapses to a single boundary term only in special/approximate
conventions (or, structurally, at the last inner boundary), which is why
the single-term OLD formula above is not a valid flux for this (correct)
dynamics at any interior boundary.

By construction (global energy conservation of the corrected nonlinearity,
verified in test_sabra_energy_conservation.py), Pi_n so defined also obeys
the telescoping identity

    dE_n/dt|_NL = Pi_{n-1} - Pi_n      for n = 0, ..., N-1

with the two SYSTEM-boundary conventions Pi_{-1} := 0 (nothing flows into
shell 0 from a nonexistent "shell -1") and Pi_{N-1} := 0 (global energy
conservation of the nonlinear term: nothing is left to flow out past the
last shell). Pi_{N-1} is a true structural zero and is intentionally NOT
part of this function's return value (there is no "N-th inner boundary");
Pi_{N-2} IS returned and is, in general, nonzero -- it is the flux INTO
the last shell, dE_{N-1}/dt|_NL = Pi_{N-2}, not a system boundary.
"""

from __future__ import annotations

import numpy as np


def sabra_energy_flux(u: np.ndarray, k_n: np.ndarray) -> np.ndarray:
    """Energy-consistent cumulative shell flux Pi_n for the Sabra model.

    Index convention: the SAME (k_n, k_{n-1}, k_{n-2}) convention as
    ``sabra_nonlinear`` (not the primary source's (k_{n+1}, k_n, k_{n-1}));
    see the module docstring, section "Index convention vs. the primary
    source", for why this is not a mere relabeling at fixed viscosity and
    forcing.

    Parameters
    ----------
    u : complex ndarray, shape (N,)
        Shell amplitudes.
    k_n : real ndarray, shape (>= N,)
        Shell wavenumbers (only k_n[:N-2] is used).

    Returns
    -------
    ndarray, shape (N-1,), real
        Pi[n] for n = 0, ..., N-2: the net nonlinear energy flow out of the
        cumulative block of shells {0, ..., n}, for EVERY inner boundary of
        the shell chain (positive = forward/downscale cascade), consistent
        with ``dE_n/dt|_NL = Pi_{n-1} - Pi_n`` for n = 0, ..., N-1, using
        the two system-boundary conventions Pi_{-1} := 0 and Pi_{N-1} := 0
        (the latter is NOT part of the returned array -- it is the true
        structural zero "past the last shell"; Pi_{N-2}, the last entry
        this function DOES return, is in general nonzero).

    Raises
    ------
    ValueError
        If N < 3 (no triad exists, flux is undefined).
    """
    N = u.shape[0]
    if N < 3:
        raise ValueError(f"sabra_energy_flux needs at least 3 shells, got N={N}")

    w = np.imag(u[: N - 2] * u[1 : N - 1] * np.conj(u[2:N]))  # W_m, m=0..N-3

    pi = np.zeros(N - 1)
    pi[: N - 2] += -2.0 * k_n[: N - 2] * w  # Pi_n "own" triad, n = 0..N-3
    pi[1 : N - 1] += -k_n[: N - 2] * w      # Pi_n "incoming" triad, n = 1..N-2
    return pi
