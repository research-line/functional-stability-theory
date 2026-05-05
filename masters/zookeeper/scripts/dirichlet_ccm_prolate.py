"""
dirichlet_ccm_prolate.py
========================

Neu-Implementation des CCM-Operators in der Prolate-Basis (orthogonale
Polynome P_n fuer das Mass dm(s) = (2 pi)^(-3/2) |Gamma(1/4 + is/2)|^2 ds),
nach CCM 2024 (arXiv:2310.18423).

Kernkonstruktion (in P_n-Basis):

    S_{n, n+1}   = i * sqrt((n + 1/2)(n + 1))          (tridiagonal, Hermitisch)
    S_{n, n}     = 0
    W_lambda     = -S^2 + 2 pi lambda^2 (4 N_grad + 1) - 1/4

wobei N_grad = diag(0, 1, 2, ..., N-1).

Fuer den vollen Operator kommt eine arithmetische Korrektur (Prim-Summe) hinzu,
die ueber den Kern q(U_n, U_m)(log k) die Primzahlen einbaut. Die einfache
Version ohne Primen sollte bereits eine Naeherung an die Riemann-Nullstellen
liefern.

Autor: LG (Opus 4.7 [1M], Session 18 beta)
Datum: 2026-04-18
Ref:   Connes-Consani-Moscovici, arXiv:2310.18423 (CCM 2024)
"""

import os
import sys
import numpy as np
from pathlib import Path

# mpmath fuer Arbitrary Precision
import mpmath as mp

# =================================================================
# Konfiguration
# =================================================================

LAMBDA_DEFAULT = float(np.sqrt(14.0))
N_DEFAULT = 50
DPS_DEFAULT = 50

LAMBDA = float(os.environ.get("LAMBDA", LAMBDA_DEFAULT))
N = int(os.environ.get("N", N_DEFAULT))
DPS = int(os.environ.get("DPS", DPS_DEFAULT))

# Bekannte Riemann-Nullstellen (LMFDB, hochgenau)
RIEMANN_ZEROS = [
    14.134725141734693,
    21.022039638771554,
    25.010857580145688,
    30.424876125859513,
    32.935061587739189,
    37.586178158825671,
    40.918719012147495,
    43.327073280914999,
    48.005150881167159,
    49.773832477672302,
]

# =================================================================
# Matrix-Konstruktion
# =================================================================


def build_S_matrix(N, dps=50):
    """
    Baut Scaling-Operator S als (N x N) Hermitesche Tridiagonalmatrix:

        S[n, n+1] = i * sqrt((n + 1/2)(n + 1))
        S[n+1, n] = conj(S[n, n+1]) = -i * sqrt((n + 1/2)(n + 1))
        S[n, n]   = 0

    Rueckgabe: mpmath.matrix (komplex).
    """
    mp.mp.dps = dps
    S = mp.matrix(N, N)
    for n in range(N - 1):
        val = mp.mpc(0, 1) * mp.sqrt(mp.mpf(n) + mp.mpf("0.5")) * mp.sqrt(mp.mpf(n) + 1)
        S[n, n + 1] = val
        S[n + 1, n] = mp.conj(val)
    return S


def build_W_lambda(N, lam, dps=50):
    """
    Baut W_lambda = -S^2 + 2 pi lambda^2 (4 N_grad + 1) - 1/4

    in der P_n-Basis, als (N x N) Hermitesche mpmath-Matrix.
    """
    mp.mp.dps = dps
    S = build_S_matrix(N, dps=dps)

    # S^2 (mpmath-Matrixmultiplikation)
    S2 = S * S

    # W = -S^2 + 2 pi lambda^2 (4 N_grad + 1) - 1/4
    W = mp.matrix(N, N)
    two_pi_lam2 = mp.mpf(2) * mp.pi * mp.mpf(lam) ** 2
    quarter = mp.mpf("0.25")

    for i in range(N):
        for j in range(N):
            W[i, j] = -S2[i, j]
        # Diagonale: + 2 pi lambda^2 (4 i + 1) - 1/4
        W[i, i] += two_pi_lam2 * (mp.mpf(4) * mp.mpf(i) + 1) - quarter

    return W


# =================================================================
# Hermitisierung-Check und Diagonalisierung
# =================================================================


def mpmat_to_numpy(M):
    """Konvertiert mpmath-Matrix zu numpy (komplex, float128 oder complex)."""
    N = M.rows
    arr = np.zeros((N, N), dtype=np.complex128)
    for i in range(N):
        for j in range(N):
            arr[i, j] = complex(M[i, j])
    return arr


def verify_hermitian(W_mp):
    """Prueft, ob die Matrix Hermitesch ist (im mpmath-Sinn)."""
    N = W_mp.rows
    max_diff = mp.mpf(0)
    for i in range(N):
        for j in range(N):
            diff = abs(W_mp[i, j] - mp.conj(W_mp[j, i]))
            if diff > max_diff:
                max_diff = diff
    return float(max_diff)


def diagonalize_mpmath(W_mp, dps=50):
    """
    Diagonalisiert eine Hermitesche mpmath-Matrix mit mp.eigsy (nach Hermitisierung
    zu rein-realer Darstellung via Realteil- und Imaginaerteil-Block).

    Alternativ: numpy.linalg.eigh (Float64) fuer schnelle Voranalyse.

    Rueckgabe: (eigenvalues, eigenvectors) als numpy-Arrays.
    """
    # Erstmal numpy-Voranalyse (schnell)
    W_np = mpmat_to_numpy(W_mp)
    # Enforce numerical Hermiticity
    W_np = 0.5 * (W_np + W_np.conj().T)
    w, V = np.linalg.eigh(W_np)
    return w, V


# =================================================================
# Verifikations-Tests
# =================================================================


def test_S_hermiticity(N=20, dps=50):
    """S muss Hermitesch sein: S[n,m] = conj(S[m,n])."""
    S = build_S_matrix(N, dps=dps)
    max_diff = verify_hermitian(S)
    print(f"[test_S_hermiticity] N={N}, max |S - S*| = {max_diff:.3e}")
    assert max_diff < 1e-40, f"S nicht Hermitesch: {max_diff}"


def test_S2_diagonal(N=10, dps=50):
    """
    S^2[n,n] = (n-1/2)*n + (n+1/2)*(n+1) = 2n^2 + n + 1/2 (analytische Formel)
    S^2[n,n+2] = -sqrt((n+1/2)(n+1)(n+3/2)(n+2))
    """
    S = build_S_matrix(N, dps=dps)
    S2 = S * S

    print(f"\n[test_S2_diagonal] N={N}")
    for n in range(min(5, N)):
        expected = 2 * n * n + n + mp.mpf("0.5")
        actual = S2[n, n]
        diff = abs(actual - expected)
        print(
            f"  n={n}: S^2[n,n] = {float(mp.re(actual)):.6f}"
            f"+{float(mp.im(actual)):.2e}i   "
            f"expected={float(expected):.6f}   diff={float(diff):.3e}"
        )

    for n in range(min(4, N - 2)):
        expected = -mp.sqrt(
            (mp.mpf(n) + mp.mpf("0.5"))
            * (mp.mpf(n) + 1)
            * (mp.mpf(n) + mp.mpf("1.5"))
            * (mp.mpf(n) + 2)
        )
        actual = S2[n, n + 2]
        diff = abs(actual - expected)
        print(
            f"  n={n}: S^2[n,n+2] = {float(mp.re(actual)):.6f}"
            f"+{float(mp.im(actual)):.2e}i   "
            f"expected={float(expected):.6f}   diff={float(diff):.3e}"
        )


def test_W_hermiticity(N=20, lam=LAMBDA_DEFAULT, dps=50):
    W = build_W_lambda(N, lam, dps=dps)
    max_diff = verify_hermitian(W)
    print(f"\n[test_W_hermiticity] N={N}, lambda={lam:.6f}, max |W-W*| = {max_diff:.3e}")


# =================================================================
# Haupt-Diagnose
# =================================================================


def diagnose_spectrum(N=N, lam=LAMBDA, dps=DPS, n_low=20):
    """
    Baut W_lambda, diagonalisiert und zeigt die niedrigsten Eigenwerte.
    """
    print("=" * 70)
    print(f"CCM-Prolate-Basis Diagnose")
    print(f"  N      = {N}")
    print(f"  lambda = {lam:.6f}  (sqrt(14) = {float(np.sqrt(14)):.6f})")
    print(f"  dps    = {dps}")
    print("=" * 70)

    print("\n1. Baue S und W_lambda...")
    W = build_W_lambda(N, lam, dps=dps)
    herm_err = verify_hermitian(W)
    print(f"   max |W - W*| = {herm_err:.3e}")

    print("\n2. Diagonalisiere...")
    w, V = diagonalize_mpmath(W, dps=dps)

    print(f"\n3. Niedrigste {n_low} Eigenwerte:")
    print(f"   {'k':>3} {'eigenwert':>20} {'sqrt(|ew|)':>15}")
    for k in range(min(n_low, len(w))):
        ew = w[k]
        sq = np.sqrt(abs(ew)) if ew > 0 else -np.sqrt(abs(ew))
        print(f"   {k:>3} {ew:>20.10f} {sq:>15.6f}")

    print(f"\n4. Vergleich mit Riemann-Nullstellen gamma_k^2:")
    print(f"   {'k':>3} {'gamma_k':>12} {'gamma_k^2':>15} "
          f"{'naechster ew':>15} {'rel-error':>12}")
    for k in range(min(5, len(RIEMANN_ZEROS))):
        gamma = RIEMANN_ZEROS[k]
        gamma2 = gamma * gamma
        # Suche naechsten Eigenwert
        idx = int(np.argmin(np.abs(w - gamma2)))
        ew = w[idx]
        rel = abs(ew - gamma2) / gamma2
        print(
            f"   {k + 1:>3} {gamma:>12.6f} {gamma2:>15.6f} "
            f"{ew:>15.6f} {rel:>12.3e}"
        )

    return w, V


# =================================================================
# Entrypoint
# =================================================================


def batch_scan(N_list, lam=LAMBDA, dps=DPS, out_dir=None):
    """
    Laeuft diagnose_spectrum fuer eine Liste von N-Werten ohne Primen.
    Speichert Ergebnisse als JSON.
    """
    import json
    import time

    results = []
    for Nv in N_list:
        print(f"\n\n{'#' * 70}")
        print(f"# BATCH: N = {Nv}")
        print(f"{'#' * 70}")
        t0 = time.time()
        w, _ = diagnose_spectrum(N=Nv, lam=lam, dps=dps, n_low=15)
        elapsed = time.time() - t0

        # Fuer jede bekannte Riemann-Nullstelle: naechster EW + rel-Fehler
        gamma_matches = []
        for gamma in RIEMANN_ZEROS[:5]:
            gamma2 = gamma * gamma
            idx = int(np.argmin(np.abs(w - gamma2)))
            ew = float(w[idx])
            rel = abs(ew - gamma2) / gamma2
            gamma_matches.append(
                {"gamma": gamma, "gamma2": gamma2, "nearest_ew": ew, "rel_error": rel}
            )

        results.append(
            {
                "N": Nv,
                "lambda": lam,
                "dps": dps,
                "elapsed_sec": elapsed,
                "eigenvalues_low15": [float(x) for x in w[:15]],
                "gamma_matches": gamma_matches,
            }
        )

        if out_dir is not None:
            out_path = Path(out_dir) / f"prolate_scan_N{Nv}_lam{lam:.6f}.json"
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text(json.dumps(results[-1], indent=2))
            print(f"\n   Saved: {out_path}")

    return results


if __name__ == "__main__":
    mp.mp.dps = DPS

    mode = os.environ.get("MODE", "diagnose")

    if mode == "batch":
        # N-Liste aus ENV, Default: 50, 100, 200
        N_list_str = os.environ.get("N_LIST", "50,100,200")
        N_list = [int(x) for x in N_list_str.split(",")]
        out_dir = Path(__file__).parent.parent / "_results"
        batch_scan(N_list, lam=LAMBDA, dps=DPS, out_dir=out_dir)
    else:
        # Verifikationen
        test_S_hermiticity(N=20)
        test_S2_diagonal(N=10)
        test_W_hermiticity(N=20, lam=LAMBDA)

        # Haupt-Diagnose
        w, V = diagnose_spectrum(N=N, lam=LAMBDA, dps=DPS, n_low=15)
