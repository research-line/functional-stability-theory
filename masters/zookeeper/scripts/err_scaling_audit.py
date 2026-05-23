"""
err-Skalierung der CCM-Cluster-Route -- praezisionsfeste Fassung.

Korrektur gegenueber v1 und gegenueber ms1_splitting_audit.py:
  * Riemann-Nullstellen via mp.zetazero (volle Praezision), nicht 16-stellige
    Python-floats.
  * err (Nullstellen-LAGE) ist skalenfrei: F -> cF laesst Nullstellen fest.
    Daher wird der L^2-eingebettete Eigenvektor verwendet (keine
    delta_N-Aufblaehung), und F-Nullstellen werden in einem POLFREIEN
    Kleinfenster um jede Riemann-Nullstelle bisektiert.

Frage: geht err mit lambda -> infinity gegen 0, oder bodet es aus?
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import mpmath as mp

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / ".LAB" / ".ZETA-ZOO" / "CORE" / "zookeeper" / "_scripts"))
from dirichlet_ccm_fourier_mp import (  # noqa: E402
    build_QW_mp, chi_trivial, project_to_parity, diagonalize_mp,
)


def F_mp(xi, N, L_mp, z):
    s = mp.mpf(0)
    for idx in range(2 * N + 1):
        d = z - 2 * mp.pi * (idx - N) / L_mp
        if d == 0:
            return mp.inf
        s += xi[idx, 0] / d
    return s


def nearest_pole_dist(gamma, N, L_mp):
    """Abstand von gamma zur naechsten Polstelle 2 pi j / L."""
    best = None
    for j in range(-N, N + 1):
        d = abs(gamma - 2 * mp.pi * j / L_mp)
        if best is None or d < best:
            best = d
    return best


def embed(Q, U_even, k, size, dim):
    """k-ter even-Eigenvektor in den vollen Raum -- L^2, KEINE Sum-Skalierung."""
    xi = mp.matrix(size, 1)
    for row in range(size):
        s = mp.mpf(0)
        for r in range(dim):
            s += U_even[row, r] * Q[r, k]
        xi[row, 0] = s
    return xi


def err_at(xi, N, L_mp, gamma):
    """|F-Nullstelle - gamma| in einem polfreien Kleinfenster um gamma.
    Fensterbreite adaptiv: stets polfrei (< 0.45 * Abstand zur naechsten
    Polstelle), gedeckelt bei 0.08. Rueckgabe None ohne Vorzeichenwechsel."""
    half = min(mp.mpf("0.08"), mp.mpf("0.45") * nearest_pole_dist(gamma, N, L_mp))
    lo, hi = gamma - half, gamma + half
    n_scan = 96
    step = (hi - lo) / n_scan
    prev_z, prev_f = lo, F_mp(xi, N, L_mp, lo)
    best = None
    for i in range(1, n_scan + 1):
        z = lo + step * i
        f = F_mp(xi, N, L_mp, z)
        if prev_f * f < 0:
            a, b, fa = prev_z, z, prev_f
            for _ in range(4 * mp.mp.dps):
                m = (a + b) / 2
                fm = F_mp(xi, N, L_mp, m)
                if fa * fm <= 0:
                    b = m
                else:
                    a, fa = m, fm
            err = abs((a + b) / 2 - gamma)
            if best is None or err < best:
                best = err
        prev_z, prev_f = z, f
    return best


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--lambda", dest="lam", type=float, required=True)
    p.add_argument("--N", type=int, required=True)
    p.add_argument("--dps", type=int, default=110)
    p.add_argument("--threshold", type=float, default=1e-10)
    p.add_argument("--num-zeros", type=int, default=5,
                   help="Anzahl der Riemann-Nullstellen gamma_j")
    p.add_argument("--max-k", type=int, default=0,
                   help="nur die untersten max_k Cluster-Vektoren (0 = alle)")
    args = p.parse_args()
    mp.mp.dps = args.dps

    gammas = [mp.zetazero(n).imag for n in range(1, args.num_zeros + 1)]
    M, L_mp = build_QW_mp(args.N, args.lam, chi_func=chi_trivial, q_mod=1,
                          include_W02=True, sign_WR=-1, parity="even",
                          conductor_correction=False, verbose=False)
    size, dim = 2 * args.N + 1, args.N + 1
    M_even, U_even = project_to_parity(M, args.N, parity="even")
    w, Q = diagonalize_mp(M_even, verbose=False)

    w_min = w[0]
    cluster = [k for k in range(dim)
               if abs(float(w[k] - w_min)) <= args.threshold]
    print(f"lambda={args.lam}  N={args.N}  dps={args.dps}  "
          f"even_cluster={len(cluster)}")

    probe = cluster if args.max_k <= 0 else cluster[:args.max_k]
    best_overall, best_k = None, None
    for k in probe:
        xi = embed(Q, U_even, k, size, dim)
        errs = [err_at(xi, args.N, L_mp, g) for g in gammas]
        valid = [e for e in errs if e is not None]
        kmax = max(valid) if len(valid) == args.num_zeros else None
        shown = " ".join(mp.nstr(e, 3) if e is not None else ">8e-2"
                         for e in errs)
        tag = mp.nstr(kmax, 4) if kmax is not None else "(unvollst.)"
        print(f"  k={k:3d}  err(g1..{args.num_zeros})= {shown}   max_err={tag}")
        if kmax is not None and (best_overall is None or kmax < best_overall):
            best_overall, best_k = kmax, k
    print(f"  ==> best vector k={best_k}  max_err over g1..{args.num_zeros} = "
          f"{mp.nstr(best_overall,5) if best_overall else '----'}")


if __name__ == "__main__":
    main()
