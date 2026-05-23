"""
eigenfunction_decay_audit.py
============================
Testet die Reduktion (ii-a) <= uniforme Lokalisierung der Eigenfunktion.

|F_lambda(x+iy)| <= int |xi_lambda(t)| e^{|y||t|} dt.

Wenn xi_lambda(t) in einem lambda-uniformen Kern nahe t=0 konzentriert ist
und davor abklingt, ist dieses Integral O(1) -> (ii-a) folgt.

Messgroessen pro lambda:
  - r_q  : Massenradius (Radius, der Anteil q der L^2-Masse enthaelt)
  - |xi(0)| vs |xi(L/2)| : Abklingen zum Rand
  - B(delta) = int |xi| e^{delta|t|} dt / ||xi||_2  (delta=0.4)
  - W(delta) = ( int e^{2 delta|t|} dt )^{1/2} / sqrt(L)  (Worst-Case-Vergleich)
Frage: bleibt r_q beschraenkt und B(delta) = O(1), waehrend W(delta) ~ lambda^delta waechst?
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import mpmath as mp

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / ".LAB" / ".ZETA-ZOO" / "CORE" / "zookeeper" / "_scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))  # Mac ~/compute: lokale Kopie
from dirichlet_ccm_fourier_mp import (  # noqa: E402
    build_QW_mp, chi_trivial, project_to_parity, diagonalize_mp,
)


def embed(Q, U_even, k, size, dim):
    xi = mp.matrix(size, 1)
    for row in range(size):
        s = mp.mpf(0)
        for r in range(dim):
            s += U_even[row, r] * Q[r, k]
        xi[row, 0] = s
    return xi


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--lambda", dest="lam", type=float, required=True)
    p.add_argument("--N", type=int, required=True)
    p.add_argument("--dps", type=int, default=60)
    p.add_argument("--delta", type=float, default=0.4)
    p.add_argument("--ngrid", type=int, default=600)
    args = p.parse_args()
    mp.mp.dps = args.dps

    M, L_mp = build_QW_mp(args.N, args.lam, chi_func=chi_trivial, q_mod=1,
                          include_W02=True, sign_WR=-1, parity="even",
                          conductor_correction=False, verbose=False)
    size, dim = 2 * args.N + 1, args.N + 1
    M_even, U_even = project_to_parity(M, args.N, parity="even")
    w, Q = diagonalize_mp(M_even, verbose=False)
    xi = embed(Q, U_even, 0, size, dim)   # Minimal-Eigenvektor, Koeffizienten

    L = L_mp
    half = L / 2
    delta = mp.mpf(args.delta)

    # Positionsraum-Rekonstruktion: xi_lambda(t) = sum_j xi_j e^{2 pi i j t / L}
    # (even -> reell, gerade). Gitter ueber [-L/2, L/2].
    ng = args.ngrid
    ts = [(-half + L * i / ng) for i in range(ng + 1)]

    def xi_of_t(t):
        # even-Eigenvektor -> reelle gerade Funktion: sum_j xi_j cos(2 pi j t / L)
        s = mp.mpf(0)
        for idx in range(size):
            j = idx - args.N
            s += xi[idx, 0] * mp.cos(2 * mp.pi * j * t / L)
        return s

    vals = [abs(xi_of_t(t)) for t in ts]

    # Trapez-Integrale
    def trap(f):
        h = L / ng
        s = (f[0] + f[-1]) / 2 + sum(f[1:-1])
        return s * h

    norm2 = mp.sqrt(trap([v * v for v in vals]))                    # ||xi||_2
    mass = [v * v for v in vals]
    total = trap(mass)
    # Massenradius r_q: kleinstes r mit int_{-r}^{r}|xi|^2 >= q*total
    def mass_radius(q):
        target = q * total
        h = L / ng
        # symmetrisch von der Mitte (Index ng/2) nach aussen
        mid = ng // 2
        acc = mass[mid] * h
        i = 1
        while acc < target and mid + i <= ng:
            acc += (mass[mid + i] + mass[mid - i]) * h
            i += 1
        return ts[mid + i] if mid + i <= ng else half

    r50 = mass_radius(mp.mpf("0.5"))
    r90 = mass_radius(mp.mpf("0.9"))
    r99 = mass_radius(mp.mpf("0.99"))

    B = trap([vals[i] * mp.e ** (delta * abs(ts[i])) for i in range(ng + 1)])
    B_ratio = B / norm2
    W = mp.sqrt(trap([mp.e ** (2 * delta * abs(t)) for t in ts])) / mp.sqrt(L)

    xi0 = abs(xi_of_t(mp.mpf(0)))
    xi_edge = abs(xi_of_t(half * mp.mpf("0.98")))

    print(f"lambda={args.lam}  N={args.N}  L={mp.nstr(L,6)}  dps={args.dps}")
    print(f"  Massenradius r50={mp.nstr(r50,5)}  r90={mp.nstr(r90,5)}  "
          f"r99={mp.nstr(r99,5)}   (Halbintervall L/2={mp.nstr(half,5)})")
    print(f"  |xi(0)|={mp.nstr(xi0,5)}  |xi(~Rand)|={mp.nstr(xi_edge,4)}  "
          f"Abkling-Verhaeltnis={mp.nstr(xi_edge/xi0,4)}")
    print(f"  B(delta={args.delta})/||xi||_2 = {mp.nstr(B_ratio,6)}   "
          f"Worst-Case W = {mp.nstr(W,6)}   (lambda^delta="
          f"{mp.nstr(mp.mpf(args.lam)**delta,5)})")


if __name__ == "__main__":
    main()
