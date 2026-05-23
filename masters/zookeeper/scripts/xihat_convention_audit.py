"""
xihat_convention_audit.py
=========================
Klaert die ξ̂-Konvention: welche Transformierte ist das Kriteriumsobjekt?

Connes-vS Thm 1.2 / K0_ERR_REDUCTION Z.36: der Operator A lebt auf
L^2([-L/2, L/2]) -- ZENTRIERT. Die Bibliothek rechnet auf [0,L].

Zwei Transformierte, Phasenfaktor e^{izL/2} auseinander:
  F(z)            = Sum_j xi_j / (z - 2 pi j / L)        (= CCM Eq 5.25, rational)
  xihat_0L(z)     = (1 - e^{-izL})/(i sqrt(L)) * F(z)     ([0,L]-Transform)
  xihat_zentr(z)  = (2/sqrt(L)) sin(zL/2) * F(z)          (zentriert; = e^{izL/2} xihat_0L)

Behauptungen, die hier numerisch geprueft werden:
  (A) Eigenvektor-Koeffizienten sind even: xi[N+j] = xi[N-j]  -> F ist ODD.
  (B) xihat_zentr ist GERADE: xihat_zentr(-z) = xihat_zentr(z).
  (C) xihat_zentr ist reell auf der reellen Achse.
  (D) xihat_zentr(z) = e^{izL/2} * xihat_0L(z)   (Phasenrelation).
  (E) xihat_zentr an den Gitterpunkten 2 pi k / L ist ENDLICH (Pole getilgt).
  (F) xihat_0L erfuellt KEINE einfache Geradheit, sondern
      xihat_0L(z) = e^{-izL} xihat_0L(-z).
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
    p.add_argument("--dps", type=int, default=50)
    args = p.parse_args()
    mp.mp.dps = args.dps

    M, L_mp = build_QW_mp(args.N, args.lam, chi_func=chi_trivial, q_mod=1,
                          include_W02=True, sign_WR=-1, parity="even",
                          conductor_correction=False, verbose=False)
    size, dim = 2 * args.N + 1, args.N + 1
    M_even, U_even = project_to_parity(M, args.N, parity="even")
    w, Q = diagonalize_mp(M_even, verbose=False)
    xi = embed(Q, U_even, 0, size, dim)
    L = L_mp
    N = args.N

    def omega(idx):
        return 2 * mp.pi * (idx - N) / L

    def F(z):
        s = mp.mpf(0)
        for idx in range(size):
            s += xi[idx, 0] / (z - omega(idx))
        return s

    def xihat_0L(z):
        return (1 - mp.e ** (-1j * z * L)) / (1j * mp.sqrt(L)) * F(z)

    def xihat_zentr(z):
        return (2 / mp.sqrt(L)) * mp.sin(z * L / 2) * F(z)

    print(f"lambda={args.lam}  N={N}  L={mp.nstr(L,8)}  dps={args.dps}")

    # (A) Koeffizienten-Symmetrie xi[N+j] = xi[N-j]
    asym = mp.mpf(0)
    for j in range(1, N + 1):
        asym = max(asym, abs(xi[N + j, 0] - xi[N - j, 0]))
    print(f"(A) max|xi[N+j]-xi[N-j]| = {mp.nstr(asym,4)}   "
          f"-> Koeffizienten {'EVEN (F odd)' if asym < mp.mpf(10)**(-args.dps+10) else 'NICHT even'}")

    # Testpunkte abseits der Gitterpunkte
    zs = [mp.mpf("1.3"), mp.mpf("4.7"), mp.mpf("2.1") + mp.mpf("0.3") * 1j,
          mp.mpf("0.9") - mp.mpf("0.25") * 1j]

    # (B) xihat_zentr gerade
    devB = max(abs(xihat_zentr(-z) - xihat_zentr(z)) for z in zs)
    print(f"(B) max|xihat_zentr(-z) - xihat_zentr(z)| = {mp.nstr(devB,4)}")

    # (C) reell auf R
    devC = max(abs(mp.im(xihat_zentr(x))) for x in [mp.mpf("1.3"), mp.mpf("4.7"), mp.mpf("7.2")])
    print(f"(C) max|Im xihat_zentr(x)| auf R = {mp.nstr(devC,4)}")

    # (D) Phasenrelation
    devD = max(abs(xihat_zentr(z) - mp.e ** (1j * z * L / 2) * xihat_0L(z)) for z in zs)
    print(f"(D) max|xihat_zentr - e^{{izL/2}} xihat_0L| = {mp.nstr(devD,4)}")

    # (E) Endlichkeit an Gitterpunkten: xihat_zentr(2 pi k / L) vs sqrt(L)(-1)^k xi_k
    print("(E) xihat_zentr an Gitterpunkten 2 pi k / L:")
    for k in [1, 2, 3]:
        zk = 2 * mp.pi * k / L
        val = xihat_zentr(zk + mp.mpf(10) ** (-args.dps + 15))  # knapp daneben
        pred = mp.sqrt(L) * ((-1) ** k) * xi[N + k, 0]
        print(f"    k={k}: xihat~={mp.nstr(val,6)}  pred sqrt(L)(-1)^k xi_k={mp.nstr(pred,6)}"
              f"  diff={mp.nstr(abs(val-pred),3)}")

    # (F) xihat_0L Funktionalgleichung
    devF = max(abs(xihat_0L(z) - mp.e ** (-1j * z * L) * xihat_0L(-z)) for z in zs)
    print(f"(F) max|xihat_0L(z) - e^{{-izL}} xihat_0L(-z)| = {mp.nstr(devF,4)}")
    # und Gegencheck: ist xihat_0L NICHT gerade?
    devF2 = max(abs(xihat_0L(-z) - xihat_0L(z)) for z in zs)
    print(f"    Gegencheck max|xihat_0L(-z)-xihat_0L(z)| = {mp.nstr(devF2,4)} (sollte gross sein)")


if __name__ == "__main__":
    main()
