"""
xihat_strip_audit.py
====================
Terminal-Ziel (ii-a), empirisch: waechst das kanonische Kriteriumsobjekt
ξ̂_λ auf einem FESTEN Streifen-Kompaktum, oder bleibt es O(1)?

ξ̂(z) = (2/√L) sin(zL/2) F(z),  F(z) = Σ_j ξ_j/(z - 2πj/L),  ξ_λ = k=0 Eigenvektor.

§E.7-Prognose (generisch): Streifen-Amplifikation A_λ(δ) ~ λ^δ.
Gemessen wird die NORMIERUNGSFREIE Groesse

    A_λ(δ) = sup_{|x|≤R} |ξ̂_λ(x + iδ)|  /  sup_{|x|≤R} |ξ̂_λ(x)|

(Verhaeltnis Streifenlinie / reelle Achse). Faellt A_λ(δ) ~ λ^δ aus, ist (ii-a)
genuin hart (Φ noetig, §E.7 bestaetigt). Bleibt A_λ(δ) = O(1), ist {ξ̂_λ}
nicht-generisch streifen-mild -> (ii-a) plausibel mit Φ = reelle Normierung.
Zusaetzlich: das rohe sup auf der reellen Achse (Wachstum vs √(logλ), §E.7 Beob 2).
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
    p.add_argument("--dps", type=int, default=30)
    p.add_argument("--R", type=float, default=10.0)
    p.add_argument("--h", type=float, default=0.02)
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
    omega = [2 * mp.pi * (idx - N) / L for idx in range(size)]

    def xihat(z):
        s = mp.mpf(0)
        for idx in range(size):
            s += xi[idx, 0] / (z - omega[idx])
        return (2 / mp.sqrt(L)) * mp.sin(z * L / 2) * s

    # Grid auf der reellen Achse + Streifenlinien; Start leicht irrational
    # versetzt, damit kein Gitterpunkt 2πj/L exakt getroffen wird.
    R, h = mp.mpf(args.R), mp.mpf(args.h)
    x0 = -R + mp.mpf("0.0013")
    npts = int((2 * R) / h)
    xs = [x0 + h * i for i in range(npts + 1)]

    def sup_on_line(y):
        m = mp.mpf(0)
        for x in xs:
            v = abs(xihat(mp.mpc(x, y)))
            if v > m:
                m = v
        return m

    sup_real = sup_on_line(mp.mpf(0))
    res = {}
    for d in [mp.mpf("0.2"), mp.mpf("0.4")]:
        # ξ̂ gerade & reell-auf-reell -> |ξ̂(x+iδ)| = |ξ̂(x-iδ)|; eine Linie genuegt
        s = sup_on_line(d)
        res[float(d)] = (s, s / sup_real)

    print(f"lambda={args.lam}  N={N}  L={mp.nstr(L,6)}  R={args.R}  h={args.h}  dps={args.dps}")
    print(f"  sup|ξ̂| auf [−R,R] (reelle Achse) = {mp.nstr(sup_real,6)}")
    for d, (s, ratio) in res.items():
        lam_d = mp.mpf(args.lam) ** d
        print(f"  δ={d}:  sup|ξ̂(x+iδ)| = {mp.nstr(s,6)}   "
              f"A_λ(δ) = {mp.nstr(ratio,6)}   (λ^δ = {mp.nstr(lam_d,5)})")


if __name__ == "__main__":
    main()
