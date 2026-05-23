"""
xihat_l2strip_audit.py
======================
Terminal-Ziel (ii-a), Versuch 2: fensterfreie L^2-Streifen-Observable.

Statt sup ueber ein festes Fenster (N-/fenster-instabil, siehe §E.17) die
GLOBALE, Plancherel-verankerte Groesse

    B_λ(δ) = ‖ξ̂_λ(· + iδ)‖_{L²(ℝ)} / ‖ξ̂_λ‖_{L²(ℝ)}.

Konturverschiebung + Plancherel:  ξ̂_λ(x+iδ) = FT von ξ_λ(t)e^{δt}, also

    B_λ(δ)^2 = ∫_{-L/2}^{L/2} |ξ_λ(t)|^2 e^{2δt} dt  /  ∫ |ξ_λ(t)|^2 dt
             = ∫ |ξ_λ(t)|^2 cosh(2δt) dt / ∫ |ξ_λ(t)|^2 dt   (ξ_λ gerade).

ξ_λ(t) = Σ_j ξ_j e^{2πijt/L}  (zentriert, ξ_j = ξ_{-j} reell). Mit der
Autokorrelation d_m = Σ_j ξ_j ξ_{m-j}  (|ξ_λ|^2 = Σ_m d_m e^{2πimt/L}) und

    ∫_{-L/2}^{L/2} cos(2πmt/L) cosh(2δt) dt = 4δ(-1)^m sinh(δL)/(4δ^2+(2πm/L)^2)

ist B_λ(δ)^2 = [ Σ_m d_m · 4δ(-1)^m sinh(δL)/(4δ^2+(2πm/L)^2) ] / (d_0 · L).
EXAKT, fensterfrei, nur ueber das kompakte Traegerintervall. Einzige Kosten:
die Diagonalisierung fuer ξ_λ. -> N-Konvergenz pro λ separat pruefbar.

Aufruf: python xihat_l2strip_audit.py --lambda L --N N [--dps 40]
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import mpmath as mp


def find_library():
    """dirichlet_ccm_fourier_mp.py finden: neben dem Skript (Mac ~/compute)
    oder ueber den Repo-Baum (Laptop/OneDrive)."""
    here = Path(__file__).resolve().parent
    cand = [here]
    p = here
    for _ in range(6):
        p = p.parent
        cand.append(p / ".LAB" / ".ZETA-ZOO" / "CORE" / "zookeeper" / "_scripts")
        cand.append(p / ".ZETA-ZOO" / "CORE" / "zookeeper" / "_scripts")
    for c in cand:
        if (c / "dirichlet_ccm_fourier_mp.py").exists():
            return str(c)
    raise FileNotFoundError("dirichlet_ccm_fourier_mp.py nicht gefunden")


sys.path.insert(0, find_library())
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
    p.add_argument("--dps", type=int, default=40)
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

    # Koeffizienten a_j, j = -N..N  (Index idx = j + N)
    a = [xi[idx, 0] for idx in range(size)]   # a[idx] = ξ_{idx-N}

    # Autokorrelation d_m = Σ_j a_j a_{m-j},  m = 0..2N  (d_m = d_{-m})
    d = [mp.mpf(0)] * (2 * N + 1)
    for m in range(2 * N + 1):
        s = mp.mpf(0)
        # j-Index (in idx-Raum): idx von max(0, m-2N).. ; a_j a_{m-j}
        # j = idx - N ; m - j = m - idx + N  -> zweiter idx = m - idx + 2N
        for idx in range(size):
            idx2 = m - (idx - N) + N
            if 0 <= idx2 < size:
                s += a[idx] * a[idx2]
        d[m] = s

    d0 = d[0]

    print(f"lambda={args.lam}  N={N}  L={mp.nstr(L,8)}  dps={args.dps}")
    print(f"  d_0 = Σ|ξ_j|^2 = {mp.nstr(d0,6)}")
    for dd in ["0.2", "0.3", "0.4"]:
        delta = mp.mpf(dd)
        sinhdL = mp.sinh(delta * L)
        num = d0 * (sinhdL / delta)            # m=0 Term: 4δ·sinh(δL)/(4δ^2)
        for m in range(1, 2 * N + 1):
            om = 2 * mp.pi * m / L
            wm = 4 * delta * ((-1) ** m) * sinhdL / (4 * delta * delta + om * om)
            num += 2 * d[m] * wm               # ±m
        B2 = num / (d0 * L)
        B = mp.sqrt(B2)
        lam_d = mp.mpf(args.lam) ** delta
        print(f"  δ={dd}:  B_λ(δ) = {mp.nstr(B,8)}   (λ^δ = {mp.nstr(lam_d,6)}  "
              f"Anteil = {mp.nstr(B/lam_d,5)})")


if __name__ == "__main__":
    main()
