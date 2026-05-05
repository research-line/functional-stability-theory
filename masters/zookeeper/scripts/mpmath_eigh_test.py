"""
mpmath_eigh_test.py
====================

Testet MS1 bei N=30, SIGN_WR=-1 mit mpmath-Praezisions-Diagonalisierung.

Frage: Ist der 21-fach entartete min-EW-Cluster in float64 ein Float-Artefakt
(echte Entartung ~10^{-10}) oder echte MS1-Verletzung?

Schritt 1: Lade W_R Cache (bereits mpmath-praezise berechnet, dann zu float64
reduziert). Berechne Wp und W02 in float64. Baue M in float64. Konvertiere zu
mpmath-Matrix. Diagonalisiere mit mp.eigsy und dps=50.

Wenn der min-EW-Cluster im mpmath-Ergebnis d=1 ist (echtes einfaches Minimum),
dann war die Entartung ein Float64-Artefakt der eigh-Routine. Wenn d>1 bleibt,
ist MS1 bei N=30 tatsaechlich verletzt.
"""
import numpy as np
import mpmath as mp
from pathlib import Path

mp.mp.dps = 50

# Import aus dem Haupt-Skript
import sys
sys.path.insert(0, str(Path(__file__).parent))
from dirichlet_ccm_mpmath import (
    LAMBDA, build_WR_matrix, build_W02_matrix, build_Wprime_matrix,
    chi_trivial
)

def main():
    N = 30
    lam = LAMBDA
    L = 2 * mp.log(lam)
    L_float = float(L)

    print(f"mpmath dps = {mp.mp.dps}")
    print(f"lambda = sqrt(14), N = {N}, L = {L_float:.10f}")

    # Lade WR aus Cache
    cache_path = Path(__file__).parent / "_cache" / f"WR_N{N}_L{L_float:.6f}_riemann.npy"
    if cache_path.exists():
        print(f"Lade WR aus Cache: {cache_path.name}")
        WR = np.load(cache_path)
    else:
        print("Baue WR...")
        WR = build_WR_matrix(N, L)

    print("Baue Wp (Prim)...")
    Wp = build_Wprime_matrix(N, L, lam, chi_trivial, 1)
    print("Baue W02 (Pol)...")
    W02 = build_W02_matrix(N, L)

    # CCM-Konvention: M = -WR - Wp + W02
    M_float = -WR - Wp + W02
    M_float = 0.5 * (M_float + M_float.T)  # Symmetrisieren

    print(f"M-Norm (float): {np.linalg.norm(M_float):.6f}")

    # Konvertiere zu mpmath-Matrix
    print(f"Konvertiere M zu mpmath (dps={mp.mp.dps})...")
    size = 2*N + 1
    M_mp = mp.matrix(size, size)
    for i in range(size):
        for j in range(size):
            M_mp[i, j] = mp.mpf(M_float[i, j])

    print(f"Diagonalisiere mit mp.eigsy ({size}x{size}, dps={mp.mp.dps})...")
    import time
    t0 = time.time()
    w_mp, V_mp = mp.eigsy(M_mp)
    t1 = time.time()
    print(f"  Fertig in {t1-t0:.1f}s")

    # Sortiere Eigenwerte (mpmath gibt unsortiert zurueck)
    w_arr = [float(w_mp[k]) for k in range(size)]
    order = sorted(range(size), key=lambda k: w_arr[k])

    print(f"\nEigenwerte sortiert (top 15, mpmath-praezise):")
    for rank in range(15):
        k = order[rank]
        w_val = w_mp[k]
        w_str_f = f"{float(w_val):+.10e}"
        # Differenz zum Minimum
        diff = float(w_mp[order[0]] - w_val) if rank > 0 else 0.0
        print(f"  rank={rank:2d}  k={k:3d}  w = {w_str_f}  delta_to_min = {-diff:+.6e}" if rank == 0
              else f"  rank={rank:2d}  k={k:3d}  w = {w_str_f}  delta_to_min = {(float(w_val)-float(w_mp[order[0]])):+.6e}")

    # Check MS1: ist der min-EW einfach?
    w_min = float(w_mp[order[0]])
    tol_range = [1e-20, 1e-15, 1e-10, 1e-6, 1e-4]
    print(f"\nMS1-Check: Cluster-Groesse bei verschiedenen Toleranzen:")
    for tol in tol_range:
        d = sum(1 for k in order if abs(float(w_mp[k]) - w_min) < tol)
        print(f"  tol = {tol:.0e}: d = {d}")

    # Extrahiere min-EV(en)
    print(f"\nMin-EV Parity-Check:")
    for rank in range(5):
        k = order[rank]
        v = [float(V_mp[i, k]) for i in range(size)]
        v_arr = np.array(v)
        is_even = np.allclose(v_arr, v_arr[::-1], atol=1e-6)
        is_odd = np.allclose(v_arr, -v_arr[::-1], atol=1e-6)
        sum_v = float(sum(v))
        sum_j = float(sum((i - N) * v[i] for i in range(size)))
        parity = 'EVEN' if is_even else ('ODD' if is_odd else 'mixed')
        print(f"  rank={rank}  parity={parity:5}  sum(xi)={sum_v:+.6e}  sum(j*xi)={sum_j:+.6e}")

    print("\n=== FAZIT ===")
    d_fine = sum(1 for k in order if abs(float(w_mp[k]) - w_min) < 1e-30)
    if d_fine == 1:
        print("MS1 ERFUELLT (Float64-Entartung war Artefakt)")
    elif d_fine <= 2:
        print("MS1 STRITTIG (2-facher Cluster bei hoher Prazision)")
    else:
        print(f"MS1 VERLETZT bei N={N} (d={d_fine} echte Entartung auch in mpmath)")

if __name__ == "__main__":
    main()
