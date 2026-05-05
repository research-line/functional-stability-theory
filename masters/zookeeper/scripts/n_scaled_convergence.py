"""
n_scaled_convergence.py — N-skalierte Weil-Positivitaets-Konvergenz.

Zentrale Frage: Konvergiert w_min -> 0 wenn N proportional zu L waechst?

Wenn ja: numerischer Beleg fuer Weil-Positivitaet (QW >= 0 im Limes).
Wenn nein: entweder Skalierung falsch, oder tieferes Problem.

N wird so gewaehlt, dass die Fourier-Aufloesung 2*pi*N/L > 50 bleibt
(genuegend Moden um gamma_5 ~ 33 aufzuloesen).
"""

import os
import sys
import time
import json
import mpmath as mp
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from dirichlet_ccm_fourier_mp import (
    build_WR_matrix_mp,
    build_W02_matrix_mp,
    build_Wprime_matrix_mp,
    chi_trivial,
    project_to_parity,
    diagonalize_mp,
    is_prime_power,
)

DPS = int(os.environ.get("DPS", 50))

GAMMA1 = 14.134725141734693
RIEMANN_ZEROS = [14.134725141734693, 21.022039638771554, 25.010857580145688,
                 30.424876125859513, 32.935061587739189]


def count_primes_up_to(k_max):
    count = 0
    for k in range(2, k_max + 1):
        p, m = is_prime_power(k)
        if p == k and m == 1:
            count += 1
    return count


def evaluate_F_full(q_k, U_even, N, L_mp, gamma):
    size = 2 * N + 1
    dim = N + 1
    xi = mp.matrix(size, 1)
    for row in range(size):
        s = mp.mpf(0)
        for r in range(dim):
            s += U_even[row, r] * q_k[r, 0]
        xi[row, 0] = s
    g = mp.mpf(gamma)
    F = mp.mpf(0)
    for idx in range(size):
        j = idx - N
        pole = 2 * mp.pi * j / L_mp
        F += xi[idx, 0] / (g - pole)
    return float(F)


def analyze(lam, N, dps):
    mp.mp.dps = dps
    L_mp = 2 * mp.log(mp.mpf(lam))
    size = 2 * N + 1
    dim = N + 1
    k_max = int(mp.exp(L_mp))
    n_primes = count_primes_up_to(k_max)

    t0 = time.time()
    print(f"  lambda={lam}, N={N}, dim={dim}, k_max={k_max}, #primes={n_primes}", flush=True)

    W02 = build_W02_matrix_mp(N, L_mp)
    WR = build_WR_matrix_mp(N, L_mp, parity="even", verbose=False)
    t_wr = time.time() - t0
    print(f"    WR fertig ({t_wr:.0f}s)", flush=True)

    Wp = build_Wprime_matrix_mp(N, L_mp, chi_trivial, 1, verbose=False)
    t_wp = time.time() - t0
    print(f"    W' fertig ({t_wp:.0f}s)", flush=True)

    dummy = mp.matrix(size, size)
    _, U_even = project_to_parity(dummy, N, parity="even")

    W02_even = U_even.T * W02 * U_even
    WR_even = U_even.T * WR * U_even
    Wp_even = U_even.T * Wp * U_even

    M_even = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            M_even[i, j] = W02_even[i, j] - WR_even[i, j] - Wp_even[i, j]
    for i in range(dim):
        for j in range(i + 1, dim):
            avg = (M_even[i, j] + M_even[j, i]) / 2
            M_even[i, j] = avg
            M_even[j, i] = avg

    w, Q = diagonalize_mp(M_even, verbose=False)
    t_diag = time.time() - t0
    print(f"    Diag fertig ({t_diag:.0f}s)", flush=True)

    w_min = float(w[0])
    w_max = float(w[dim - 1])
    n_neg = sum(1 for k in range(dim) if float(w[k]) < -1e-15)
    cluster = [k for k in range(dim) if abs(float(w[k]) - w_min) < 1e-10]
    n_cluster = len(cluster)

    w02_evals, _ = diagonalize_mp(W02_even, verbose=False)
    C_tilde = float(w02_evals[dim - 1])

    # F(gamma_1) fuer besten Cluster-EV
    k_best = cluster[-1]
    q_best = mp.matrix(dim, 1)
    for r in range(dim):
        q_best[r, 0] = Q[r, k_best]
    F_g1 = evaluate_F_full(q_best, U_even, N, L_mp, GAMMA1)

    # Cancellation-Analyse fuer k=0
    inv_sqrt2 = 1 / mp.sqrt(mp.mpf(2))
    u_tilde = mp.matrix(dim, 1)
    max_k_w02 = max(range(dim), key=lambda i: float(w02_evals[i]))
    w02_evecs_all = diagonalize_mp(W02_even, verbose=False)[1]
    for r in range(dim):
        u_tilde[r, 0] = w02_evecs_all[r, max_k_w02]

    q0 = mp.matrix(dim, 1)
    for r in range(dim):
        q0[r, 0] = Q[r, 0]
    uq0 = sum(float(u_tilde[r, 0] * q0[r, 0]) for r in range(dim))

    # Eval vector at gamma_1 (corrected with sqrt2)
    g = mp.mpf(GAMMA1)
    ue_g1 = float(u_tilde[0, 0]) / float(g)
    for n in range(1, N + 1):
        pole = 2 * mp.pi * n / L_mp
        ue_g1 += float(u_tilde[n, 0]) * float(inv_sqrt2 * 2 * g / (g**2 - pole**2))
    pol_term = ue_g1 * uq0

    F_g1_k0 = evaluate_F_full(q0, U_even, N, L_mp, GAMMA1)
    null_term = F_g1_k0 - pol_term
    if abs(F_g1_k0) > 1e-60 and abs(pol_term) > 1e-60:
        cancel_digits = -np.log10(abs(F_g1_k0) / max(abs(pol_term), abs(null_term)))
    else:
        cancel_digits = 0

    dt = time.time() - t0

    return {
        "lam": lam, "N": N, "dim": dim, "L": float(L_mp),
        "k_max": k_max, "n_primes": n_primes,
        "w_min": w_min, "w_max": w_max, "C_tilde": C_tilde,
        "n_neg": n_neg, "n_cluster": n_cluster,
        "F_g1": F_g1, "F_g1_k0": F_g1_k0,
        "pol_term": pol_term, "cancel_digits": cancel_digits,
        "k_best": k_best, "dt": dt,
    }


def main():
    # N skaliert mit L: N = max(30, round(14 * L))
    # So dass Fourier-Aufloesung 2*pi*N/L ~ 88 (konstant)
    configs = [
        (3, 30),       # L=2.2, Referenz
        (5, 45),       # L=3.2
        (7, 55),       # L=3.9
        (10, 65),      # L=4.6
        (15, 76),      # L=5.4
        (20, 84),      # L=6.0
        (30, 96),      # L=6.8
        (50, 110),     # L=7.8
    ]

    if os.environ.get("EXTENDED"):
        configs.extend([
            (75, 122),   # L=8.6
            (100, 130),  # L=9.2
            (150, 140),  # L=10.0
            (200, 148),  # L=10.6
        ])

    print(f"N-skalierte Weil-Positivitaets-Konvergenz")
    print(f"  DPS={DPS}, Skalierung: N ~ 14*L (Fourier-Aufloesung ~ 88)")
    print(f"  Frage: Konvergiert w_min -> 0 wenn N ~ L?")
    print("=" * 110)
    print(f"  {'lam':>5s} {'N':>4s} {'dim':>4s} {'L':>6s} {'#p':>5s}  "
          f"{'w_min':>16s} {'C_tilde':>10s} {'w_min/C':>9s} "
          f"{'n_neg':>5s} {'clust':>5s} {'|F(g1)|':>11s} {'cancel':>7s} {'time':>7s}")
    print("-" * 110)

    results = []
    for lam, N in configs:
        r = analyze(lam, N, DPS)
        results.append(r)
        ratio = r["w_min"] / r["C_tilde"] if r["C_tilde"] != 0 else 0
        print(f"  {r['lam']:5.0f} {r['N']:4d} {r['dim']:4d} {r['L']:6.3f} {r['n_primes']:5d}  "
              f"{r['w_min']:+16.10e} {r['C_tilde']:+10.4e} {ratio:+9.5f} "
              f"{r['n_neg']:5d} {r['n_cluster']:5d} {abs(r['F_g1']):11.3e} "
              f"{r['cancel_digits']:7.1f} {r['dt']:6.0f}s")

    print("\n" + "=" * 110)
    print("\nKonvergenz-Tabelle (N skaliert mit L):")
    print(f"  {'lam':>5s} {'N':>4s} {'N/L':>6s} {'w_min':>16s} {'w_min/C':>10s} {'cancel':>8s}")
    for r in results:
        ratio = r["w_min"] / r["C_tilde"]
        n_over_l = r["N"] / r["L"]
        print(f"  {r['lam']:5.0f} {r['N']:4d} {n_over_l:6.1f} "
              f"{r['w_min']:+16.10e} {ratio:+10.6f} {r['cancel_digits']:8.1f}")

    w_first = results[0]["w_min"]
    w_last = results[-1]["w_min"]
    print(f"\n  w_min: {w_first:+.6e} (lam={results[0]['lam']}) -> {w_last:+.6e} (lam={results[-1]['lam']})")
    if w_last > w_first:
        print(f"  TREND: WENIGER NEGATIV => Konvergenz Richtung QW >= 0")
    elif abs(w_last - w_first) < 0.1:
        print(f"  TREND: STABIL (kein klarer Trend)")
    else:
        print(f"  TREND: NEGATIVER => N-Skalierung reicht nicht")

    # JSON fuer spaetere Analyse
    out_path = Path(__file__).parent.parent / "_results" / "batch_beta" / "n_scaled_convergence.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  JSON gespeichert: {out_path}")


if __name__ == "__main__":
    main()
