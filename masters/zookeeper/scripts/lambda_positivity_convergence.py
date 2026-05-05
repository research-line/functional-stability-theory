"""
lambda_positivity_convergence.py — Konvergiert QW zu positiver Semidefinitheit?

Test der "Konsensus-Hypothese": Mehr Primzahlen (groesseres lambda)
sollten QW weniger negativ machen, d.h. w_min -> 0 fuer lambda -> inf.

Weil-Kriterium: QW >= 0  <=>  RH.

Fuer jedes lambda:
  - Anzahl Primzahlen im Konsensus (bis lambda^2)
  - w_min (kleinster Eigenwert von QW_even)
  - Anzahl negativer Eigenwerte
  - Clustergroesse
  - F(gamma_1) fuer besten Cluster-EV
"""

import os
import sys
import time
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

N = int(os.environ.get("N", 30))
DPS = int(os.environ.get("DPS", 50))

GAMMA1 = 14.134725141734693


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


def analyze_lambda(lam, N, dps):
    mp.mp.dps = dps
    L_mp = 2 * mp.log(mp.mpf(lam))
    size = 2 * N + 1
    dim = N + 1
    k_max = int(mp.exp(L_mp))
    n_primes = count_primes_up_to(k_max)

    t0 = time.time()

    W02 = build_W02_matrix_mp(N, L_mp)
    WR = build_WR_matrix_mp(N, L_mp, parity="even", verbose=False)
    Wp = build_Wprime_matrix_mp(N, L_mp, chi_trivial, 1, verbose=False)

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

    w_min = float(w[0])
    w_max = float(w[dim - 1])
    n_neg = sum(1 for k in range(dim) if float(w[k]) < -1e-15)

    cluster = [k for k in range(dim) if abs(float(w[k]) - w_min) < 1e-10]
    n_cluster = len(cluster)

    # W02 Eigenwert
    w02_evals, _ = diagonalize_mp(W02_even, verbose=False)
    C_tilde = float(w02_evals[dim - 1])

    # F(gamma_1) fuer besten Cluster-EV (letzter im Cluster = groesster Gap)
    k_best = cluster[-1]
    q_best = mp.matrix(dim, 1)
    for r in range(dim):
        q_best[r, 0] = Q[r, k_best]
    F_g1 = evaluate_F_full(q_best, U_even, N, L_mp, GAMMA1)

    dt = time.time() - t0

    return {
        "lam": lam,
        "L": float(L_mp),
        "k_max": k_max,
        "n_primes": n_primes,
        "w_min": w_min,
        "w_max": w_max,
        "C_tilde": C_tilde,
        "n_neg": n_neg,
        "n_cluster": n_cluster,
        "F_g1": F_g1,
        "k_best": k_best,
        "dt": dt,
    }


def main():
    lambdas = [3, 5, 7, 10, 15, 20, 30, 50]

    # Optionale Erweiterung
    if os.environ.get("EXTENDED"):
        lambdas.extend([75, 100])

    print(f"Lambda-Positivitaets-Konvergenz: N={N}, dps={DPS}")
    print(f"Frage: Konvergiert w_min -> 0 fuer lambda -> inf?")
    print(f"Weil-Kriterium: QW >= 0  <=>  RH")
    print("=" * 100)
    print(f"  {'lam':>5s}  {'L':>7s}  {'k_max':>6s}  {'#primes':>7s}  "
          f"{'w_min':>16s}  {'C_tilde':>10s}  {'n_neg':>5s}  "
          f"{'cluster':>7s}  {'|F(g1)|':>12s}  {'time':>6s}")
    print("-" * 100)

    results = []
    for lam in lambdas:
        r = analyze_lambda(lam, N, DPS)
        results.append(r)
        print(f"  {r['lam']:5.0f}  {r['L']:7.4f}  {r['k_max']:6d}  {r['n_primes']:7d}  "
              f"{r['w_min']:+16.10e}  {r['C_tilde']:+10.4e}  {r['n_neg']:5d}  "
              f"{r['n_cluster']:7d}  {abs(r['F_g1']):12.3e}  {r['dt']:5.1f}s")

    print("\n" + "=" * 100)
    print("\nKonvergenz-Analyse:")

    # w_min Trend
    print(f"\n  {'lam':>5s}  {'w_min':>16s}  {'w_min/C_tilde':>14s}  {'Delta w_min':>14s}")
    for i, r in enumerate(results):
        ratio = r["w_min"] / r["C_tilde"] if r["C_tilde"] != 0 else 0
        delta = r["w_min"] - results[i - 1]["w_min"] if i > 0 else 0
        print(f"  {r['lam']:5.0f}  {r['w_min']:+16.10e}  {ratio:+14.6f}  "
              f"{delta:+14.6e}")

    # Negativitaets-Trend
    print(f"\n  Negativitaets-Trend:")
    print(f"  {'lam':>5s}  {'n_neg':>5s}  {'n_neg/dim':>9s}  {'w_min':>16s}")
    for r in results:
        frac = r["n_neg"] / (N + 1)
        print(f"  {r['lam']:5.0f}  {r['n_neg']:5d}  {frac:9.2%}  {r['w_min']:+16.10e}")

    # Zusammenfassung
    w_min_first = results[0]["w_min"]
    w_min_last = results[-1]["w_min"]
    direction = "WENIGER NEGATIV" if w_min_last > w_min_first else "NEGATIVER"
    print(f"\n  Ergebnis: w_min ging von {w_min_first:+.6e} (lam={results[0]['lam']}) "
          f"auf {w_min_last:+.6e} (lam={results[-1]['lam']})")
    print(f"  Trend: {direction}")
    if w_min_last > w_min_first:
        print(f"  => Konsistent mit Konvergenz zu QW >= 0 (Weil-Positivitaet / RH)")
    else:
        print(f"  => NICHT konsistent mit einfacher Konvergenz")
        print(f"     (moeglicherweise nicht-monoton oder N-abhaengig)")


if __name__ == "__main__":
    main()
