"""
c2cd_abel_partial_sum_diagnostic.py - Direkte Berechnung der Abel-Teilsummen

Wir beobachten in c2cd_prime_resonance_test.py:
  - Hoehere Off-Cluster Eigenwerte u_n (n>=3) sitzen NAHE bei log(p^k)/L Positionen
  - Aber b_n ist DOMINIERT von Cluster-Leckage (Moden 0,1,2)
  - Hoehere Moden haben |b_n| < 1e-7

Frage: Werden die Soll-Aussage Lambda_m^{off} und Theta_m^{off} wirklich nahe
bei 1, oder sind sie strukturell viel kleiner?

Diese Diagnostik berechnet:
  B^L_{m-1} := sum_{n<=m-1} b_n
  C^R_m    := sum_{n>=m} b_n
  Lambda_m := (s_m - u_0) * B^L_{m-1} / (s_m - u_{m-1})
  Theta_m  := (s_m - u_0) * C^R_m   / (u_m - s_m)

und vergleicht mit:
  - Cluster-Leckage Anteil: nur Moden mit d_n < g_*
  - Echte Off-Cluster (d_n > g_*): das was die Mertens-Diagnostik wuerde

Das C2by-Microcluster-Konzept besagt: die echten Off-Cluster-Moden tragen
viel weniger als die Soll-Aussage erlaubt, sobald die Cluster-Leckage
abgespalten ist.

Autor: LG (Opus 4.7, scheduled-task "kreativer-forscher")
Datum: 2026-05-01
"""

import os
import sys
import time
import json

import numpy as np
from pathlib import Path

sys.stdout.reconfigure(line_buffering=True)
sys.path.insert(0, str(Path(__file__).parent))

DPS = int(os.environ.get("DPS", 30))
os.environ["DPS"] = str(DPS)


def primes_up_to(n_max):
    if n_max < 2:
        return []
    sieve = np.ones(n_max + 1, dtype=bool)
    sieve[0] = sieve[1] = False
    for i in range(2, int(np.sqrt(n_max)) + 1):
        if sieve[i]:
            sieve[i*i::i] = False
    return np.where(sieve)[0].tolist()


def prime_powers_up_to(n_max, k_max=8):
    primes = primes_up_to(int(n_max))
    out = []
    for p in primes:
        pk = p
        k = 1
        while pk <= n_max:
            out.append((pk, p, k, np.log(p)))
            k += 1
            if k > k_max:
                break
            pk = p ** k
            if pk > n_max:
                break
    out.sort()
    return out


def diagnose_partial_sums(lam, N):
    """Compute Abel partial sums and bound diagnostics."""
    import mpmath as mp
    mp.mp.dps = DPS
    print(f"\n{'='*78}")
    print(f"  C2cd ABEL PARTIAL SUM DIAGNOSTIC")
    print(f"  lambda = {lam},  N = {N},  DPS = {DPS}")
    print(f"{'='*78}")

    from c2q_first_shell_dominance import build_operators
    t0 = time.time()
    ops = build_operators(lam, N)
    print(f"  build_operators: {time.time()-t0:.0f}s")

    alpha = ops["alpha"]
    w_min = ops["w_min"]
    w_arr = ops["w_arr"]
    c_arr = ops["c_arr"]
    alpha_a = ops["alpha_a"]
    cl_A = ops["cl_A"]
    noncl_A = ops["noncl_A"]
    t_vals = ops["t_vals"]  # T-Eigenwerte (= s_m Brunnen-Positionen)

    L = 2 * np.log(lam)
    print(f"\n  L = 2 log(lambda) = {L:.6f}")
    print(f"  u_0 (Pol-Cluster) = {w_min:.10f}")

    # Off-Cluster sortiert nach Eigenwert
    noncl_sorted = sorted(noncl_A, key=lambda a: w_arr[a])
    u_off = np.array([w_arr[a] for a in noncl_sorted])
    b_off = np.array([c_arr[a] * alpha_a[a] / alpha for a in noncl_sorted])
    M = len(u_off)

    # Distanz d_n = u_n - u_0
    d_off = u_off - w_min

    # T-Eigenwerte (Brunnen-Positionen s_m)
    t_sorted = np.sort(t_vals)
    print(f"\n  M off-cluster eigenvalues: {M}")
    print(f"  T eigenvalues range: [{t_sorted.min():.4f}, {t_sorted.max():.4f}]")

    # Bestimme Microcluster via C2by-Definition
    # alpha_a^2 = beta_n. Cluster-Massen:
    beta_off = np.array([alpha_a[a] ** 2 for a in noncl_sorted])
    cl_mass_off = np.array([alpha_a[a] ** 2 for a in cl_A]).sum()

    print(f"\n  Massenverteilung:")
    print(f"    Cluster (cl_A): mass = {cl_mass_off:.6f}")
    print(f"    First 5 off-cluster beta_n: {beta_off[:5]}")
    print(f"    Sum off-cluster beta_n: {beta_off.sum():.6e}")

    # Microcluster: erste Off-Cluster Moden mit d_n < g_eff (z.B. d < 1)
    micro_threshold = 0.5  # u-coord, anpassbar
    micro_idx = [n for n in range(M) if d_off[n] < micro_threshold]
    print(f"\n  Microcluster (d_n < {micro_threshold}): {len(micro_idx)} modes")
    for n in micro_idx:
        print(f"    n={n}: d_n={d_off[n]:.6e}, b_n={b_off[n]:.4e}")

    # Berechne Abel-Partialsummen
    cumsum_left = np.cumsum(b_off)  # B^L_m = sum_{n<=m} b_n
    sum_total = cumsum_left[-1]
    cumsum_right = sum_total - cumsum_left + b_off  # C^R_m = sum_{n>=m} b_n

    print(f"\n  Abel-Partialsummen (full):")
    print(f"    Total sum b_n: {sum_total:.6e}")
    print(f"    max |B^L|:     {np.max(np.abs(cumsum_left)):.6e}")
    print(f"    max |C^R|:     {np.max(np.abs(cumsum_right)):.6e}")

    # Now apply C2by-microcluster cleanup: b_n nur fuer d_n >= g_*
    g_star = 0.5  # threshold (anpassbar; in C2by ~5-6 fuer lambda=3)
    bulk_idx = [n for n in range(M) if d_off[n] >= g_star]
    leak_idx = [n for n in range(M) if d_off[n] < g_star]
    b_bulk = b_off.copy()
    b_leak = b_off.copy()
    for n in leak_idx:
        b_bulk[n] = 0.0
    for n in bulk_idx:
        b_leak[n] = 0.0

    print(f"\n  Microcluster-Bereinigung (g_* = {g_star}):")
    print(f"    Bulk modes (d_n >= g_*): {len(bulk_idx)}")
    print(f"    Leakage modes (d_n < g_*): {len(leak_idx)}")
    print(f"    Sum |b_bulk|: {np.abs(b_bulk).sum():.6e}")
    print(f"    Sum |b_leak|: {np.abs(b_leak).sum():.6e}")
    print(f"    Sum |B_bulk|: {np.max(np.abs(np.cumsum(b_bulk))):.6e}")
    print(f"    Sum |B_leak|: {np.max(np.abs(np.cumsum(b_leak))):.6e}")

    # Berechne Lambda_m und Theta_m fuer jedes Brunnen s_m
    print(f"\n  Soll-Aussage Diagnostik (Lambda, Theta) pro Brunnen-Position s_m:")
    print(f"  {'m':>3} {'s_m':>10} {'B_left':>14} {'C_right':>14} {'Lambda':>14} {'Theta':>14} {'max':>14}")
    print(f"  {'-'*100}")

    # s_m sind die Eigenwerte von T zwischen u-Werten
    # Filter: nur s_m im "regulaeren" Bereich (zwischen u_min+epsilon und u_max-epsilon)
    s_in_range = [s for s in t_sorted if s > u_off.min() and s < u_off.max()]

    table_rows = []
    max_lambda = 0.0
    max_theta = 0.0
    bulk_max_lambda = 0.0
    bulk_max_theta = 0.0

    for s in s_in_range:
        # Finde m: s zwischen u_{m-1} und u_m
        # u_off ist sortiert. Finde Index, sodass u_off[m-1] < s < u_off[m]
        idx_below = int(np.searchsorted(u_off, s) - 1)
        idx_above = int(np.searchsorted(u_off, s))
        if idx_below < 0 or idx_above >= M:
            continue
        u_below = u_off[idx_below]
        u_above = u_off[idx_above]
        gap_L = s - u_below
        gap_R = u_above - s
        d_m = s - w_min  # (s_m - u_0)

        # Volle Berechnung
        B_left = float(np.sum(b_off[:idx_above]))  # B^L_{m-1} = sum_{n<m} b_n
        C_right = float(np.sum(b_off[idx_above:]))  # C^R_m = sum_{n>=m} b_n
        Lambda_m = abs(d_m * B_left / gap_L) if gap_L > 0 else float('inf')
        Theta_m = abs(d_m * C_right / gap_R) if gap_R > 0 else float('inf')
        max_m = max(Lambda_m, Theta_m)

        # Bulk-only (Microcluster cleaned)
        B_left_bulk = float(np.sum(b_bulk[:idx_above]))
        C_right_bulk = float(np.sum(b_bulk[idx_above:]))
        Lambda_bulk = abs(d_m * B_left_bulk / gap_L) if gap_L > 0 else float('inf')
        Theta_bulk = abs(d_m * C_right_bulk / gap_R) if gap_R > 0 else float('inf')

        table_rows.append({
            "s_m": float(s), "m": idx_above,
            "u_below": float(u_below), "u_above": float(u_above),
            "gap_L": float(gap_L), "gap_R": float(gap_R), "d_m": float(d_m),
            "B_left": B_left, "C_right": C_right,
            "Lambda_m": Lambda_m, "Theta_m": Theta_m,
            "max_m": max_m,
            "Lambda_bulk": Lambda_bulk, "Theta_bulk": Theta_bulk,
        })
        max_lambda = max(max_lambda, Lambda_m)
        max_theta = max(max_theta, Theta_m)
        bulk_max_lambda = max(bulk_max_lambda, Lambda_bulk)
        bulk_max_theta = max(bulk_max_theta, Theta_bulk)

        if max_m > 1e-3 or len(table_rows) <= 5:
            print(f"  {idx_above:3d} {s:10.4f} {B_left:14.4e} {C_right:14.4e} "
                  f"{Lambda_m:14.4e} {Theta_m:14.4e} {max_m:14.4e}")

    print(f"\n  GESAMT-DIAGNOSTIK:")
    print(f"    Max Lambda_m (full):     {max_lambda:.4e}")
    print(f"    Max Theta_m (full):      {max_theta:.4e}")
    print(f"    Max max(Lambda,Theta):   {max(max_lambda, max_theta):.4e}")
    print(f"")
    print(f"    Max Lambda_m (bulk only):  {bulk_max_lambda:.4e}")
    print(f"    Max Theta_m (bulk only):   {bulk_max_theta:.4e}")
    print(f"    Max max(Lambda,Theta) bulk: {max(bulk_max_lambda, bulk_max_theta):.4e}")
    print(f"")
    print(f"    Soll: max(Lambda,Theta) < 1.")
    print(f"    Voll: {max(max_lambda, max_theta):.4e} {'OK' if max(max_lambda, max_theta) < 1 else 'FAIL'}")
    print(f"    Bulk: {max(bulk_max_lambda, bulk_max_theta):.4e} {'OK' if max(bulk_max_lambda, bulk_max_theta) < 1 else 'FAIL'}")

    # NEUE STRUKTURELLE BEOBACHTUNG: Quotient bulk/full
    if max(max_lambda, max_theta) > 0:
        ratio = max(bulk_max_lambda, bulk_max_theta) / max(max_lambda, max_theta)
        print(f"\n    Microcluster-Reduktion: {ratio:.4f}x ({(1-ratio)*100:.1f}% durch Cluster-Leckage)")

    # Prime-power-Lokalisierung: was tragen die "primnahen" Moden bei?
    n_max_search = max(100, int(2 * np.exp(d_off.max() * L)))
    prime_powers = prime_powers_up_to(n_max_search, k_max=5)
    if prime_powers:
        ppos_arr = np.array([np.log(pk)/L for (pk, _, _, _) in prime_powers])

        # Fuer jede Off-Cluster-Mode (idx >= 3): finde naechste Primzahlpotenz
        prime_aligned_b = 0.0
        prime_misaligned_b = 0.0
        prime_align_threshold = 0.05  # in u-coord
        for n in range(M):
            if d_off[n] < g_star:
                continue
            diffs = np.abs(ppos_arr - d_off[n])
            min_diff = float(np.min(diffs))
            if min_diff < prime_align_threshold:
                prime_aligned_b += abs(b_off[n])
            else:
                prime_misaligned_b += abs(b_off[n])
        print(f"\n  Prime-power-Lokalisierung im Bulk:")
        print(f"    |b_n| sum (within {prime_align_threshold} of prime power): {prime_aligned_b:.4e}")
        print(f"    |b_n| sum (outside):                                      {prime_misaligned_b:.4e}")

    return {
        "lam": lam, "N": N, "L": L, "u_0": float(w_min), "M": M,
        "max_lambda": float(max_lambda),
        "max_theta": float(max_theta),
        "bulk_max_lambda": float(bulk_max_lambda),
        "bulk_max_theta": float(bulk_max_theta),
        "n_microcluster": len(leak_idx),
        "n_bulk": len(bulk_idx),
        "table_rows": table_rows[:30],  # save first 30 only
    }


def main():
    all_results = []
    configs = [{"lam": 3.0, "N": 20}]

    for cfg in configs:
        try:
            res = diagnose_partial_sums(cfg["lam"], cfg["N"])
            all_results.append(res)
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"\n  FEHLER bei lambda={cfg['lam']}: {e}")

    out_path = Path(__file__).parent.parent / "_results" / "C2CD_ABEL_DIAGNOSTIC_2026-05-01.json"
    out_path.parent.mkdir(exist_ok=True, parents=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(all_results, f, indent=2, ensure_ascii=False)
    print(f"\n  Saved: {out_path}")
    return all_results


if __name__ == "__main__":
    main()
