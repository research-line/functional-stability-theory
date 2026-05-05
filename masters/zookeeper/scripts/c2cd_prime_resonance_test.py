"""
c2cd_prime_resonance_test.py - Primzahl-Resonanz-Hypothese (PR1+PR2)

Folgefrage aus C2cc (kreativer-forscher 2026-05-01):
  Die Massenstruktur ist arithmetisch (primzahl-getrieben), nicht analytisch
  glatt. Frage: Sind die b_n-Spitzen explizit an Indizes mit
    u_n = log(prime)/L  (oder log(prime^k)/L)
  lokalisiert?

Hypothesen:

  PR0a (Prime Localization, schwach):
    Fuer jeden Off-Cluster-Modus u_n existiert ein "naechstes" prime power
    p^k mit log(p^k)/L moeglichst nahe bei u_n. Das relative Residuum
       r_n := |u_n - log(p_n^{k_n})/L| / Delta_typical
    ist klein fuer hohe-Masse-Moden.

  PR0b (Prime Mass Profile):
    Die b_n-Werte folgen einer Mertens-Struktur:
       |b_n| ~ const * Lambda(p_n^{k_n}) / sqrt(p_n^{k_n})
    (von Mangoldt-gewichtete Primzahl-Reihe).

  PR1 (Mertens Partial Sum Bound):
    Falls PR0a+b gelten, dann gilt fuer Off-Cluster Abel-Teilsummen
       B_{m-1}^{off} = sum_{u_n < s_m} b_n
       ~ C * sum_{p^k <= e^{s_m * L}} Lambda(p^k)/sqrt(p^k)
       ~ C * 2*sqrt(e^{s_m L}) / log(e^{s_m L})  (Mertens 2nd theorem),
    und die Soll-Aussage |delta_m| < 1 folgt aus klassischer Zahlentheorie.

Methodik:
  1. build_operators(lambda, N) - berechne CCM-Spektrum
  2. Fuer jedes off-cluster u_n: finde naechstes prime power, Residuum r_n
  3. Korreliere Top-|b_n|-Moden mit Top-Lambda(p^k)/sqrt(p^k)-Werten
  4. Plot/Tabelle: u_n vs. log(p)/L, b_n vs. Lambda/sqrt(p)
  5. Falls Korrelation hoch: Mertens-Anker-Test direkt auf prime-indizierter Reihe.

Status: EXPLORATIV. Skript prueft strukturell, ob die b_n als
        primzahl-gewichtete Reihe geschrieben werden koennen.
        Falsifizierbar: wenn die u_n zufaellig verteilt sind oder
        die b_n nicht mit Lambda(p^k)/sqrt(p^k) korrelieren.

Autor: LG (Opus 4.7, scheduled-task "kreativer-forscher" Session 2)
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

CONFIGS = [
    {"lam": 3.0, "N": 20},  # Schneller initial test
]


def primes_up_to(n_max):
    """Sieb des Eratosthenes bis n_max."""
    if n_max < 2:
        return []
    sieve = np.ones(n_max + 1, dtype=bool)
    sieve[0] = sieve[1] = False
    for i in range(2, int(np.sqrt(n_max)) + 1):
        if sieve[i]:
            sieve[i*i::i] = False
    return np.where(sieve)[0].tolist()


def prime_powers_up_to(n_max, k_max=10):
    """Gibt Liste (n, p, k, Lambda(n)) fuer alle Primzahlpotenzen p^k <= n_max."""
    primes = primes_up_to(int(n_max))
    out = []
    for p in primes:
        pk = p
        k = 1
        while pk <= n_max:
            out.append((pk, p, k, np.log(p)))  # Lambda(p^k) = log p
            k += 1
            if k > k_max:
                break
            pk = p ** k
            if pk > n_max:
                break
    out.sort()
    return out  # list of (n, p, k, log_p)


def find_nearest_prime_power(u_target, L, prime_powers, log_p_table):
    """Findet das prime power n mit log(n)/L am naechsten zu u_target."""
    if not log_p_table:
        return None, None, None, None, None

    log_p_arr = np.array([log(L) for _, _, _, log in [(0,0,0,lp) for lp in log_p_table]])
    # log(n)/L fuer alle n
    u_candidates = np.array([np.log(n) / L for (n, _, _, _) in prime_powers])
    if len(u_candidates) == 0:
        return None, None, None, None, None

    diffs = np.abs(u_candidates - u_target)
    idx = int(np.argmin(diffs))
    n, p, k, lambda_val = prime_powers[idx]
    return n, p, k, lambda_val, float(diffs[idx])


def analyze_prime_resonance(lam, N):
    """Pruefe PR0a (Lokalisierung) und PR0b (Massenprofil)."""
    import mpmath as mp
    mp.mp.dps = DPS
    print(f"\n{'='*78}")
    print(f"  C2cd PRIME RESONANCE ANALYSIS")
    print(f"  lambda = {lam},  N = {N},  DPS = {DPS}")
    print(f"{'='*78}")

    from c2q_first_shell_dominance import build_operators
    t0 = time.time()
    ops = build_operators(lam, N)
    build_time = time.time() - t0
    print(f"  build_operators: {build_time:.0f}s")

    alpha = ops["alpha"]
    w_min = ops["w_min"]
    w_arr = ops["w_arr"]
    c_arr = ops["c_arr"]
    alpha_a = ops["alpha_a"]
    cl_A = ops["cl_A"]
    noncl_A = ops["noncl_A"]

    L = 2 * np.log(lam)
    print(f"\n  L = 2 log(lambda) = {L:.6f}")
    print(f"  u_0 (Pol-Cluster) = {w_min:.10f}")
    print(f"  exp(u_0 * L)      = {np.exp(w_min * L):.6f}")

    # Off-Cluster Eigenwerte
    noncl_sorted = sorted(noncl_A, key=lambda a: w_arr[a])
    u_off = np.array([w_arr[a] for a in noncl_sorted])
    b_off = np.array([c_arr[a] * alpha_a[a] / alpha for a in noncl_sorted])
    beta_n = np.array([alpha_a[a] ** 2 for a in noncl_sorted])

    M = len(u_off)
    print(f"\n  M = {M} off-cluster modes")

    # Wandle u_n in "natuerliche Skala" um: x_n = exp(u_n * L)
    x_off = np.exp(u_off * L)
    print(f"\n  exp(u_n * L) range: [{x_off.min():.6f}, {x_off.max():.6f}]")

    # Erweiterte natuerliche Skala: u_n koennten log(p^k)/L positiv sein
    # ODER auch durch Spiegel-Symmetrie auf der negativen Seite (Pol bei u_0 negativ)
    # Pruefe beide Skalen:
    # Option A: x_n = exp(u_n * L)  (direkt)
    # Option B: x_n = exp((u_n - u_0) * L) (relativ zum Pol)

    x_off_rel = np.exp((u_off - w_min) * L)
    print(f"  exp((u_n - u_0) * L) range: [{x_off_rel.min():.6f}, {x_off_rel.max():.6f}]")

    # Generiere prime power Tabelle hinreichend gross
    n_max_search = max(100, int(2 * x_off.max()), int(2 * x_off_rel.max()))
    prime_powers = prime_powers_up_to(n_max_search, k_max=5)
    print(f"  Prime powers <= {n_max_search}: {len(prime_powers)}")

    # Fuer jeden u_n: finde naechste Primzahlpotenz in beiden Skalen
    print(f"\n  Tabelle: u_n -> naechste Primzahlpotenz (Skala A: direkt)")
    print(f"  {'n':>3} {'u_n':>12} {'x_n':>12} {'p^k':>10} {'log(p^k)/L':>12} {'res':>10} "
          f"{'b_n':>14} {'Lam/sqrt(p^k)':>14}")
    print(f"  {'-'*98}")

    results_A = []
    for n in range(min(M, 20)):
        u_n = u_off[n]
        x_n = x_off[n]

        # Suche unter prime powers
        u_candidates = np.array([np.log(pk) / L for (pk, _, _, _) in prime_powers])
        diffs = np.abs(u_candidates - u_n)
        idx = int(np.argmin(diffs))
        pk, p, k, log_p = prime_powers[idx]
        u_pk = np.log(pk) / L
        residual = float(diffs[idx])
        lam_per_sqrt = log_p / np.sqrt(pk)

        results_A.append({
            "n": n, "u_n": float(u_n), "x_n": float(x_n),
            "pk": int(pk), "p": int(p), "k": int(k),
            "u_pk": float(u_pk), "residual": residual,
            "b_n": float(b_off[n]), "lambda_over_sqrt": float(lam_per_sqrt),
            "beta_n": float(beta_n[n]),
        })

        marker = ""
        if abs(b_off[n]) > 1e-5:
            marker = " <- HIGH MASS"

        print(f"  {n:3d} {u_n:12.6f} {x_n:12.4f} {pk:10d} {u_pk:12.6f} {residual:10.4e} "
              f"{b_off[n]:14.6e} {lam_per_sqrt:14.6e}{marker}")

    # Skala B: relativ zum Pol u_0
    print(f"\n  Tabelle: u_n -> naechste Primzahlpotenz (Skala B: rel. zum Pol)")
    print(f"  {'n':>3} {'u_n-u_0':>12} {'x_n_rel':>12} {'p^k':>10} {'log(p^k)/L':>12} {'res':>10} "
          f"{'b_n':>14}")
    print(f"  {'-'*84}")

    results_B = []
    for n in range(min(M, 20)):
        u_n_rel = u_off[n] - w_min
        if u_n_rel <= 0:
            continue
        x_n_rel = x_off_rel[n]

        u_candidates = np.array([np.log(pk) / L for (pk, _, _, _) in prime_powers])
        diffs = np.abs(u_candidates - u_n_rel)
        idx = int(np.argmin(diffs))
        pk, p, k, log_p = prime_powers[idx]
        u_pk = np.log(pk) / L
        residual = float(diffs[idx])

        results_B.append({
            "n": n, "u_n_rel": float(u_n_rel), "x_n_rel": float(x_n_rel),
            "pk": int(pk), "p": int(p), "k": int(k),
            "u_pk": float(u_pk), "residual": residual,
            "b_n": float(b_off[n]),
        })

        marker = ""
        if abs(b_off[n]) > 1e-5:
            marker = " <- HIGH MASS"

        print(f"  {n:3d} {u_n_rel:12.6f} {x_n_rel:12.4f} {pk:10d} {u_pk:12.6f} {residual:10.4e} "
              f"{b_off[n]:14.6e}{marker}")

    # Korrelations-Analyse: |b_n| vs. Lambda(p^k)/sqrt(p^k)
    print(f"\n  Korrelations-Analyse (Skala A):")
    bn_arr = np.array([abs(r["b_n"]) for r in results_A])
    lam_arr = np.array([r["lambda_over_sqrt"] for r in results_A])
    res_arr = np.array([r["residual"] for r in results_A])

    if np.std(bn_arr) > 0 and np.std(lam_arr) > 0:
        # Pearson, Spearman
        from scipy.stats import pearsonr, spearmanr
        try:
            r_p, p_p = pearsonr(bn_arr, lam_arr)
            r_s, p_s = spearmanr(bn_arr, lam_arr)
            print(f"    Pearson  |b_n| vs Lambda/sqrt(p^k): r = {r_p:.4f}, p = {p_p:.4f}")
            print(f"    Spearman |b_n| vs Lambda/sqrt(p^k): r = {r_s:.4f}, p = {p_s:.4f}")
        except Exception as e:
            print(f"    scipy.stats Korrelation fehlgeschlagen: {e}")

    # Top |b_n| Moden: bei welchen Primzahlpotenzen?
    print(f"\n  Top-5 |b_n| Moden (Skala A):")
    sorted_by_bn = sorted(results_A, key=lambda r: -abs(r["b_n"]))[:5]
    for r in sorted_by_bn:
        print(f"    n={r['n']:3d}  u_n={r['u_n']:.4f}  p^k={r['pk']}  res={r['residual']:.3e}  "
              f"|b_n|={abs(r['b_n']):.3e}  Lambda/sqrt(p^k)={r['lambda_over_sqrt']:.3e}")

    # Top |b_n| Moden (Skala B):
    print(f"\n  Top-5 |b_n| Moden (Skala B, rel. zum Pol):")
    sorted_by_bn_B = sorted(results_B, key=lambda r: -abs(r["b_n"]))[:5]
    for r in sorted_by_bn_B:
        print(f"    n={r['n']:3d}  u_n-u_0={r['u_n_rel']:.4f}  p^k={r['pk']}  res={r['residual']:.3e}  "
              f"|b_n|={abs(r['b_n']):.3e}")

    # Histogramm der Residuen
    print(f"\n  Residuen-Statistik (Skala A):")
    print(f"    median residual: {np.median(res_arr):.4e}")
    print(f"    max residual:    {np.max(res_arr):.4e}")
    print(f"    min residual:    {np.min(res_arr):.4e}")

    # Charakteristische Skala: typischer Abstand zwischen aufeinanderfolgenden
    # Primzahlpotenzen in der u-Skala
    if len(prime_powers) >= 2:
        u_pp = np.array([np.log(pk) / L for (pk, _, _, _) in prime_powers if pk > 1])
        u_pp.sort()
        spacings = np.diff(u_pp)
        print(f"    typ. spacing log(p^k)/L: median={np.median(spacings):.4f}, "
              f"min={spacings.min():.4f}")

    return {
        "lam": lam, "N": N, "L": L, "u_0": float(w_min),
        "M": M, "build_time": build_time,
        "results_A": results_A,
        "results_B": results_B,
        "bn_lam_corr_pearson": (
            float(pearsonr(bn_arr, lam_arr)[0]) if len(bn_arr) >= 3
            and np.std(bn_arr) > 0 and np.std(lam_arr) > 0 else None
        ),
        "residual_stats": {
            "median": float(np.median(res_arr)),
            "max": float(np.max(res_arr)),
            "min": float(np.min(res_arr)),
        },
    }


def main():
    all_results = []
    for cfg in CONFIGS:
        try:
            res = analyze_prime_resonance(cfg["lam"], cfg["N"])
            all_results.append(res)
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"\n  FEHLER bei lambda={cfg['lam']}: {e}")

    # JSON-Output
    out_path = Path(__file__).parent.parent / "_results" / "C2CD_PRIME_RESONANCE_2026-05-01.json"
    out_path.parent.mkdir(exist_ok=True, parents=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(all_results, f, indent=2, ensure_ascii=False)
    print(f"\n  Saved to: {out_path}")

    return all_results


if __name__ == "__main__":
    try:
        from scipy.stats import pearsonr  # noqa
    except ImportError:
        print("scipy.stats nicht verfuegbar, verwende numpy-corr")
        # fallback definition
        def pearsonr(x, y):
            return (float(np.corrcoef(x, y)[0,1]), 0.0)
    main()
