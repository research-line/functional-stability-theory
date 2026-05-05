"""
c2ce_prime_square_dominance_test.py — Prime-Square-Channel Dominance (C2ce)

Folgesession nach C2cd (kreativer-forscher 2026-05-02).

Hypothese C2ce.1 (Prime-Square Channel Law):
  Fuer Off-Cluster-Eigenwerte u_n im Bulk-Bereich (d_n > g_*) erfuellen die
  Massen |b_n| die strukturelle Skalierung
      |b_n| approx G(N,L) / sqrt(p_n^{k_n})  fuer k_n >= 2,
      |b_n| << G/sqrt(p_n)                    fuer k_n = 1
  wobei (p_n, k_n) das nahegelegenste Prime-Power zu u_n ist und G(N,L)
  eine Galerkin-Trunkierungs-Konstante.

Datenquelle:
  - Existierender JSON: _results/C2CD_PRIME_RESONANCE_2026-05-01.json
  - Optional: Erneuter build_operators-Lauf bei groesserem lambda/N

Tests:
  T1: Klassifiziere bulk modes nach k.
  T2: Fitte log|b_n| vs (d_n, log(p^k), 1) per Multilineare Regression
       und vergleiche Modelle:
         M0: log|b| = const + alpha*d                   (Distance-Decay)
         M1: log|b| = const + alpha*d + beta*log(p^k)   (Distance + Power)
         M2: log|b| = const + beta*log(p^k)             (Pure Power)
         M3: log|b| = const, sep. fuer k=1, k=2, k>=3   (Multiplicity-Group)
  T3: Verifikation: |b_{p^k}|*sqrt(p^k) ~ const fuer k>=2.
  T4: Mertens-Anker: M_off vs G * sum_{p^k <= X, k>=2} 1/sqrt(p^k).
  T5: Suppressions-Faktor: |b_{p}|/(G/sqrt(p)) fuer k=1.

Autor: LG (Opus 4.7, scheduled-task "kreativer-forscher" Session 3, 2026-05-02)
"""

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(line_buffering=True)
sys.path.insert(0, str(Path(__file__).parent))

DPS = int(os.environ.get("DPS", 30))
os.environ["DPS"] = str(DPS)

# Default: load existing C2CD JSON
DEFAULT_JSON = Path(__file__).parent.parent / "_results" / "C2CD_PRIME_RESONANCE_2026-05-01.json"


def primes_up_to(n_max):
    if n_max < 2:
        return []
    sieve = np.ones(n_max + 1, dtype=bool)
    sieve[0] = sieve[1] = False
    for i in range(2, int(np.sqrt(n_max)) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    return np.where(sieve)[0].tolist()


def mertens_partial_sum(x_max, k):
    """Berechne sum_{p^k <= x_max} 1/sqrt(p^k) for given k."""
    if k == 1:
        prime_bound = int(x_max)
    else:
        prime_bound = int(x_max ** (1.0 / k)) + 1
    primes = primes_up_to(prime_bound)
    return sum(1.0 / np.sqrt(p ** k) for p in primes if p ** k <= x_max)


def analyze_bulk_modes(results_B, threshold_d=0.5):
    """Klassifiziere bulk modes (d > threshold) nach k und teste Hypothesen."""
    bulk = [r for r in results_B if r["u_n_rel"] > threshold_d and r["b_n"] != 0]
    print(f"\n  Bulk modes (d > {threshold_d}): {len(bulk)}")

    by_k = {}
    for r in bulk:
        by_k.setdefault(r["k"], []).append(r)

    print(f"\n  Verteilung nach k:")
    total_mass = sum(abs(r["b_n"]) for r in bulk)
    for k in sorted(by_k):
        modes = by_k[k]
        mass = sum(abs(r["b_n"]) for r in modes)
        print(
            f"    k={k}: {len(modes):>2} modes, sum |b_n| = {mass:.3e} ({100*mass/total_mass:5.1f}%)"
        )
    print(f"    total: {total_mass:.3e}")

    # T3: |b|*sqrt(p^k) const for k>=2?
    print(f"\n  T3: Skalierungsgesetz |b|*sqrt(p^k) fuer k>=2:")
    for k in sorted(by_k):
        if k < 2:
            continue
        vals = np.array([abs(r["b_n"]) * np.sqrt(r["pk"]) for r in by_k[k]])
        if len(vals) >= 2:
            log_std = np.std(np.log(vals))
            ratio = vals.max() / vals.min()
            mean_g = float(np.mean(vals))
            print(
                f"    k={k}: G_emp = mean(|b|*sqrt(pk)) = {mean_g:.3e}, "
                f"log_std={log_std:.3f}, max/min={ratio:.2f}"
            )

    # Aggregated G across k>=2
    arr_pp = [r for r in bulk if r["k"] >= 2]
    if arr_pp:
        G_arr = np.array([abs(r["b_n"]) * np.sqrt(r["pk"]) for r in arr_pp])
        G_mean = float(np.mean(G_arr))
        G_std = float(np.std(G_arr))
        print(
            f"\n    Aggregat (k>=2): G = {G_mean:.3e} +/- {G_std:.3e}, "
            f"CV = {G_std/G_mean:.3f}, log_std = {np.std(np.log(G_arr)):.3f}"
        )

    # T5: Suppression for k=1
    print(f"\n  T5: Suppressions-Faktor fuer k=1 (Primzahlen):")
    if arr_pp:
        G_mean = float(np.mean([abs(r["b_n"]) * np.sqrt(r["pk"]) for r in arr_pp]))
    else:
        G_mean = None
    for r in by_k.get(1, []):
        if G_mean:
            G_pred = G_mean / np.sqrt(r["pk"])
            actual = abs(r["b_n"])
            suppression = G_pred / actual if actual > 0 else float("inf")
            print(
                f"    p={r['p']:>5}: |b|={actual:.3e}, G/sqrt(p)={G_pred:.3e}, "
                f"suppression = {suppression:.1f}x"
            )

    # T4: Mertens-Anker
    print(f"\n  T4: Mertens-Anker — Vorhersage M_off via Sum 1/sqrt(p^k) fuer k>=2:")
    if arr_pp:
        x_max = max(r["pk"] for r in bulk)
        G = float(np.mean([abs(r["b_n"]) * np.sqrt(r["pk"]) for r in arr_pp]))
        M_pred_k2 = G * mertens_partial_sum(x_max, 2)
        M_pred_k3 = G * mertens_partial_sum(x_max, 3)
        M_pred_k4 = G * mertens_partial_sum(x_max, 4)
        M_pred = M_pred_k2 + M_pred_k3 + M_pred_k4
        M_obs = total_mass
        print(f"    x_max = {x_max}")
        print(f"    k=2 contrib: {M_pred_k2:.3e}")
        print(f"    k=3 contrib: {M_pred_k3:.3e}")
        print(f"    k=4 contrib: {M_pred_k4:.3e}")
        print(f"    M_pred (sum): {M_pred:.3e}")
        print(f"    M_obs       : {M_obs:.3e}")
        if M_pred > 0:
            print(f"    Ratio pred/obs: {M_pred/M_obs:.2f}")

    return {
        "total_mass": total_mass,
        "by_k_mass": {
            k: sum(abs(r["b_n"]) for r in by_k[k]) for k in sorted(by_k)
        },
        "G_emp": G_mean if arr_pp else None,
        "G_log_std_kpp": float(
            np.std(
                np.log(
                    [abs(r["b_n"]) * np.sqrt(r["pk"]) for r in bulk if r["k"] >= 2]
                )
            )
        )
        if any(r["k"] >= 2 for r in bulk)
        else None,
    }


def regression_analysis(results_B, threshold_d=0.5):
    """T2: Multilineare Regression."""
    bulk = [r for r in results_B if r["u_n_rel"] > threshold_d and r["b_n"] != 0]
    if len(bulk) < 4:
        return None

    d = np.array([r["u_n_rel"] for r in bulk])
    log_pk = np.log(np.array([r["pk"] for r in bulk]))
    k_arr = np.array([r["k"] for r in bulk])
    log_b = np.log(np.array([abs(r["b_n"]) for r in bulk]))

    print(f"\n  T2: Modellvergleich (R^2):")

    # M0: log|b| = a*d + c
    if len(d) >= 2:
        slope, intercept = np.polyfit(d, log_b, 1)
        pred = slope * d + intercept
        r2_m0 = 1 - np.sum((log_b - pred) ** 2) / np.sum((log_b - np.mean(log_b)) ** 2)
        print(f"    M0 (log|b| = a*d + c):                R^2 = {r2_m0:.3f}")

    # M2: log|b| = b*log(p^k) + c
    if np.std(log_pk) > 0:
        slope2, intercept2 = np.polyfit(log_pk, log_b, 1)
        pred2 = slope2 * log_pk + intercept2
        r2_m2 = 1 - np.sum((log_b - pred2) ** 2) / np.sum(
            (log_b - np.mean(log_b)) ** 2
        )
        print(
            f"    M2 (log|b| = b*log(p^k) + c):         R^2 = {r2_m2:.3f}, "
            f"slope = {slope2:.3f} (alpha = {-slope2:.3f})"
        )

    # M1: log|b| = a*d + b*log(p^k) + c (multilineare Regression)
    X = np.column_stack([d, log_pk, np.ones_like(d)])
    coef, residuals, rank, _ = np.linalg.lstsq(X, log_b, rcond=None)
    pred1 = X @ coef
    r2_m1 = 1 - np.sum((log_b - pred1) ** 2) / np.sum((log_b - np.mean(log_b)) ** 2)
    print(
        f"    M1 (log|b| = a*d + b*log(p^k) + c):  R^2 = {r2_m1:.3f}, "
        f"a={coef[0]:.3f}, b={coef[1]:.3f}, c={coef[2]:.3f}"
    )

    # M3: separate const fuer k=1, k=2, k>=3
    print(f"    M3 (Multiplicity-Group):")
    for k_target in sorted(set(k_arr.tolist())):
        mask = k_arr == k_target
        if mask.sum() >= 1:
            mean = float(np.mean(log_b[mask]))
            std = float(np.std(log_b[mask])) if mask.sum() >= 2 else 0.0
            print(
                f"      k={k_target}: n={mask.sum()}, mean log|b|={mean:.3f}, std={std:.3f}, "
                f"geom_mean |b|={np.exp(mean):.3e}"
            )

    # M4 (best known): log|b| = -0.5*log(p^k) + const, applied to k>=2
    mask_pp = k_arr >= 2
    if mask_pp.sum() >= 2:
        log_b_pp = log_b[mask_pp]
        log_pk_pp = log_pk[mask_pp]
        # Force slope = -0.5 (theoretical Mertens-style), fit only intercept
        residual = log_b_pp + 0.5 * log_pk_pp  # = log(|b|*sqrt(p^k))
        c_emp = float(np.mean(residual))
        std_residual = float(np.std(residual))
        print(
            f"    M4 (k>=2: log|b| = -0.5*log(p^k) + c, fixed slope):"
            f" c = {c_emp:.3f}, residual_std = {std_residual:.3f}"
        )
        print(
            f"        => G = exp(c) = {np.exp(c_emp):.3e}, |b|*sqrt(p^k) tightness {std_residual:.3f}"
        )


def main():
    json_path = DEFAULT_JSON
    print(f"Loading {json_path}")
    with open(json_path) as f:
        data = json.load(f)

    out = []
    for cfg in data:
        print(f"\n{'='*78}")
        print(f"  C2ce ANALYSIS: lambda={cfg['lam']}, N={cfg['N']}, L={cfg['L']:.4f}")
        print(f"{'='*78}")

        result = analyze_bulk_modes(cfg["results_B"])
        regression_analysis(cfg["results_B"])

        out.append(
            {
                "lam": cfg["lam"],
                "N": cfg["N"],
                "L": cfg["L"],
                **result,
            }
        )

    out_path = Path(__file__).parent.parent / "_results" / "C2CE_PRIME_SQUARE_DOMINANCE.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\n  -> Results: {out_path}")


if __name__ == "__main__":
    main()
