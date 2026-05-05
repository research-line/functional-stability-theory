"""
c2cc_log_concavity_test.py — Log-Konkavitaet & Detailed-Balance Diagnose

Testet drei neue strukturelle Hypothesen, die im bisherigen Beweislog NICHT
abgedeckt sind, fuer das aktive Endproblem C2ao.1 / C2as.1:

  H1 (LOG-KONKAVITAET):    log beta_n^2 >= log beta_{n-1} + log beta_{n+1}
                          (klassisch: Newton-Ungleichung => Partialsum-Schranke)

  H1.2 (GEOMETRISCHER DECAY): beta_n <= C * r^n mit r < 1
                              (per linearem Fit von log beta_n gegen n)

  H3 (MERTENS-ANKER):     sum b_n * log(1/d_n) <= K * (-log beta_0)

Gegen welche Profile testen wir?
  beta_n := alpha_{a_n}^2     (kanonisches Pol-Gewicht aus C2at.3)
  b_n    := c_{a_n} alpha_{a_n} / alpha   (reales Off-Cluster-Profil aus C2y/C2z)

Eingangsdaten kommen aus build_operators() (siehe c2q_first_shell_dominance.py).

Status: EXPLORATIV. Skript falsifiziert oder verifiziert numerisch eine
Hypothese, beweist nichts.

Autor: LG (Opus 4.7, scheduled-task "kreativer-forscher")
Datum: 2026-05-01
"""

import os
import sys
import time

import numpy as np
from pathlib import Path

sys.stdout.reconfigure(line_buffering=True)

sys.path.insert(0, str(Path(__file__).parent))
from c2q_first_shell_dominance import build_operators

DPS = int(os.environ.get("DPS", 50))
CONFIGS = [
    {"lam": 3.0, "N": 30},
    {"lam": 5.0, "N": 55},
    {"lam": 7.0, "N": 77},
]


def analyze_profile(label, profile, n_show=12):
    """
    Analysiert eine Off-Cluster-Sequenz (positive Eintraege, sortiert) auf:
      - Log-Konkavitaet rho_n = 2 log(p_n) - log(p_{n-1}) - log(p_{n+1})
      - Geometrischer Decay (linearer Fit von log p_n)
      - Konzentration N_eff = exp(H) wo H = -sum p log p (normiert)
    """
    p = np.asarray(profile, dtype=float)
    p = p[p > 0]  # entferne nullwertige Eintraege
    M = len(p)

    if M < 3:
        return None

    log_p = np.log(p)

    # Log-Konkavitaet (rho_n >= 0 fuer 1 <= n <= M-2 in 0-indexed)
    rho = np.zeros(M)
    rho[0] = np.nan
    rho[-1] = np.nan
    for n in range(1, M - 1):
        rho[n] = 2 * log_p[n] - log_p[n - 1] - log_p[n + 1]

    rho_inner = rho[1:-1]
    n_concave = int(np.sum(rho_inner >= 0))
    n_convex = int(np.sum(rho_inner < 0))
    rho_min = float(np.min(rho_inner))
    rho_min_idx = int(np.argmin(rho_inner)) + 1  # zurueck zu n-Index

    # Geometrischer Decay-Fit: log p_n = a + b*n
    n_arr = np.arange(M, dtype=float)
    coef = np.polyfit(n_arr, log_p, 1)
    slope = float(coef[0])
    intercept = float(coef[1])
    fit_resid = log_p - (slope * n_arr + intercept)
    rms_resid = float(np.sqrt(np.mean(fit_resid ** 2)))
    r_geom = float(np.exp(slope))

    # Konzentration
    p_norm = p / p.sum()
    H = float(-np.sum(p_norm * np.log(p_norm)))
    N_eff = float(np.exp(H))

    print(f"\n{'-'*78}")
    print(f"  Profil: {label}  (M={M} Eintraege, alle > 0)")
    print(f"{'-'*78}")
    print(f"  Konzentration:    H = {H:.4f}, N_eff = {N_eff:.3f} ({100*N_eff/M:.1f}% von M)")
    print(f"  Geom. Decay:      r = exp(slope) = {r_geom:.4e}, intercept = {intercept:.3f}")
    print(f"                    Fit-Residuum (rms log) = {rms_resid:.3f}")
    print(f"  Log-Konkavitaet:  rho_n inneres Min = {rho_min:.4e} bei n = {rho_min_idx}")
    print(f"                    rho_n >= 0: {n_concave} / {n_concave + n_convex} (innere Indizes)")
    print(f"                    H1 voll log-konkav?  {'JA' if n_convex == 0 else 'NEIN'}")

    # Tail-Log-Konkavitaet: ab welchem n_0 ist rho_n >= 0 fuer alle n >= n_0?
    n0_tail = None
    for n_start in range(1, M - 1):
        if np.all(rho[n_start:M - 1] >= 0):
            n0_tail = n_start
            break
    if n0_tail is not None:
        print(f"  Tail-Log-Konkav:  ab n_0 = {n0_tail} (alle weiteren rho_n >= 0)")
    else:
        print(f"  Tail-Log-Konkav:  KEIN konsistenter Tail")

    # Erste paar Eintraege ausgeben
    show = min(n_show, M)
    print(f"\n  Erste {show} Eintraege:")
    print(f"  {'n':>3}  {'p_n':>14}  {'log p_n':>10}  {'rho_n':>12}  {'fit log':>10}")
    for n in range(show):
        rho_str = f"{rho[n]:12.4e}" if np.isfinite(rho[n]) else f"{'(rand)':>12}"
        fit_val = slope * n + intercept
        print(f"  {n:3d}  {p[n]:14.6e}  {log_p[n]:10.3f}  {rho_str}  {fit_val:10.3f}")

    return {
        "label": label,
        "M": M,
        "H": H,
        "N_eff": N_eff,
        "r_geom": r_geom,
        "intercept": intercept,
        "rms_resid": rms_resid,
        "rho_min": rho_min,
        "rho_min_idx": rho_min_idx,
        "n_concave": n_concave,
        "n_convex": n_convex,
        "fully_log_concave": (n_convex == 0),
        "tail_concave_from": n0_tail,
        "rho": rho.tolist(),
        "log_p": log_p.tolist(),
    }


def mertens_anchor_test(b_off, d_arr, beta_0):
    """
    Testet H3: sum |b_n| * log(1/d_n) <= K * (-log beta_0).

    Hier d_n = u_n - u_0 = Off-Cluster Pol-Distance.
    """
    mask = (np.abs(b_off) > 0) & (d_arr > 0)
    weighted_sum = float(np.sum(np.abs(b_off[mask]) * np.log(1.0 / d_arr[mask])))
    target = -np.log(beta_0) if beta_0 > 0 else float("inf")
    ratio = weighted_sum / target if target > 0 else float("inf")

    print(f"\n  H3 Mertens-Anker:")
    print(f"    sum |b_n| log(1/d_n)  = {weighted_sum:.6e}")
    print(f"    -log(beta_0)          = {target:.6e}  (beta_0 = {beta_0:.6f})")
    print(f"    Ratio K              = {ratio:.4e}")

    return {
        "weighted_sum": weighted_sum,
        "log_inv_beta0": target,
        "ratio_K": ratio,
    }


def detailed_balance_test(p_off, q_forward):
    """
    Testet H2: ist (p_n) stationaer fuer eine reversible Markov-Kette mit
    Vorwaertsraten q_n = q_forward[n] (z.B. g_n^R / Delta_n)?

    Detailed balance: p_n * q_n = p_{n+1} * q_{rev,n+1}.
    Fuer reversibilitaet aus den Daten: definiere q_{rev,n+1} := p_n*q_n / p_{n+1}.
    Pruefe Konsistenz: q_{rev,n} sollte glatt verlaufen (kein Sprung).
    """
    M = len(p_off)
    if M < 3 or len(q_forward) != M:
        return None
    q_rev = np.zeros(M)
    for n in range(1, M):
        if p_off[n] > 0:
            q_rev[n] = p_off[n - 1] * q_forward[n - 1] / p_off[n]
        else:
            q_rev[n] = float("inf")

    finite_mask = np.isfinite(q_rev[1:])
    if finite_mask.sum() < 2:
        return None

    log_q_rev = np.log(q_rev[1:][finite_mask])
    median = float(np.median(log_q_rev))
    spread = float(np.std(log_q_rev))

    print(f"\n  H2 Detailed Balance:")
    print(f"    median log q_rev = {median:.4f}  (= {np.exp(median):.4e})")
    print(f"    std log q_rev    = {spread:.4f}  (kleiner = glatter)")
    print(f"    Interpretation:  spread klein => approx. detailed balance")

    return {
        "median_log_q_rev": median,
        "std_log_q_rev": spread,
        "q_rev": q_rev.tolist(),
    }


def matched_check(b_off, alpha_a_arr, w_arr, w_min, t_vals, j_reg):
    """
    Falls H1 oder H1.2 gilt, leite daraus eine matched-side-mass-Schranke ab
    und vergleiche mit existierender C2AJ-Schranke.
    """
    noncl_sorted = sorted(range(len(alpha_a_arr)), key=lambda a: w_arr[a])
    # b_off ist bereits in noncl_sorted-Reihenfolge gegeben
    pass  # Detail-Vergleich ist bereits in c2ao_matched_sidemass_test.py


def run(lam, N):
    import mpmath as mp
    mp.mp.dps = DPS
    print(f"\n{'='*78}")
    print(f"  C2cc LOG-CONKAVITAET / DETAILED BALANCE / MERTENS-ANKER")
    print(f"  lambda = {lam},  N = {N}")
    print(f"{'='*78}")
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

    # Cluster-Residue beta_0 = sum alpha_a^2 ueber Cluster
    beta_0 = float(np.sum([alpha_a[a] ** 2 for a in cl_A]))

    # Sortiere Off-Cluster-Indizes nach Eigenwert
    noncl_sorted = sorted(noncl_A, key=lambda a: w_arr[a])

    # Profile aufbauen
    beta_arr = np.array([alpha_a[a] ** 2 for a in noncl_sorted])
    b_off = np.array([c_arr[a] * alpha_a[a] / alpha for a in noncl_sorted])
    u_off = np.array([w_arr[a] for a in noncl_sorted])
    d_off = u_off - w_min  # Delta_n = u_n - u_0

    # Wrapper: nur positive Eintraege im b-Profil (die meisten sind positiv)
    b_pos_mask = b_off > 0
    n_neg = int(np.sum(~b_pos_mask))

    print(f"\n  Off-Cluster-Diagnostik:")
    print(f"    Anzahl Off-Cluster-Moden M = {len(beta_arr)}")
    print(f"    beta_0 (Cluster-Residue)    = {beta_0:.6f}")
    print(f"    sum beta_n (Off-Cluster)    = {beta_arr.sum():.6f}  (sollte 1-beta_0 = {1-beta_0:.6f})")
    print(f"    sum b_n                     = {b_off.sum():.6e}")
    print(f"    sum |b_n|                   = {np.abs(b_off).sum():.6e}")
    print(f"    n_neg im b-Profil           = {n_neg} / {len(b_off)} ({100*n_neg/len(b_off):.1f}%)")

    # H1 / H1.2 / Konzentration: Profile analysieren
    res_beta = analyze_profile("beta_n = alpha_a^2 (kanonisch)", beta_arr)
    res_b_abs = analyze_profile("|b_n| = |c_a alpha_a / alpha|", np.abs(b_off))

    # H3 (Mertens-Anker)
    mertens = mertens_anchor_test(b_off, d_off, beta_0)

    # H2 (Detailed Balance) — auf normierten beta_n
    p_norm = beta_arr / beta_arr.sum()
    # Vorwaertsrate q_n = (u_{n+1} - u_n) / Delta_n_to_normalize
    # Vereinfacht: q_n = 1 / Delta_n als Surrogat fuer Hazard-Rate
    q_forward = 1.0 / np.where(d_off > 0, d_off, 1.0)
    db_res = detailed_balance_test(p_norm, q_forward)

    out = {
        "lam": lam, "N": N,
        "beta_0": beta_0,
        "M": len(beta_arr),
        "n_neg_b": n_neg,
        "profile_beta": res_beta,
        "profile_b_abs": res_b_abs,
        "mertens": mertens,
        "detailed_balance": db_res,
    }
    return out


def main():
    all_results = {}
    for cfg in CONFIGS:
        try:
            res = run(cfg["lam"], cfg["N"])
            all_results[f"lam_{cfg['lam']}"] = res
        except Exception as e:
            print(f"\n  FEHLER bei lambda={cfg['lam']}: {e}")
            import traceback
            traceback.print_exc()

    # Zusammenfassung
    print(f"\n\n{'='*78}")
    print(f"  ZUSAMMENFASSUNG")
    print(f"{'='*78}")
    print(f"\n  Hypothese H1 (Volle Log-Konkavitaet von beta_n):")
    for key, res in all_results.items():
        if res and res.get("profile_beta"):
            pb = res["profile_beta"]
            label = "JA" if pb["fully_log_concave"] else "NEIN"
            tail = pb["tail_concave_from"]
            tail_str = f"Tail ab n0={tail}" if tail else "kein Tail"
            print(f"    {key}: {label}  (rho_min={pb['rho_min']:.2e} @ n={pb['rho_min_idx']}, {tail_str})")

    print(f"\n  Hypothese H1.2 (Geometrischer Decay r < 1):")
    for key, res in all_results.items():
        if res and res.get("profile_beta"):
            pb = res["profile_beta"]
            print(f"    {key}: r = {pb['r_geom']:.4e}, fit-rms (log) = {pb['rms_resid']:.3f}")

    print(f"\n  Hypothese H3 (Mertens-Anker, K = sum b_n log(1/d_n) / -log beta_0):")
    for key, res in all_results.items():
        if res and res.get("mertens"):
            m = res["mertens"]
            print(f"    {key}: K = {m['ratio_K']:.4e}")

    return all_results


if __name__ == "__main__":
    results = main()
    # Optional: dump als JSON
    try:
        import json
        out_path = Path(__file__).parent.parent / "_results" / "C2CC_LOG_CONCAVITY_2026-05-01.json"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        # Numpy-Werte konvertieren
        def conv(o):
            if isinstance(o, (np.integer,)):
                return int(o)
            if isinstance(o, (np.floating,)):
                return float(o)
            if isinstance(o, np.ndarray):
                return o.tolist()
            raise TypeError
        with open(out_path, "w") as f:
            json.dump(results, f, indent=2, default=conv)
        print(f"\n  JSON gespeichert: {out_path}")
    except Exception as e:
        print(f"\n  JSON-Dump uebersprungen: {e}")
