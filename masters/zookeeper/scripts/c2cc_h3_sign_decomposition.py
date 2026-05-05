"""
c2cc_h3_sign_decomposition.py — Robustheits-Diagnostik fuer Mertens-Anker (H3)

Pruef-Frage (advisor-getrieben):
  Ist H3 ("sum |b_n| log(1/d_n) ~ 0.4% von -log beta_0") ein robuster
  Strukturzwang oder ein Cancellation-Artefakt?

Diagnostik:
  S_signed = sum_n |b_n| * log(1/d_n)              # original, signed kernel
  S_abs    = sum_n |b_n| * |log(1/d_n)|            # echte L1-Norm
  ratio    = S_abs / |S_signed|

  ratio ~ 1   => robust (alle Beitraege gleiches Vorzeichen)
  ratio >> 1  => Cancellation, K aus Glueck nicht Struktur

Spaltung in Beitraege mit positivem vs. negativem Kernel:
  S_pos = sum_{d_n < 1} |b_n| * log(1/d_n)         (positive)
  S_neg = sum_{d_n > 1} |b_n| * log(1/d_n)         (negative)

Konfiguration: lambda=3, N=30, DPS=30 (schneller als DPS=50, fuer
Diagnostik ausreichend; betrifft nur Eigenwert-Schaerfe).

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

# DPS=30 statt 50 fuer schnelleren Build, ausreichend fuer Diagnostik
os.environ["DPS"] = "30"
from c2q_first_shell_dominance import build_operators


def diagnostic(lam, N):
    print(f"\n{'='*78}")
    print(f"  H3 SIGN-DECOMPOSITION DIAGNOSTIK")
    print(f"  lambda = {lam},  N = {N},  DPS = 30")
    print(f"{'='*78}")
    t0 = time.time()
    ops = build_operators(lam, N)
    print(f"  build_operators: {time.time()-t0:.0f}s\n")

    alpha = ops["alpha"]
    w_min = ops["w_min"]
    w_arr = ops["w_arr"]
    c_arr = ops["c_arr"]
    alpha_a = ops["alpha_a"]
    cl_A = ops["cl_A"]
    noncl_A = ops["noncl_A"]

    beta_0 = float(np.sum([alpha_a[a] ** 2 for a in cl_A]))
    target = -np.log(beta_0)

    noncl_sorted = sorted(noncl_A, key=lambda a: w_arr[a])
    b_off = np.array([c_arr[a] * alpha_a[a] / alpha for a in noncl_sorted])
    u_off = np.array([w_arr[a] for a in noncl_sorted])
    d_off = u_off - w_min

    # Mask: nur Off-Cluster mit d > 0
    mask = d_off > 0
    b_abs = np.abs(b_off[mask])
    d = d_off[mask]
    log_inv_d = np.log(1.0 / d)
    abs_log_inv_d = np.abs(log_inv_d)

    # S_signed = original signed sum
    S_signed = float(np.sum(b_abs * log_inv_d))
    # S_abs = L1 norm
    S_abs = float(np.sum(b_abs * abs_log_inv_d))
    # S_pos / S_neg split by kernel sign
    pos_mask = log_inv_d > 0  # d < 1
    neg_mask = log_inv_d < 0  # d > 1
    S_pos = float(np.sum(b_abs[pos_mask] * log_inv_d[pos_mask]))
    S_neg = float(np.sum(b_abs[neg_mask] * log_inv_d[neg_mask]))  # negative

    n_pos = int(pos_mask.sum())
    n_neg = int(neg_mask.sum())
    M = len(b_abs)

    print(f"  Off-Cluster-Statistik:")
    print(f"    M (Off-Cluster Moden, d>0)   = {M}")
    print(f"    w_min                          = {w_min:.6f}")
    print(f"    d_min, d_max                   = {d.min():.4f}, {d.max():.4f}")
    print(f"    Kernel pos (d_n < 1):          {n_pos} Moden")
    print(f"    Kernel neg (d_n > 1):          {n_neg} Moden")

    print(f"\n  Mertens-Sum Decomposition:")
    print(f"    S_signed = sum |b_n| log(1/d_n)        = {S_signed:+.6e}")
    print(f"    S_abs    = sum |b_n| |log(1/d_n)|      = {S_abs:+.6e}")
    print(f"    S_pos    = sum_{{d<1}} |b_n| log(1/d_n) = {S_pos:+.6e}")
    print(f"    S_neg    = sum_{{d>1}} |b_n| log(1/d_n) = {S_neg:+.6e}")
    print(f"    Identitaet S_pos + S_neg = {S_pos + S_neg:+.6e} (sollte = S_signed)")

    if abs(S_signed) > 0:
        ratio_robust = S_abs / abs(S_signed)
        print(f"\n  Robustheits-Diagnose:")
        print(f"    ratio S_abs / |S_signed|   = {ratio_robust:.4f}")
        if ratio_robust < 1.5:
            print(f"    => H3 ist NICHT durch Cancellation getrieben")
        elif ratio_robust < 5:
            print(f"    => H3 hat moderate Cancellation")
        else:
            print(f"    => H3 ist STARK durch Cancellation getrieben")
            print(f"       (S_pos und S_neg heben sich gegenseitig auf)")

    print(f"\n  Vergleich gegen -log beta_0:")
    print(f"    -log(beta_0)              = {target:+.6e}  (beta_0 = {beta_0:.6f})")
    print(f"    K_signed = S_signed/(-log beta_0) = {S_signed/target if target>0 else float('nan'):+.4e}")
    print(f"    K_abs    = S_abs/(-log beta_0)    = {S_abs/target if target>0 else float('nan'):+.4e}")

    # Weiterer Test: was, wenn wir d_n statt log(1/d_n) verwenden?
    S_d = float(np.sum(b_abs * d))
    S_inv_d = float(np.sum(b_abs / d))
    print(f"\n  Alternative Skalen:")
    print(f"    sum |b_n| * d_n            = {S_d:.6e}")
    print(f"    sum |b_n| * (1/d_n)        = {S_inv_d:.6e}")

    # Ist die "wahre" Summe Sum |b_n|/d_n besser bei -log beta_0?
    if S_inv_d > 0:
        ratio_inv = S_inv_d / target
        print(f"    K_inv = (sum |b|/d) / (-log beta_0) = {ratio_inv:+.4e}")

    # Mehrere Top-Beitraege ausgeben
    contrib = b_abs * abs_log_inv_d
    top5_idx = np.argsort(-contrib)[:5]
    print(f"\n  Top 5 Beitraege zu S_abs:")
    print(f"    {'n':>3}  {'u_n':>10}  {'d_n':>10}  {'log(1/d_n)':>12}  {'|b_n|':>14}  {'|b_n||log..|':>14}")
    for n in top5_idx:
        print(f"    {n:3d}  {u_off[mask][n]:10.4f}  {d[n]:10.4f}  {log_inv_d[n]:+12.4f}  "
              f"{b_abs[n]:14.6e}  {contrib[n]:14.6e}")

    return {
        "lam": lam, "N": N,
        "M": M,
        "beta_0": beta_0,
        "target_minus_log_beta0": target,
        "S_signed": S_signed,
        "S_abs": S_abs,
        "S_pos": S_pos,
        "S_neg": S_neg,
        "n_pos_kernel": n_pos,
        "n_neg_kernel": n_neg,
        "ratio_robust_S_abs_over_signed": ratio_robust if abs(S_signed) > 0 else None,
        "K_signed": S_signed/target if target > 0 else None,
        "K_abs": S_abs/target if target > 0 else None,
        "S_d": S_d,
        "S_inv_d": S_inv_d,
    }


if __name__ == "__main__":
    print(f"H3 SIGN-DECOMPOSITION (advisor-getrieben)")
    print(f"Frage: Ist H3 robust oder Cancellation-Artefakt?\n")

    res = diagnostic(lam=3.0, N=30)

    print(f"\n{'='*78}")
    print(f"  ZUSAMMENFASSUNG")
    print(f"{'='*78}")
    print(f"  ratio S_abs / |S_signed| = {res['ratio_robust_S_abs_over_signed']:.4f}")
    if res['ratio_robust_S_abs_over_signed'] is not None:
        if res['ratio_robust_S_abs_over_signed'] < 1.5:
            print(f"  => H3 ist STRUKTURELL ROBUST")
        elif res['ratio_robust_S_abs_over_signed'] < 5:
            print(f"  => H3 hat MODERATE Cancellation")
        else:
            print(f"  => H3 ist CANCELLATION-getrieben (NICHT strukturell robust)")
    print(f"  K_signed = {res['K_signed']:+.4e}")
    print(f"  K_abs    = {res['K_abs']:+.4e}")

    try:
        import json
        out_path = Path(__file__).parent.parent / "_results" / "C2CC_H3_SIGN_DECOMP_2026-05-01.json"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as f:
            json.dump(res, f, indent=2, default=lambda o: float(o) if hasattr(o, '__float__') else str(o))
        print(f"\n  JSON gespeichert: {out_path}")
    except Exception as e:
        print(f"\n  JSON skipped: {e}")
