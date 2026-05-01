#!/usr/bin/env python3
# coding: utf-8
"""
Session 8 Follow-up: Spectral-gap-index probe.

Hypothesis (Meta-Framework "Schatz" #2): The ratio
    I_chi := gamma_chi^(1) / sqrt(L),   L = log(lambda)
is a useful character-specific index that correlates with the
normalised empirical gap
    C_chi := gap_chi(N=600) / sqrt(L).

Two separate questions:
  (a) MAGNITUDE: does |C_chi| correlate with a simple function of
      gamma_chi^(1) (e.g. 1/gamma, 1/gamma^2, 1/gamma^3)?
  (b) SIGN: can sign(C_chi) be recovered from simple index data,
      possibly combined with L(1,chi) or root-number info?

Data inputs:
  _results/zeros_all_chars.json           - Riemann zero ordinates gamma_k^chi
  _results/ARCH_TERM_N600_SERVER.json     - empirical gaps (N=600, lambda=20000)
  _results/char_features_all.json         - L(1), -L'/L, etc.

Output:
  _results/SPECTRAL_GAP_INDEX.json        - raw correlations per index
  _results/SPECTRAL_GAP_INDEX.md          - human summary with ranking
"""
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RES = ROOT / "_results"

CHARS = ["chi_5", "chi_8", "chi_12", "chi_13", "chi_17",
         "chi_21", "chi_24", "chi_29", "chi_33", "chi_60"]


def pearson(x, y):
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    mx = x.mean(); my = y.mean()
    sx = x.std(); sy = y.std()
    if sx == 0 or sy == 0:
        return 0.0
    return float(((x - mx) * (y - my)).mean() / (sx * sy))


def load_data():
    with (RES / "zeros_all_chars.json").open(encoding="utf-8") as f:
        zeros = json.load(f)
    with (RES / "ARCH_TERM_N600_SERVER.json").open(encoding="utf-8") as f:
        arch = json.load(f)
    with (RES / "char_features_all.json").open(encoding="utf-8") as f:
        feats = {entry["name"]: entry for entry in json.load(f)}
    gaps = {entry["chi"]: entry["gap_galerkin"] for entry in arch["results"]}
    return zeros, gaps, feats


def main():
    lam = 20000
    L = math.log(lam)
    sqrtL = math.sqrt(L)
    print(f"[setup] lambda={lam}, L={L:.4f}, sqrt(L)={sqrtL:.4f}")

    zeros, gaps, feats = load_data()

    rows = []
    for name in CHARS:
        g1 = zeros[name][0]   # gamma_chi^(1)
        g2 = zeros[name][1] if len(zeros[name]) >= 2 else None
        gap = gaps[name]
        C = gap / sqrtL
        L1 = feats[name]["L1"]
        logD = feats[name]["log_D"]
        rows.append({
            "chi": name,
            "D": feats[name]["D"],
            "gamma1": g1,
            "gamma2": g2,
            "gap": gap,
            "C_chi": C,
            "abs_C": abs(C),
            "sign_C": 1 if C > 0 else (-1 if C < 0 else 0),
            "L1": L1,
            "log_D": logD,
            "I_ratio": g1 / sqrtL,
            "inv_gamma": 1.0 / g1,
            "inv_gamma2": 1.0 / (g1 * g1),
            "inv_gamma3": 1.0 / (g1 ** 3),
            "R_chi": 1.0 / (g1 ** 3) + (1.0 / (g2 ** 3) if g2 else 0.0),
        })

    print(f"\n{'chi':<8} {'D':>3} {'gamma1':>7} {'gap':>10} {'C_chi':>10} {'I_ratio':>8} {'R_chi':>9}")
    for r in rows:
        print(f"{r['chi']:<8} {r['D']:>3} {r['gamma1']:>7.3f} {r['gap']:>+10.5f} "
              f"{r['C_chi']:>+10.5f} {r['I_ratio']:>8.3f} {r['R_chi']:>9.5f}")

    # ---- Correlations ----
    C = np.array([r["C_chi"] for r in rows])
    abs_C = np.array([r["abs_C"] for r in rows])
    sign_C = np.array([r["sign_C"] for r in rows])

    idx_features = {
        "gamma1": "gamma_chi^(1)",
        "inv_gamma": "1/gamma^(1)",
        "inv_gamma2": "1/gamma^(1)^2",
        "inv_gamma3": "1/gamma^(1)^3 (Session 8 R_chi single-zero)",
        "R_chi": "R_chi = sum 1/gamma^(k)^3 (k<=2)",
        "I_ratio": "gamma^(1)/sqrt(L)",
        "L1": "L(1,chi)",
        "log_D": "log D",
    }

    print("\n--- Correlations vs C_chi (signed) ---")
    signed_corr = {}
    for key, label in idx_features.items():
        vals = np.array([r[key] for r in rows])
        R = pearson(vals, C)
        signed_corr[key] = R
        print(f"  R(C_chi, {label:<40}) = {R:+.4f}  R^2 = {R*R:.4f}")

    print("\n--- Correlations vs |C_chi| (magnitude) ---")
    abs_corr = {}
    for key, label in idx_features.items():
        vals = np.array([r[key] for r in rows])
        R = pearson(vals, abs_C)
        abs_corr[key] = R
        print(f"  R(|C_chi|, {label:<40}) = {R:+.4f}  R^2 = {R*R:.4f}")

    # ---- Test combined features for sign ----
    print("\n--- Combined features for sign recovery ---")
    # Hypothesis: sign(C) could depend on (gamma1 * L1) or similar
    for combo_name, combo_fn in [
        ("gamma1 * L1", lambda r: r["gamma1"] * r["L1"]),
        ("L1 - gamma1/5", lambda r: r["L1"] - r["gamma1"] / 5.0),
        ("gamma1 / L1", lambda r: r["gamma1"] / r["L1"]),
        ("log(gamma1)", lambda r: math.log(r["gamma1"])),
        ("(gamma1 - 3.0)", lambda r: r["gamma1"] - 3.0),
        ("inv_gamma3 * sign(L1-0.7)",
         lambda r: (1.0 / r["gamma1"] ** 3) * (1 if r["L1"] > 0.7 else -1)),
    ]:
        vals = np.array([combo_fn(r) for r in rows])
        R = pearson(vals, C)
        sign_match = sum(1 for r, v in zip(rows, vals)
                         if np.sign(v) == np.sign(r["C_chi"])
                         and r["C_chi"] != 0)
        print(f"  R(C_chi, {combo_name:<40}) = {R:+.4f}, sign match: {sign_match}/10")

    # ---- Ranking ----
    ranked = sorted(abs_corr.items(), key=lambda kv: abs(kv[1]), reverse=True)
    print("\n--- Top predictors of |C_chi| ---")
    for key, R in ranked[:4]:
        print(f"  {idx_features[key]:<42} R = {R:+.4f}, R^2 = {R*R:.4f}")

    # ---- Save JSON and MD ----
    summary = {
        "parameters": {"lambda": lam, "L": L, "sqrt_L": sqrtL},
        "rows": rows,
        "signed_correlations": signed_corr,
        "abs_correlations": abs_corr,
    }
    with (RES / "SPECTRAL_GAP_INDEX.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    md = RES / "SPECTRAL_GAP_INDEX.md"
    with md.open("w", encoding="utf-8") as f:
        f.write("# Spectral-Gap-Index Probe\n\n")
        f.write("**Datum:** 2026-04-16\n")
        f.write("**Skript:** `_scripts/spectral_gap_index.py`\n")
        f.write("**Motivation:** Meta-Framework 'verborgene Schoenheit' #2 - "
                "Verhaeltnis gamma^(1)/sqrt(L) als charakter-spezifischer "
                "Gap-Index.\n\n")
        f.write(f"## Parameter\n- lambda = {lam}, L = {L:.4f}, "
                f"sqrt(L) = {sqrtL:.4f}\n\n")
        f.write("## Daten pro Charakter\n\n")
        f.write("| chi | D | gamma^(1) | gap | C_chi | I_ratio=g/sqrtL | R_chi=1/g^3+1/g2^3 |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|\n")
        for r in rows:
            f.write(f"| {r['chi']} | {r['D']} | {r['gamma1']:.3f} | "
                    f"{r['gap']:+.5f} | {r['C_chi']:+.5f} | "
                    f"{r['I_ratio']:.3f} | {r['R_chi']:.5f} |\n")
        f.write("\n## Signed correlations (C_chi)\n\n")
        f.write("| Feature | Pearson R | R^2 |\n|---|---:|---:|\n")
        for key, label in idx_features.items():
            R = signed_corr[key]
            f.write(f"| {label} | {R:+.4f} | {R*R:.4f} |\n")
        f.write("\n## Magnitude correlations (|C_chi|)\n\n")
        f.write("| Feature | Pearson R | R^2 |\n|---|---:|---:|\n")
        for key, label in idx_features.items():
            R = abs_corr[key]
            f.write(f"| {label} | {R:+.4f} | {R*R:.4f} |\n")
        f.write("\n## Top-4 Predictors of |C_chi|\n\n")
        f.write("| Rank | Feature | R | R^2 |\n|:---:|---|---:|---:|\n")
        for i, (key, R) in enumerate(ranked[:4], 1):
            f.write(f"| {i} | {idx_features[key]} | {R:+.4f} | {R*R:.4f} |\n")

    print(f"\n[write] {md}")
    print(f"[write] {RES / 'SPECTRAL_GAP_INDEX.json'}")


if __name__ == "__main__":
    main()
