#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MS1 — Newton-Power-Sum-Bilanz der grad-gepaarten Normalisierung (Prüfschritt 2026-06-10).

Kontext (Zookeeper, Route Z / Determinanten-Kanal, Programmschritt 2):
  Die Bulk-Quellen-Triangulation (Proof-Note MS1_DET_ZERO_COUNTING_HARNESS §5.1) hat gezeigt,
  dass die z-Abhängigkeit des Bulk-Faktors ausschließlich vom Pol-Kompensations-Polynom und
  der Xi-Hüllkurve getragen wird. Prüfverdikt (2026-06-10): kanonische
  grad-gepaarte Normalisierung mit ERZWUNGENEM K (keine Wahl-Marge, Lambda=0-rein) plus
  Newton-Power-Sum-Diagnose, welche Momente die Paarung verfehlt.

Objekt:
  P(z) := F_N(z) * prod_{j=-N..N} (z - p_j)   — Polynom Grad 2N (exakt, RH-frei),
  Wurzeln = saemtliche Nullstellen der Saekularfunktion F_N (inkl. komplexer).
  F odd & Kamm symmetrisch  =>  P even (nur gerade Koeffizienten; wird geprueft).

Erzwungene Truncation (KEIN Fit-Parameter):
  K := #{gamma_k : gamma_k <= p_N}   (durch das Pol-Gitter eindeutig bestimmt).

Gemessen:
  (0) Grad-/Paritaets-Lage von P (Fernfeld-Ordnung des Quotienten Q = P/(Pi*Xi_K)).
  (1) ALLE Wurzeln von P (mp.polyroots) + Klassifikation:
        gamma-Treffer (nahe ±gamma_k), Kamm-Schatten (nahe p_j), komplex, sonstig.
  (2) Schatten-Verschiebungen delta_j = w_schatten - p_j (endlich-dim. spektrale Shift-Folge).
  (3) Newton-Bilanzen c_m = S_m(Wurzeln P) - S_m(Pole) - S_m(±gamma_{k<=K}), m = 1..6.
        odd m: strukturell 0 (Paritaet); m=2,4,6: die nicht-trivialen Momente.

SCOPE: Numerische Identitaets-Diagnose. KEIN Beweis, KEINE MS1/MS2-Schliessung, KEIN RH-Claim.

Aufruf (Laptop, < 15 min):
  PYTHONIOENCODING=utf-8 python ms1_newton_power_sum_balance.py --c 14 --N 16 --T 200 --dps 64

Autor: L.G., 2026-06-10.
"""
from __future__ import annotations
import sys, os, json, time, argparse

import mpmath as mp

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from ms1_det_zero_counting_harness import (  # noqa: E402
    _import_ccv, ground_state_even, make_poles, F_secular, GAMMA_KNOWN,
    real_zeros_of_F, RESULTS_DIR,
)


def poly_from_roots(roots):
    """Monisches Polynom-Koeffizientenliste (absteigend) aus Wurzelliste."""
    coeffs = [mp.mpf(1)]
    for r in roots:
        new = [mp.mpf(0)] * (len(coeffs) + 1)
        for i, c in enumerate(coeffs):
            new[i] += c
            new[i + 1] -= c * r
        coeffs = new
    return coeffs


def synth_div(coeffs, r):
    """Synthetische Division: coeffs / (z - r). Gibt Quotienten-Koeffizienten zurueck."""
    out = [coeffs[0]]
    for c in coeffs[1:-1]:
        out.append(c + out[-1] * r)
    # Rest = coeffs[-1] + out[-1]*r  (muss ~0 sein wenn r Wurzel ist)
    return out


def build_P(xi_full, poles):
    """P(z) = sum_j xi_j * prod_{i != j} (z - p_i) als Koeffizientenliste (absteigend)."""
    Pi = poly_from_roots(poles)            # Grad 2N+1, monisch
    n = len(Pi) - 1                        # 2N+1
    acc = [mp.mpf(0)] * n                  # Grad 2N -> n Koeffizienten
    for j, p in enumerate(poles):
        q = synth_div(Pi, p)               # Grad 2N
        w = xi_full[j]
        for i in range(n):
            acc[i] += w * q[i]
    return acc


def power_sums_from_roots(roots, m_max):
    return [sum(r ** m for r in roots) for m in range(1, m_max + 1)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--c", type=float, default=14.0)
    ap.add_argument("--N", type=int, default=16)
    ap.add_argument("--T", type=int, default=200)
    ap.add_argument("--dps", type=int, default=64)
    ap.add_argument("--m-max", type=int, default=6)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    build = _import_ccv()
    mp.mp.dps = args.dps
    t0 = time.time()
    print(f"=== Newton-Power-Sum-Bilanz: c={args.c}, N={args.N}, T={args.T}, dps={args.dps} ===")
    Q = build(c=args.c, N=args.N, T=args.T, dps=args.dps)
    print(f"  Galerkin-Matrix gebaut in {time.time()-t0:.1f}s")
    lam1, xi_full = ground_state_even(Q)
    print(f"  even-Grundzustand: lambda_1 = {mp.nstr(lam1, 8)}")

    poles = make_poles(args.N, args.c, factor=2)   # CCM-Konvention (Selbstcheck im Harness)
    p_N = max(poles)
    # gamma-Liste dynamisch via mp.zetazero erweitern, bis sie 1.5*p_N abdeckt
    # (fuer Klassifikation der Tail-Wurzeln und erzwungenes K bei grossem N)
    gammas = list(GAMMA_KNOWN)
    k_next = len(gammas) + 1
    while gammas[-1] < mp.mpf("1.5") * p_N and k_next <= 200:
        gammas.append(mp.zetazero(k_next).imag)
        k_next += 1
    K = sum(1 for g in gammas if g <= p_N)
    print(f"  p_N = {float(p_N):.4f}  =>  ERZWUNGENES K = #{{gamma_k <= p_N}} = {K} "
          f"(gamma-Liste: {len(gammas)} Eintraege bis {float(gammas[-1]):.2f})")

    # --- P(z) bauen und Struktur pruefen --------------------------------
    P = build_P(xi_full, poles)
    degP = len(P) - 1
    lead = P[0]
    sum_xi = sum(xi_full[j] for j in range(len(poles)))
    # Paritaets-Check: odd-Koeffizienten von P relativ zum fuehrenden
    odd_rel = max(abs(P[i]) for i in range(1, len(P), 2)) / abs(lead)
    print(f"  deg P = {degP} (erwartet 2N = {2*args.N});  lead = {mp.nstr(lead, 6)} "
          f"(= Sum xi_j: {mp.nstr(sum_xi, 6)});  odd-Koeff./lead = {mp.nstr(odd_rel, 3)}")

    # --- alle Wurzeln ----------------------------------------------------
    t1 = time.time()
    monic = [c / lead for c in P]
    roots = mp.polyroots(monic, maxsteps=200, extraprec=4 * args.dps)
    print(f"  polyroots: {len(roots)} Wurzeln in {time.time()-t1:.1f}s")

    # --- Klassifikation --------------------------------------------------
    comb_gap = float(2 * mp.pi / mp.log(args.c))
    tol_im = mp.mpf(10) ** (-args.dps // 4)
    cls = {"gamma_hit": [], "comb_shadow": [], "complex": [], "other": []}
    shadow_shifts = []
    for w in roots:
        if abs(w.imag) > tol_im * max(1, abs(w.real)):
            cls["complex"].append(w)
            continue
        x = w.real
        dg = min((abs(abs(x) - g), k + 1) for k, g in enumerate(gammas))
        dp = min((abs(x - p), j) for j, p in enumerate(poles))
        if dg[0] < mp.mpf("0.5"):
            cls["gamma_hit"].append((x, dg[1], dg[0]))
        elif dp[0] < mp.mpf("0.35") * comb_gap:
            cls["comb_shadow"].append((x, float(poles[dp[1]]), dp[0]))
            shadow_shifts.append((float(poles[dp[1]]), float(x - poles[dp[1]])))
        else:
            cls["other"].append(x)
    print(f"  Klassifikation: {len(cls['gamma_hit'])} gamma-Treffer, "
          f"{len(cls['comb_shadow'])} Kamm-Schatten, {len(cls['complex'])} komplex, "
          f"{len(cls['other'])} sonstige")
    for x, k, d in sorted(cls["gamma_hit"], key=lambda t: abs(t[0])):
        print(f"    gamma-Treffer: w = {mp.nstr(x, 12)}  ~ ±gamma_{k}  |diff| = {mp.nstr(d, 3)}")
    for x, p, d in sorted(cls["comb_shadow"], key=lambda t: t[0])[:50]:
        print(f"    Kamm-Schatten: w = {mp.nstr(x, 10)}  Pol p = {p:.4f}  shift = {mp.nstr(mp.mpf(x) - mp.mpf(p), 4)}")
    for x in cls["other"]:
        print(f"    SONSTIGE Wurzel: w = {mp.nstr(x, 12)}")
    if cls["complex"]:
        for w in cls["complex"][:10]:
            print(f"    komplexe Wurzel: {mp.nstr(w, 10)}")

    # --- Newton-Bilanzen --------------------------------------------------
    m_max = args.m_max
    S_W = power_sums_from_roots(roots, m_max)
    S_p = power_sums_from_roots([mp.mpc(p) for p in poles], m_max)
    gammas_K = [mp.mpc(g) for g in gammas[:K]] + [mp.mpc(-g) for g in gammas[:K]]
    S_g = power_sums_from_roots(gammas_K, m_max)
    print(f"\n  Fernfeld-Ordnung Q = P/(Pi*Xi_K):  degP - (2N+1) - 2K = "
          f"{degP - (2*args.N + 1) - 2*K}")
    print(f"  {'m':>3} {'S_m(P-Wurzeln)':>24} {'S_m(Pole)':>24} {'S_m(±gamma<=K)':>24} "
          f"{'c_m = W - p - g':>24} {'c_m / S_m(g)':>14}")
    balances = []
    for m in range(1, m_max + 1):
        c_m = S_W[m - 1] - S_p[m - 1] - S_g[m - 1]
        rel = abs(c_m) / abs(S_g[m - 1]) if abs(S_g[m - 1]) > 0 else mp.inf
        balances.append({"m": m, "S_W": mp.nstr(S_W[m-1], 12), "S_p": mp.nstr(S_p[m-1], 12),
                         "S_g": mp.nstr(S_g[m-1], 12), "c_m": mp.nstr(c_m, 12),
                         "c_m_rel": float(rel) if rel != mp.inf else None})
        print(f"  {m:>3} {mp.nstr(S_W[m-1], 10):>24} {mp.nstr(S_p[m-1], 10):>24} "
              f"{mp.nstr(S_g[m-1], 10):>24} {mp.nstr(c_m, 10):>24} {mp.nstr(rel, 4):>14}")

    # Zerlegte Bilanz: Beitrag der Klassen zu c_2
    if m_max >= 2:
        s2_hit = sum(mp.mpf(x) ** 2 for x, _, _ in cls["gamma_hit"])
        s2_shadow = sum(mp.mpf(x) ** 2 for x, _, _ in cls["comb_shadow"])
        s2_shadow_poles = sum(mp.mpf(p) ** 2 for _, (x, p, d) in
                              zip(range(len(cls["comb_shadow"])), cls["comb_shadow"]))
        s2_g = S_g[1]
        print(f"\n  c_2-Zerlegung: S2(gamma-Treffer) - S2(±gamma<=K) = "
              f"{mp.nstr(s2_hit - s2_g, 6)}  (Treffer-Defekt)")
        print(f"                 S2(Kamm-Schatten) - S2(zugeh. Pole) = "
              f"{mp.nstr(s2_shadow - s2_shadow_poles, 6)}  (Schatten-Defekt)")
        unmatched_poles = [p for p in poles
                           if not any(abs(mp.mpf(ps) - p) < 1e-9 for ps, _ in shadow_shifts)]
        s2_unmatched = sum(mp.mpf(p) ** 2 for p in unmatched_poles)
        print(f"                 S2(Pole OHNE Schatten) = {mp.nstr(s2_unmatched, 6)} "
              f"({len(unmatched_poles)} Pole)")

    out = args.out or os.path.join(
        RESULTS_DIR, f"MS1_NEWTON_BALANCE_c{int(args.c)}_N{args.N}_{time.strftime('%Y-%m-%d')}.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump({
            "script": os.path.basename(__file__),
            "purpose": "Newton-Power-Sum-Bilanz der grad-gepaarten MS1-Normalisierung (Prüfschritt)",
            "params": {"c": args.c, "N": args.N, "T": args.T, "dps": args.dps},
            "p_N": float(p_N), "K_forced": K,
            "degP": degP, "odd_coeff_rel": float(odd_rel),
            "lambda_1_even": mp.nstr(lam1, 10),
            "classification": {
                "gamma_hits": [{"w": mp.nstr(x, 20), "k": k, "absdiff": mp.nstr(d, 6)}
                               for x, k, d in cls["gamma_hit"]],
                "comb_shadows": [{"w": mp.nstr(mp.mpf(x), 16), "pole": p, "shift": s}
                                 for (x, p, d), (_, s) in zip(cls["comb_shadow"], shadow_shifts)],
                "complex": [mp.nstr(w, 12) for w in cls["complex"]],
                "other": [mp.nstr(x, 12) for x in cls["other"]],
            },
            "balances": balances,
            "runtime_s": round(time.time() - t0, 1),
        }, f, indent=2)
    print(f"\nRohdaten: {out}")


if __name__ == "__main__":
    main()
