#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MS1 — Cartwright-Ledger (Proof-Note MS1_CARTWRIGHT_LEDGER_2026-06-11.md).

Exakte Nahfeld-Buchhaltung der Paket-Paarung auf der {γ₁,γ₂}-Kontur:

    Q(z) = P(z) / [Ξ(z)·W(z)],   log Q = const + R_hit + Tail − γRest − log W

mit P aus den gespeicherten Newton-Wurzeln (hits + other), Ξ via mpmath (exakt),
W = endlicher Kamm (j≤N, ∉I) oder Sinc-Vollkamm sin(z·log c/2)/(z·log c/2)/∏_{j∈I}.

Misst margin/C-spread direkt UND rekonstruiert sie aus der Block-Zerlegung
(V1-Gate); berichtet den log-Spread-Anteil jedes Blocks (V2) und den N-Trend (V3).
Nur vorhandene Daten — keine Galerkin-Läufe.

Aufruf:
  PYTHONIOENCODING=utf-8 python ms1_cartwright_ledger.py \
      --inputs results/ms1_campaign_2026-06/MS1_NEWTON_BALANCE_c14_N40_T800_2026-06-10.json ...

Autor: L.G., 2026-06-11.
"""
from __future__ import annotations
import sys, os, json, math, argparse

import mpmath as mp

mp.mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "..", "results", "ms1_campaign_2026-06")
CACHE = os.path.join(RESULTS_DIR, "gamma_zeros_cache.json")
M_PTS = 48
R_CONTOUR = mp.mpf("4.444")


def gammas_cached():
    with open(CACHE, "r", encoding="utf-8") as f:
        return [mp.mpf(g) for g in json.load(f)]


def xi(z):
    s = mp.mpf("0.5") + mp.mpc(0, 1) * z
    return mp.mpf("0.5") * s * (s - 1) * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def sigma_rest_integral(G, m):
    """Σ_{γ>G} γ^(−2m) via Riemann–von-Mangoldt-Dichte (1/2π)·log(t/2π)."""
    G = float(G)
    return (1.0 / (2 * math.pi)) * G ** (1 - 2 * m) / (2 * m - 1) * (
        math.log(G / (2 * math.pi)) + 1.0 / (2 * m - 1))


def log_gamma_rest_tail(z, G):
    """−Σ_{γ>G} log(1−z²/γ²) ≈ Σ_m z^(2m)/m · σ_m^rest (Integral-Dichte), m=1..8."""
    acc = mp.mpc(0)
    z2 = z * z
    for m in range(1, 9):
        acc += z2 ** m / m * mp.mpf(sigma_rest_integral(G, m))
    return acc


def contour(z0, r, n=M_PTS):
    return [z0 + r * mp.e ** (mp.mpc(0, 1) * 2 * mp.pi * k / n) for k in range(n)]


def margin_from_logq(logqs):
    """C = geometrisches Mittel; margin = max|Q/C − 1|; spread = max|Q|/min|Q|."""
    # Imaginärteile (Phasen) entlang der Kontur entwickeln (winding 0 ⟹ schließt)
    re = [mp.re(l) for l in logqs]
    im = [mp.im(l) for l in logqs]
    unw = [im[0]]
    for v in im[1:]:
        d = v - unw[-1]
        while d > mp.pi: d -= 2 * mp.pi
        while d < -mp.pi: d += 2 * mp.pi
        unw.append(unw[-1] + d)
    mre = sum(re) / len(re); mim = sum(unw) / len(unw)
    margin = mp.mpf(0)
    for r_, i_ in zip(re, unw):
        margin = max(margin, abs(mp.e ** mp.mpc(r_ - mre, i_ - mim) - 1))
    spread = mp.e ** (max(re) - min(re))
    return float(margin), float(spread), float(max(re) - min(re))


def analyze(path, gammas, mode, window=None):
    """mode: 'finite' (Nenner-Kamm = alle j≤N ∉ S), 'sinc' (zusätzlich j>N via Sinus),
    'none' (S = alle j≤N: kein Kamm-Nenner). window w (Harness-Konvention,
    dist zum Konturrand < w) erweitert S über die inneren Pole hinaus."""
    with open(path, "r", encoding="utf-8") as f:
        d = json.load(f)
    c = mp.mpf(d["params"]["c"]); N = int(d["params"]["N"])
    delta = 2 * mp.pi / mp.log(c)
    cl = d["classification"]
    hits = [(int(h["k"]), mp.mpf(h["w"])) for h in cl["gamma_hits"] if mp.mpf(h["w"]) > 0]
    tail = [mp.mpf(t) for t in sorted({abs(float(mp.mpf(w))) for w in cl["other"]})]
    hit_ks = {k for k, _ in hits}
    G = gammas[-1]
    z0 = (gammas[0] + gammas[1]) / 2
    pts = contour(z0, R_CONTOUR)
    inner = [j for j in range(1, N + 1) if abs(j * delta - z0) < R_CONTOUR]
    S = set(inner)
    if window is not None:
        S |= {j for j in range(1, N + 1)
              if abs(j * delta - z0) - R_CONTOUR < mp.mpf(window)}
    if mode == "none":
        S = set(range(1, N + 1))

    def blocks(z):
        b = {}
        b["hit"] = sum(mp.log(1 - z * z / (w * w)) - mp.log(1 - z * z / (gammas[k - 1] ** 2))
                       for k, w in hits)
        b["tail"] = sum(mp.log(1 - z * z / (t * t)) for t in tail)
        # grest = −Σ_{k≤kmax,∉hits} log(1−z²/γ²) − Σ_{γ>G} log(1−z²/γ²) (Integral-Approx)
        b["grest"] = (-sum(mp.log(1 - z * z / (g * g))
                           for k, g in enumerate(gammas, start=1) if k not in hit_ks)
                      + log_gamma_rest_tail(z, G))
        if mode == "sinc":
            arg = z * mp.log(c) / 2
            b["comb"] = -(mp.log(mp.sin(arg) / arg)
                          - sum(mp.log(1 - z * z / ((j * delta) ** 2)) for j in S))
        else:  # 'finite'/'none'/Fenster: Nenner-Kamm = j≤N außerhalb S
            b["comb"] = -sum(mp.log(1 - z * z / ((j * delta) ** 2))
                             for j in range(1, N + 1) if j not in S)
        return b

    # direkte Messung: Q = P/(Ξ·W)
    logq_meas, logq_pred, per_block = [], [], {k: [] for k in ("hit", "tail", "grest", "comb")}
    for z in pts:
        b = blocks(z)
        for k in per_block: per_block[k].append(mp.re(b[k]))
        logq_pred.append(b["hit"] + b["tail"] + b["grest"] + b["comb"])
        # Messung: P exakt aus Wurzeln, Ξ via mpmath, W je Modus
        logP = sum(mp.log(z * z - w * w) for _, w in hits) + \
               sum(mp.log(z * z - t * t) for t in tail)
        if mode == "sinc":
            arg = z * mp.log(c) / 2
            logW = mp.log(mp.sin(arg) / arg) - sum(
                mp.log(1 - z * z / ((j * delta) ** 2)) for j in S)
        else:
            logW = sum(mp.log(1 - z * z / ((j * delta) ** 2))
                       for j in range(1, N + 1) if j not in S)
        logq_meas.append(logP - mp.log(xi(z)) - logW)

    mar_m, spr_m, lsp_m = margin_from_logq(logq_meas)
    mar_p, spr_p, lsp_p = margin_from_logq(logq_pred)
    spreads = {k: float(max(v) - min(v)) for k, v in per_block.items()}
    return dict(source=os.path.basename(path), c=float(c), N=N, mode=mode,
                window=window, n_S=len(S), inner_poles=len(inner),
                n_hits=len(hits), n_tail=len(tail),
                margin_measured=mar_m, margin_predicted=mar_p,
                logspread_measured=lsp_m, logspread_predicted=lsp_p,
                block_logspread=spreads)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--scan", action="store_true",
                    help="Fenster-Scan w=1.5..96 + Endpunkte none/finite/sinc")
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()
    gammas = gammas_cached()
    print(f"γ-Cache: {len(gammas)} Nullstellen bis {float(gammas[-1]):.1f}")
    results = []
    if args.scan:
        jobs = [("finite", None)] + [("finite", w) for w in (1.5, 3, 6, 12, 24, 48, 96)] \
               + [("none", None), ("sinc", None)]
    else:
        jobs = [("finite", None), ("sinc", None)]
    for path in args.inputs:
        print(f"\n=== {os.path.basename(path)} ===")
        print(f"  {'Modus':>7} {'w':>5} {'#S':>4} {'margin':>9} {'pred':>9} {'logSpread':>9}"
              f"  {'R_hit':>7} {'Tail':>7} {'γRest':>7} {'Kamm':>7}")
        for mode, w in jobs:
            r = analyze(path, gammas, mode, window=w)
            results.append(r)
            bs = r["block_logspread"]
            wtxt = f"{w:g}" if w is not None else "—"
            print(f"  {mode:>7} {wtxt:>5} {r['n_S']:>4} {r['margin_measured']:>9.4g}"
                  f" {r['margin_predicted']:>9.4g} {r['logspread_measured']:>9.4f}"
                  f"  {bs['hit']:>7.4f} {bs['tail']:>7.4f} {bs['grest']:>7.4f} {bs['comb']:>7.4f}")
    out = args.out or os.path.join(RESULTS_DIR, "MS1_CARTWRIGHT_LEDGER_2026-06-11.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(dict(script="ms1_cartwright_ledger.py", results=results), f, indent=1)
    print(f"\nRohdaten: {out}")


if __name__ == "__main__":
    main()
