"""
c2e_anchor_flatness_diagnostic.py - C2E diagnostics for D2

Diagnoses the D2 anchor-flatness program:
1. Anchor: |m_lambda(w_*) - alpha_lambda|
2. Flatness: sup_{k notin Cl} |m_lambda(w_k) - m_lambda(w_*)|
3. Cross-measure structure nu_j = <e_j,f><e_j,g>
4. Direct-vs-factorized check for m_lambda(w_k)

The core scalar kernel is
    m_lambda(w) = <(T_lambda - w I)^(-1) f_lambda, g_lambda>
with
    T_lambda = P0 H P0
    f_lambda = P0 H u_tilde
    g_lambda = P0 k_full

Author: Codex
Date: 2026-04-19
"""

import json
import os
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from dirichlet_ccm_fourier_mp import (  # noqa: E402
    build_QW_mp,
    chi_trivial,
    diagonalize_mp,
    project_to_parity,
)
from c2_approximation_test import (  # noqa: E402
    inner_product,
    k_lambda_value,
    norm,
)
from c2_poisson_decomposition import project_to_fourier  # noqa: E402


DPS = int(os.environ.get("DPS", 50))
CLUSTER_DELTA = float(os.environ.get("CLUSTER_DELTA", "1e-10"))
CONFIGS_JSON = os.environ.get("CONFIGS_JSON", "")
OUTPUT_JSON = os.environ.get("OUTPUT_JSON", "")

DEFAULT_CONFIGS = [
    {"lam": 3.0, "N": 30},
    {"lam": 5.0, "N": 55},
    {"lam": 7.0, "N": 85},
]


def parse_configs():
    if CONFIGS_JSON:
        return json.loads(CONFIGS_JSON)
    return DEFAULT_CONFIGS


def exp_fit(vals, lams):
    pairs = [(l, v) for l, v in zip(lams, vals) if abs(v) > 0]
    if len(pairs) < 2:
        return float("nan")
    xs = [np.log(l) for l, _ in pairs]
    ys = [np.log(abs(v)) for _, v in pairs]
    n = len(xs)
    sx = sum(xs)
    sy = sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    den = n * sxx - sx * sx
    if den == 0:
        return float("nan")
    return (n * sxy - sx * sy) / den


def mat_vec_mul(M, v, dim):
    out = mp.matrix(dim, 1)
    for i in range(dim):
        s = mp.mpf(0)
        for j in range(dim):
            s += M[i, j] * v[j, 0]
        out[i, 0] = s
    return out


def project_qw_components(lam, N):
    dim = N + 1
    lam_mp = mp.mpf(lam)
    L_mp = 2 * mp.log(lam_mp)

    # QW (with W02) and H (without W02), both in even sector
    Mq_full, _ = build_QW_mp(N, lam, chi_trivial, 1, True, -1, verbose=False)
    Mh_full, _ = build_QW_mp(N, lam, chi_trivial, 1, False, -1, verbose=False)
    Aq, _ = project_to_parity(Mq_full, N, parity="even")
    Ah, _ = project_to_parity(Mh_full, N, parity="even")

    W02 = mp.matrix(dim, dim)
    H = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            W02[i, j] = Aq[i, j] - Ah[i, j]
            H[i, j] = Ah[i, j]

    ws, Qs = diagonalize_mp(Aq, verbose=False)
    wmin = float(ws[0])

    # Extract u_tilde from rank-one W02
    col0 = mp.matrix(dim, 1)
    for i in range(dim):
        col0[i, 0] = W02[i, 0]
    cn = norm(col0, dim)
    ut = mp.matrix(dim, 1)
    for i in range(dim):
        ut[i, 0] = col0[i, 0] / cn

    # P0 = I - |u><u|
    P0 = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            P0[i, j] = (mp.mpf(1) if i == j else mp.mpf(0)) - ut[i, 0] * ut[j, 0]

    # T = P0 H P0
    HP0 = H * P0
    T = P0 * HP0

    # f = P0 H u
    Hut = mat_vec_mul(H, ut, dim)
    f = mat_vec_mul(P0, Hut, dim)

    # k_full projected into even Fourier basis
    def kfull_of_x(x):
        u = mp.exp(x) / lam_mp
        return k_lambda_value(float(u), lam)

    cf = project_to_fourier(kfull_of_x, L_mp, N)
    nf = norm(cf, dim)
    kn = mp.matrix(dim, 1)
    for i in range(dim):
        kn[i, 0] = cf[i, 0] / nf

    alpha = float(inner_product(ut, kn, dim))
    g = mat_vec_mul(P0, kn, dim)

    # Tight cluster by eigenvalue proximity to w_min
    cl_end = 0
    for k in range(dim):
        if float(ws[k]) - wmin < CLUSTER_DELTA:
            cl_end = k

    return {
        "dim": dim,
        "L_mp": L_mp,
        "Aq": Aq,
        "H": H,
        "T": T,
        "P0": P0,
        "ut": ut,
        "f": f,
        "g": g,
        "kn": kn,
        "alpha": alpha,
        "ws": [float(w) for w in ws],
        "Qs": Qs,
        "cl_end": cl_end,
    }


def diagonalize_T(T):
    ts, Es = diagonalize_mp(T, verbose=False)
    return [float(t) for t in ts], Es


def spectral_data(ts, Es, f, g, dim):
    nus = []
    overlaps_f = []
    overlaps_g = []
    for j in range(dim):
        ej = mp.matrix(dim, 1)
        for i in range(dim):
            ej[i, 0] = Es[i, j]
        of = float(inner_product(ej, f, dim))
        og = float(inner_product(ej, g, dim))
        nu = of * og
        overlaps_f.append(of)
        overlaps_g.append(og)
        nus.append(nu)
    return overlaps_f, overlaps_g, nus


def m_from_spectral(ts, nus, w, pole_tol=1e-24):
    s = mp.mpf(0)
    for tj, nuj in zip(ts, nus):
        den = mp.mpf(tj) - mp.mpf(w)
        if abs(den) < pole_tol:
            return float("nan")
        s += mp.mpf(nuj) / den
    return float(s)


def mprime_from_spectral(ts, nus, w, pole_tol=1e-24):
    s = mp.mpf(0)
    for tj, nuj in zip(ts, nus):
        den = mp.mpf(tj) - mp.mpf(w)
        if abs(den) < pole_tol:
            return float("nan")
        s += mp.mpf(nuj) / (den * den)
    return float(s)


def weighted_center(ts, nus):
    weights = [abs(x) for x in nus]
    total = sum(weights)
    if total == 0:
        return 0.0
    return sum(w * t for w, t in zip(weights, ts)) / total


def cross_measure_stats(ts, nus, w_star):
    tv = sum(abs(x) for x in nus)
    pos = sum(abs(x) for x in nus if x > 0)
    neg = sum(abs(x) for x in nus if x < 0)
    coherence = abs(sum(nus)) / tv if tv > 0 else 0.0
    center = weighted_center(ts, nus)
    M0 = sum(nus)
    M1 = sum((t - center) * nu for t, nu in zip(ts, nus))
    M2_signed = sum(((t - center) ** 2) * nu for t, nu in zip(ts, nus))
    M2_abs = sum(((t - center) ** 2) * abs(nu) for t, nu in zip(ts, nus))
    invdist = sum(abs(nu) / max(abs(t - w_star), 1e-30) for t, nu in zip(ts, nus))
    return {
        "tv": tv,
        "pos_mass": pos,
        "neg_mass": neg,
        "sign_coherence": coherence,
        "center_abs": center,
        "M0": M0,
        "M1": M1,
        "M2_signed": M2_signed,
        "M2_abs": M2_abs,
        "invdist_mass_at_wstar": invdist,
    }


def analyze_config(lam, N):
    mp.mp.dps = DPS
    t0 = time.time()

    pieces = project_qw_components(lam, N)
    dim = pieces["dim"]
    alpha = pieces["alpha"]
    ws = pieces["ws"]
    Qs = pieces["Qs"]
    ut = pieces["ut"]
    kn = pieces["kn"]
    T = pieces["T"]
    f = pieces["f"]
    g = pieces["g"]
    cl_end = pieces["cl_end"]
    w_star = ws[cl_end + 1] if cl_end + 1 < dim else ws[cl_end]

    ts, Es = diagonalize_T(T)
    overlaps_f, overlaps_g, nus = spectral_data(ts, Es, f, g, dim)
    cross_stats = cross_measure_stats(ts, nus, w_star)

    rows = []
    for k in range(dim):
        qk = mp.matrix(dim, 1)
        for i in range(dim):
            qk[i, 0] = Qs[i, k]

        beta = float(inner_product(ut, qk, dim))
        c = float(inner_product(kn, qk, dim))
        w = ws[k]
        is_cl = k <= cl_end

        if abs(beta) > 1e-18:
            delta_fact = c / beta
            m_fact = alpha - delta_fact
        else:
            delta_fact = float("nan")
            m_fact = float("nan")

        m_direct = m_from_spectral(ts, nus, w)
        mprime = mprime_from_spectral(ts, nus, w)
        delta_direct = alpha - m_direct
        formula_err = abs(m_fact - m_direct) if not np.isnan(m_fact) else float("nan")

        row = {
            "k": k,
            "w": w,
            "beta": beta,
            "c": c,
            "m_fact": m_fact,
            "m_direct": m_direct,
            "delta_fact": delta_fact,
            "delta_direct": delta_direct,
            "formula_err": formula_err,
            "mprime": mprime,
            "is_cl": is_cl,
        }
        rows.append(row)

    off_rows = [r for r in rows if (not r["is_cl"]) and (not np.isnan(r["m_direct"]))]
    anchor_abs = 0.0
    flat_sup = 0.0
    flat_rms = 0.0
    dprime_sup = 0.0
    if off_rows:
        m_star = off_rows[0]["m_direct"]
        anchor_abs = abs(m_star - alpha)
        flat_vals = [abs(r["m_direct"] - m_star) for r in off_rows]
        flat_sup = max(flat_vals)
        flat_rms = float(np.sqrt(np.mean([x * x for x in flat_vals])))
        finite_mprimes = [abs(r["mprime"]) for r in off_rows if not np.isnan(r["mprime"])]
        dprime_sup = max(finite_mprimes) if finite_mprimes else 0.0

    delta_sup = max(abs(r["delta_direct"]) for r in off_rows) if off_rows else 0.0
    delta_rms = float(np.sqrt(np.mean([r["delta_direct"] ** 2 for r in off_rows]))) if off_rows else 0.0
    beta2_sum = sum(r["beta"] ** 2 for r in off_rows)
    c2_sum = sum(r["c"] ** 2 for r in off_rows)
    eff_delta2 = c2_sum / beta2_sum if beta2_sum > 0 else 0.0

    # finite-difference flatness on off-cluster band
    fdiffs = []
    for i in range(len(off_rows) - 1):
        dw = off_rows[i + 1]["w"] - off_rows[i]["w"]
        if abs(dw) > 1e-30:
            fd = (off_rows[i + 1]["m_direct"] - off_rows[i]["m_direct"]) / dw
            fdiffs.append(fd)
    fdiff_sup = max(abs(x) for x in fdiffs) if fdiffs else 0.0

    return {
        "lam": lam,
        "N": N,
        "alpha": alpha,
        "dim": dim,
        "cl_size": cl_end + 1,
        "w_star": w_star,
        "anchor_abs": anchor_abs or 0.0,
        "flat_sup": flat_sup,
        "flat_rms": flat_rms,
        "delta_sup": delta_sup,
        "delta_rms": delta_rms,
        "beta2_sum": beta2_sum,
        "c2_sum": c2_sum,
        "eff_delta2": eff_delta2,
        "mprime_sup": dprime_sup,
        "fdiff_sup": fdiff_sup,
        "cross_stats": cross_stats,
        "rows": rows,
        "off_rows": off_rows,
        "ts": ts,
        "nus": nus,
        "overlaps_f": overlaps_f,
        "overlaps_g": overlaps_g,
        "time_s": time.time() - t0,
    }


def print_config_report(r):
    print(f"\n{'=' * 88}")
    print(f"  C2E: lambda={r['lam']}, N={r['N']}, cluster={r['cl_size']}/{r['dim']}")
    print(f"{'=' * 88}")
    print(f"  alpha               = {r['alpha']:.8f}")
    print(f"  w_star              = {r['w_star']:.8f}")
    print(f"  anchor |m(w*)-a|    = {r['anchor_abs']:.6e}")
    print(f"  flat sup            = {r['flat_sup']:.6e}")
    print(f"  flat RMS            = {r['flat_rms']:.6e}")
    print(f"  sup |delta_k|       = {r['delta_sup']:.6e}")
    print(f"  RMS |delta_k|       = {r['delta_rms']:.6e}")
    print(f"  eff delta^2         = {r['eff_delta2']:.6e}")
    print(f"  sum |beta|^2        = {r['beta2_sum']:.6e}")
    print(f"  sum |c|^2           = {r['c2_sum']:.6e}")
    print(f"  sup |m'(w_k)|       = {r['mprime_sup']:.6e}")
    print(f"  sup |FD slope|      = {r['fdiff_sup']:.6e}")
    print(f"  time                = {r['time_s']:.1f}s")

    cs = r["cross_stats"]
    print("\n  Kreuzmass:")
    print(f"    total variation   = {cs['tv']:.6e}")
    print(f"    positive mass     = {cs['pos_mass']:.6e}")
    print(f"    negative mass     = {cs['neg_mass']:.6e}")
    print(f"    sign coherence    = {cs['sign_coherence']:.6e}")
    print(f"    center_abs        = {cs['center_abs']:.8f}")
    print(f"    M0                = {cs['M0']:.6e}")
    print(f"    M1                = {cs['M1']:.6e}")
    print(f"    M2_signed         = {cs['M2_signed']:.6e}")
    print(f"    M2_abs            = {cs['M2_abs']:.6e}")
    print(f"    invdist@w_star    = {cs['invdist_mass_at_wstar']:.6e}")

    print("\n  Off-cluster profile:")
    print(f"  {'k':>4s}  {'m_dir':>12s}  {'|m-a|':>12s}  {'|m-m*|':>12s}  {'m_fact err':>12s}  {'|mprime|':>12s}")
    m_star = r["off_rows"][0]["m_direct"] if r["off_rows"] else float("nan")
    for row in r["off_rows"][:8]:
        print(
            f"  {row['k']:4d}  {row['m_direct']:12.8f}  "
            f"{abs(row['m_direct'] - r['alpha']):12.4e}  "
            f"{abs(row['m_direct'] - m_star):12.4e}  "
            f"{row['formula_err']:12.4e}  {abs(row['mprime']):12.4e}"
        )

    print("\n  Largest |nu_j| modes:")
    ranked = sorted(
        enumerate(zip(r["ts"], r["nus"], r["overlaps_f"], r["overlaps_g"])),
        key=lambda x: abs(x[1][1]),
        reverse=True,
    )
    print(f"  {'j':>4s}  {'t_j':>12s}  {'nu_j':>12s}  {'<e,f>':>12s}  {'<e,g>':>12s}")
    for j, (tj, nuj, of, og) in ranked[:8]:
        print(f"  {j:4d}  {tj:12.8f}  {nuj:12.4e}  {of:12.4e}  {og:12.4e}")


def print_scaling_report(results):
    lams = [r["lam"] for r in results]
    metrics = [
        ("anchor^2", [r["anchor_abs"] ** 2 for r in results]),
        ("flat_sup^2", [r["flat_sup"] ** 2 for r in results]),
        ("flat_rms^2", [r["flat_rms"] ** 2 for r in results]),
        ("delta_sup^2", [r["delta_sup"] ** 2 for r in results]),
        ("delta_rms^2", [r["delta_rms"] ** 2 for r in results]),
        ("eff_delta^2", [r["eff_delta2"] for r in results]),
        ("sum|beta|^2", [r["beta2_sum"] for r in results]),
        ("sum|c|^2", [r["c2_sum"] for r in results]),
        ("mprime_sup", [r["mprime_sup"] for r in results]),
        ("fdiff_sup", [r["fdiff_sup"] for r in results]),
        ("1-sign_coh", [1.0 - r["cross_stats"]["sign_coherence"] for r in results]),
    ]

    print(f"\n{'=' * 88}")
    print("  C2E MULTI-LAMBDA SCALING")
    print(f"{'=' * 88}")
    header = f"  {'metric':<16s}"
    for r in results:
        header += f"  {'lam=' + str(r['lam']):>12s}"
    header += f"  {'exp(fit)':>10s}"
    print(header)
    for name, vals in metrics:
        line = f"  {name:<16s}"
        for v in vals:
            line += f"  {v:12.6e}"
        line += f"  {exp_fit(vals, lams):+10.4f}"
        print(line)


def maybe_write_json(results):
    if not OUTPUT_JSON:
        return
    payload = {
        "cluster_delta": CLUSTER_DELTA,
        "results": results,
    }
    with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    print(f"\nSaved JSON to: {OUTPUT_JSON}")


def main():
    mp.mp.dps = DPS
    configs = parse_configs()
    results = []
    for cfg in configs:
        r = analyze_config(cfg["lam"], cfg["N"])
        results.append(r)
        print_config_report(r)

    if len(results) >= 2:
        print_scaling_report(results)

    maybe_write_json(results)


if __name__ == "__main__":
    main()
