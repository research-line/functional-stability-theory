#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MS1 — Determinanten-Zero-Counting-Harness (Programmschritt 3, 2026-06-03).

Zweck (Zookeeper, Z-Route / MS1-Reformulierung):
  Empirische Erst-Kartierung des Rouché-Kriteriums aus dem Analyseprogramm
  vom 2026-06-03 (interne Notiz, nicht veröffentlicht; Punkt 3):

      |d_lambda(z) - C_lambda * Xi(z)| < |C_lambda * Xi(z)|   auf Konturen um
      Nullstellen-PAKETE (nicht einzelne Nullstellen).

  Als Cluster-Determinante dient die CCM-Eq-5.25-Säkularfunktion des
  Rang-1-gestörten Skalierungsoperators (vgl. _scripts/dirichlet_ccm_operator.py):

      F_N(z) = sum_{j=-N}^{N} xi_j / (z - p_j),   p_j = Basis-Frequenz j,

  deren Nullstellen die Eigenwerte von D_log^(lambda,N) sind und die niedrigen
  Riemann-Nullstellen gamma_k approximieren. Windungs-konsistente Rouché-Paarung:

      H_S(z) = F_N(z) * prod_{p in S} (z - p)    vs.   C * Xi(z),

  wobei S die Pole im Konturinneren enthaelt (Pflicht: Pol-Kompensation, macht
  H holomorph mit winding(H) = #F-Nullstellen innen = winding(Xi)) sowie optional
  kontur-nahe AEUSSERE Pole (Konditionierung: multipliziert den dominanten
  Bulk-Pol-Gradienten heraus, fuegt innen KEINE Nullstelle hinzu). Der restliche
  O(1)-Bulk wird in die Konstante C ausfaktoriert. Das ist die lokale Version
  von "det_reg(I - z K_lambda) nach Ausfaktorieren des O(1)-Bulks".

GEMESSEN wird pro Kontur:
  - margin  = max_{z in K} |F/(C*Xi) - 1|  (bzw. G-Version)  -> Rouché zertifiziert
              gleiche Nullstellenzahl im Inneren gdw. margin < 1.
  - C-Streuung über die Kontur (max/min von |F/Xi|) -> wie konstant ist die
    Bulk-Normalisierung wirklich (Diagnose für die Cartwright-Limes-Hypothese).
  - winding number von F bzw. G (Phasen-Tracking) als Argument-Prinzip-Cross-Check.
  - Nullstellen-Fehler |z_k - gamma_k| (Brent auf F entlang der reellen Achse).
  - N-Sweep: Trend margin(N) -- schrumpft der Rouché-Defekt mit der Trunkierung?

SCOPE-GRENZE (kein Overclaim):
  - Rein numerische Erst-Kartierung des Pivots. KEIN Beweis von Schritt (2)
    (Trace-Class-Konvergenz), KEINE MS1-Schließung, KEINE RH-Aussage.
  - Die analytische Fortsetzungs-Wand (C+-Pol-Kamm, Programm-Note
    interne Programmnotiz vom 2026-06-03, §2) wird hiervon
    NICHT berührt: die Matrix ist zeta-seitig gebaut (CvS-Weil-Form), der Test
    quantifiziert nur, ob die TRUNKIERTE Determinante Zero-Counting im Sinne
    des reformulierten MS1 trägt und mit welchem Defekt.

Selbst-Check Pol-Konvention: Groskins Basis e_k(t)=exp(2*pi*i*k*t/(2*log c))
  hat Frequenzen pi*k/log(c); die CCM-Variante (dirichlet_ccm_operator.py)
  nutzt 2*pi*k/L mit L=2*log(lambda)=log(c) -- numerisch identisch zu 2*pi*k/log(c).
  Faktor-2-Ambiguität! Der Harness testet BEIDE Pol-Gitter und wählt das, dessen
  F-Nullstellen die gamma_k treffen (Report dokumentiert die Wahl).

Aufruf (Laptop-tauglich, < 1 h):
  PYTHONIOENCODING=utf-8 python ms1_det_zero_counting_harness.py \
      --c 14 --N-list 16,24,32 --T 200 --dps 60

Autor: L.G., 2026-06-09.
Refs: CCM arXiv:2511.22755 Eq 5.25; CvS arXiv:2511.23257 Prop 4.1;
      Groskin arXiv:2605.20224 (connes-cvs, https://github.com/akivag613/connes-cvs-); Prüfung 2026-06-03.
"""
from __future__ import annotations
import sys, os, json, time, argparse

import mpmath as mp

# CCV-Quelle: Umgebungsvariable CCV_SRC = Pfad zum geklonten connes-cvs-Repo (MIT)
_CCV_CANDIDATES = [
    os.environ.get("CCV_SRC", ""),
]

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "ms1_campaign_2026-06")

GAMMA_KNOWN = [
    mp.mpf("14.134725141734693790"),
    mp.mpf("21.022039638771554993"),
    mp.mpf("25.010857580145688763"),
    mp.mpf("30.424876125859513210"),
    mp.mpf("32.935061587739189691"),
    mp.mpf("37.586178158825671257"),
    mp.mpf("40.918719012147495187"),
    mp.mpf("43.327073280914999519"),
    mp.mpf("48.005150881167159727"),
    mp.mpf("49.773832477672302181"),
]


def _import_ccv():
    for cand in _CCV_CANDIDATES:
        if cand and os.path.isdir(cand):
            if cand not in sys.path:
                sys.path.insert(0, cand)
            try:
                from connes_cvs import build_galerkin_matrix  # noqa
                return build_galerkin_matrix
            except ImportError:
                continue
    raise ImportError("connes_cvs nicht gefunden — CCV_SRC setzen")


# ----------------------------------------------------------------------------
# even-Sektor (wie z1_aib_pollution_audit.py)
# ----------------------------------------------------------------------------

def even_sector_matrix_and_map(Q: "mp.matrix"):
    """Projiziere (2N+1)x(2N+1) auf even-Sektor; gib (A_even, V) zurueck,
    V = (2N+1)x(N+1)-Isometrie zum Zurueckmappen der Eigenvektoren."""
    DIM = Q.rows
    N = (DIM - 1) // 2
    V = mp.matrix(DIM, N + 1)
    V[N, 0] = mp.mpf(1)
    inv = 1 / mp.sqrt(2)
    for k in range(1, N + 1):
        V[N + k, k] = inv
        V[N - k, k] = inv
    return V.T * Q * V, V


def ground_state_even(Q: "mp.matrix"):
    """Unterster even-Eigenwert + Eigenvektor (full-Basis-Koordinaten xi_j, j=-N..N)."""
    A, V = even_sector_matrix_and_map(Q)
    E, EV = mp.eigsy(A)
    pairs = sorted((mp.mpf(E[i]), i) for i in range(A.rows))
    lam1, i1 = pairs[0]
    v = mp.matrix([EV[r, i1] for r in range(A.rows)])
    xi_full = V * v
    return lam1, xi_full


# ----------------------------------------------------------------------------
# Säkularfunktion / -polynom und Xi
# ----------------------------------------------------------------------------

def make_poles(N: int, c, factor: int):
    """Pol-Gitter p_j = factor*pi*j/log(c), j=-N..N (factor=1: Groskin-Frequenzen,
    factor=2: CCM-2pi/L-Konvention)."""
    logc = mp.log(c)
    return [factor * mp.pi * j / logc for j in range(-N, N + 1)]


def F_secular(z, xi_full, poles):
    total = mp.mpc(0)
    for j, p in enumerate(poles):
        total += xi_full[j] / (z - p)
    return total


def H_adjusted(z, xi_full, poles, S):
    """H_S(z) = F(z) * prod_{p in S} (z - p). Fuer S = Pole im Konturinneren
    (+ optional nahe aeussere Pole) ist H im Inneren holomorph mit
    winding(H) = Anzahl der F-Nullstellen im Inneren."""
    prod = mp.mpc(1)
    for p in S:
        prod *= (z - p)
    return F_secular(z, xi_full, poles) * prod


def Xi_completed(z):
    """Xi(z) = xi(1/2 + i z), xi(s) = s(s-1)/2 * pi^(-s/2) * Gamma(s/2) * zeta(s)."""
    s = mp.mpf("0.5") + mp.mpc(0, 1) * z
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


# ----------------------------------------------------------------------------
# Nullstellen von F auf der reellen Achse (Pol-Konventions-Selbstcheck)
# ----------------------------------------------------------------------------

def real_zeros_of_F(xi_full, poles, lo, hi, samples=4000):
    """Vorzeichenwechsel-Scan + Bisektion fuer F auf (lo,hi), Pole ausgespart."""
    pts = [lo + (hi - lo) * k / samples for k in range(samples + 1)]
    pset = sorted(p for p in poles if lo < p < hi)
    vals = []
    for x in pts:
        if any(abs(x - p) < mp.mpf("1e-6") for p in pset):
            vals.append(None)
            continue
        vals.append(F_secular(mp.mpf(x), xi_full, poles).real)
    zeros = []
    for i in range(len(pts) - 1):
        a, b = vals[i], vals[i + 1]
        if a is None or b is None or a * b >= 0:
            continue
        if any(pts[i] < p < pts[i + 1] for p in pset):
            continue
        f = lambda x: F_secular(mp.mpf(x), xi_full, poles).real
        try:
            zeros.append(mp.findroot(f, (mp.mpf(pts[i]), mp.mpf(pts[i + 1])), solver="bisect", tol=mp.mpf("1e-30")))
        except Exception:
            pass
    return zeros


def pick_pole_convention(xi_full, N, c, verbose=True):
    """Teste factor=1 und factor=2; waehle das Gitter, dessen F-Nullstellen
    gamma_1 am besten treffen. Scan bis min(zmax, 55) — liefert auch die
    aeusseren Bulk-Wurzeln fuer das B(z)-Modell (Folgeschritt 2)."""
    best = None
    for factor in (1, 2):
        poles = make_poles(N, c, factor)
        zmax = float(max(poles))
        hi = min(zmax, 55.0)
        zs = real_zeros_of_F(xi_full, poles, mp.mpf("0.5"), mp.mpf(hi), samples=4500)
        if not zs:
            err = mp.inf
        else:
            err = min(abs(z - GAMMA_KNOWN[0]) for z in zs)
        if verbose:
            print(f"  [pole-check] factor={factor}: zeros in (0.5,{hi:.1f}): "
                  f"{len(zs)}, best |z-gamma_1| = {mp.nstr(err, 6)}")
        if best is None or err < best[1]:
            best = (factor, err, zs)
    return best


# ----------------------------------------------------------------------------
# Rouché-Kontur-Messung
# ----------------------------------------------------------------------------

def contour_circle(center, radius, M):
    return [center + radius * mp.exp(mp.mpc(0, 2) * mp.pi * k / M) for k in range(M)]


def winding_number(vals):
    """Numerische Windungszahl aus Phasen-Tracking einer geschlossenen Werteliste."""
    total = mp.mpf(0)
    for i in range(len(vals)):
        a = vals[i]
        b = vals[(i + 1) % len(vals)]
        dphi = mp.arg(b / a)
        total += dphi
    return total / (2 * mp.pi)


def rouche_on_contour(xi_full, poles, center, radius, M=48, pole_adjust=None,
                      bulk_model=None):
    """Messe Rouché-Margin auf Kreis um center.
    pole_adjust=None: vergleiche F(z) gegen C*Xi(z)   (Kontur muss polfrei sein)
    pole_adjust=S:    vergleiche H_S(z)=F*prod_{p in S}(z-p) gegen C*Xi(z);
                      S MUSS alle inneren Pole enthalten, darf nahe aeussere enthalten.
    bulk_model=(w_out, g_out): zusaetzlicher expliziter Bulk-Faktor
                      B(z) = prod_k (z^2-w_k^2)/(z^2-gamma_k^2) auf der Xi-Seite
                      (aeussere F-Wurzeln gegen aeussere Zeta-Nullstellen gepaart;
                      alle ausserhalb der Kontur -> windungsneutral).
    C := geometrisches Mittel des Quotienten ueber die Kontur."""
    K = contour_circle(center, radius, M)
    num_vals, den_vals = [], []
    for z in K:
        if pole_adjust is not None:
            num = H_adjusted(z, xi_full, poles, pole_adjust)
        else:
            num = F_secular(z, xi_full, poles)
        den = Xi_completed(z)
        if bulk_model is not None:
            w_out, g_out = bulk_model
            B = mp.mpc(1)
            for w, g in zip(w_out, g_out):
                B *= (z * z - w * w) / (z * z - g * g)
            den = den * B
        num_vals.append(num)
        den_vals.append(den)
    # geometrisches Mittel des Quotienten (log-Mittel), nur Betrag + mittlere Phase
    logs = [mp.log(num_vals[i] / den_vals[i]) for i in range(M)]
    mean_log = sum(l.real for l in logs) / M
    # Phase: zirkular mitteln
    mean_phase = mp.atan2(sum(mp.sin(l.imag) for l in logs) / M,
                          sum(mp.cos(l.imag) for l in logs) / M)
    C = mp.exp(mp.mpc(mean_log, mean_phase))
    ratios = [num_vals[i] / (C * den_vals[i]) for i in range(M)]
    margin = max(abs(r - 1) for r in ratios)
    spread = max(abs(r) for r in ratios) / min(abs(r) for r in ratios)
    wind_num = winding_number(num_vals)
    wind_den = winding_number(den_vals)
    return {
        "C_abs": float(abs(C)),
        "C_phase": float(mp.arg(C)),
        "margin": float(margin),
        "ratio_spread": float(spread),
        "winding_F_side": float(wind_num),
        "winding_Xi_side": float(wind_den),
    }


# ----------------------------------------------------------------------------
# Hauptlauf
# ----------------------------------------------------------------------------

def run_for_N(c, N, T, dps, M=48):
    build = _import_ccv()
    mp.mp.dps = dps
    t0 = time.time()
    print(f"\n=== c={c}, N={N}, T={T}, dps={dps} ===")
    Q = build(c=c, N=N, T=T, dps=dps)
    t1 = time.time()
    print(f"  Galerkin-Matrix gebaut in {t1 - t0:.1f}s")
    lam1, xi_full = ground_state_even(Q)
    t2 = time.time()
    print(f"  even-Grundzustand: lambda_1 = {mp.nstr(lam1, 8)}  ({t2 - t1:.1f}s)")

    factor, err1, zs = pick_pole_convention(xi_full, N, c)
    poles = make_poles(N, c, factor)
    print(f"  Pol-Konvention: factor={factor} (|z1-gamma_1|={mp.nstr(err1, 6)})")

    zero_errs = []
    for k, g in enumerate(GAMMA_KNOWN[:3]):
        cand = [z for z in zs if abs(z - g) < 1]
        zero_errs.append(float(min(abs(z - g) for z in cand)) if cand else None)

    result = {
        "c": float(c), "N": N, "T": T, "dps": dps,
        "lambda_1_even": mp.nstr(lam1, 10),
        "pole_factor": factor,
        "zero_errors_gamma123": zero_errs,
        "contours": {},
    }

    # Einzel-Konturen um gamma_1..gamma_3: roh (F vs C*Xi) und pol-adjustiert
    # (H_S vs C*Xi mit S = die zwei kontur-naechsten Pole, beide aussen)
    for k, g in enumerate(GAMMA_KNOWN[:3]):
        dist_pole = min(abs(g - p) for p in poles)
        neighbors = [GAMMA_KNOWN[i] for i in (k - 1, k + 1) if 0 <= i < len(GAMMA_KNOWN)]
        dist_g = min(abs(g - h) for h in neighbors)
        radius = mp.mpf("0.6") * min(dist_pole, dist_g)
        r_raw = rouche_on_contour(xi_full, poles, g, radius, M=M)
        near = sorted(poles, key=lambda p: abs(p - g))[:2]
        r_adj = rouche_on_contour(xi_full, poles, g, radius, M=M, pole_adjust=near)
        rec = {
            "center": float(g), "radius": float(radius),
            "dist_nearest_pole": float(dist_pole),
            "raw": r_raw,
            "pole_adjusted": r_adj,
            "adjust_poles": [float(p) for p in near],
        }
        result["contours"][f"single_gamma{k + 1}"] = rec
        print(f"  Kontur gamma_{k + 1} (r={float(radius):.3f}): "
              f"margin raw={r_raw['margin']:.4g} adj={r_adj['margin']:.4g}, "
              f"winding F={r_raw['winding_F_side']:.2f}/Xi={r_raw['winding_Xi_side']:.2f}, "
              f"C-spread raw={r_raw['ratio_spread']:.4g} adj={r_adj['ratio_spread']:.4g}")

    # Paket-Kontur um {gamma_1, gamma_2}: H_S vs C*Xi, S = innere Pole + nahe aeussere
    center = (GAMMA_KNOWN[0] + GAMMA_KNOWN[1]) / 2
    radius = (GAMMA_KNOWN[1] - GAMMA_KNOWN[0]) / 2 + mp.mpf("1.0")
    # Radius ggf. schrumpfen bis kein Pol AUF der Kontur (Abstand > 0.15)
    for _ in range(40):
        if all(abs(abs(p - center) - radius) > mp.mpf("0.15") for p in poles) \
           and all(abs(g - center) < radius - mp.mpf("0.3") for g in GAMMA_KNOWN[:2]) \
           and abs(GAMMA_KNOWN[2] - center) > radius + mp.mpf("0.3"):
            break
        radius -= mp.mpf("0.05")
    poles_inside = [p for p in poles if abs(p - center) < radius]
    best = None
    best_S = None
    for win in ("1.5", "3.0", "6.0", "12.0"):
        poles_near_out = [p for p in poles
                          if radius < abs(p - center) < radius + mp.mpf(win)]
        S = poles_inside + poles_near_out
        r = rouche_on_contour(xi_full, poles, center, radius, M=M * 2, pole_adjust=S)
        r["near_out_window"] = float(mp.mpf(win))
        r["n_poles_near_outside"] = len(poles_near_out)
        print(f"  Paket-Kontur g1+g2 (r={float(radius):.3f}, {len(poles_inside)} Pole innen, "
              f"Fenster {win}: {len(poles_near_out)} aussen kompensiert): "
              f"margin={r['margin']:.4g}, winding H={r['winding_F_side']:.2f} (erwartet 2), "
              f"winding Xi={r['winding_Xi_side']:.2f}, C-spread={r['ratio_spread']:.4g}")
        if best is None or r["margin"] < best["margin"]:
            best = r
            best_S = S
    best["center"] = float(center); best["radius"] = float(radius)
    best["poles_inside"] = [float(p) for p in poles_inside]
    result["contours"]["packet_gamma12"] = best

    # ------------------------------------------------------------------
    # Bulk-Faktor-Modell B(z) (Proof-Note 2026-06-09 §4 Folgeschritt 2):
    # paare aeussere F-Wurzeln w_k mit aeusseren Zeta-Nullstellen gamma_k
    # (k>=3) und messe kumulativ, ob H_S/(C * B_K * Xi) auf der Paket-
    # Kontur konstant wird. B_K = prod_{k=3..K} (z^2-w_k^2)/(z^2-g_k^2).
    # Alle Faktoren windungsneutral (Wurzeln+Pole ausserhalb der Kontur).
    # ------------------------------------------------------------------
    pairs_out = []
    for k in range(2, len(GAMMA_KNOWN)):
        g = GAMMA_KNOWN[k]
        cand = [z for z in zs if abs(z - g) < mp.mpf("1.0")]
        if cand:
            w = min(cand, key=lambda z: abs(z - g))
            pairs_out.append((k + 1, w, g))
    if pairs_out and best_S is not None:
        bulk_records = []
        print(f"  Bulk-Modell B(z): {len(pairs_out)} aeussere Wurzel/Nullstellen-Paare "
              f"(k={pairs_out[0][0]}..{pairs_out[-1][0]}), S=bestes Fenster "
              f"({best.get('near_out_window')})")
        for K_keep in range(1, len(pairs_out) + 1):
            w_out = [p[1] for p in pairs_out[:K_keep]]
            g_out = [p[2] for p in pairs_out[:K_keep]]
            r = rouche_on_contour(xi_full, poles, center, radius, M=M * 2,
                                  pole_adjust=best_S, bulk_model=(w_out, g_out))
            rec = {"K_outer": K_keep,
                   "k_indices": [p[0] for p in pairs_out[:K_keep]],
                   "margin": r["margin"], "ratio_spread": r["ratio_spread"],
                   "C_abs": r["C_abs"],
                   "winding_H": r["winding_F_side"], "winding_Xi": r["winding_Xi_side"]}
            bulk_records.append(rec)
            print(f"    B_K mit K={K_keep:2d} aeusseren Paaren: margin={r['margin']:.4g}, "
                  f"C-spread={r['ratio_spread']:.4g}")
        result["contours"]["packet_gamma12_bulkmodel"] = {
            "pairs_outer": [{"k": p[0], "w": mp.nstr(p[1], 20), "gamma": mp.nstr(p[2], 20),
                             "w_minus_gamma": float(p[1] - p[2])} for p in pairs_out],
            "S_window": best.get("near_out_window"),
            "baseline_margin_constC": best["margin"],
            "sweep": bulk_records,
        }
    else:
        print("  Bulk-Modell-Test uebersprungen (keine aeusseren Paare gefunden)")

    # Expliziter Bulk-Faktor-Test (Proof-Note Folgeschritt 2): vergleiche das
    # lokale Wurzel-Polynom (z-w1)(z-w2) der beiden Paket-Wurzeln von F direkt
    # gegen C*Xi(z). Alles andere steckt per Definition exakt im Bulk-Faktor
    # B = H/((z-w1)(z-w2)); margin_polyB < 1 = Paket-Rouché mit explizitem B.
    w_packet = []
    for g in GAMMA_KNOWN[:2]:
        cand = [z for z in zs if abs(z - g) < 1]
        if cand:
            w_packet.append(min(cand, key=lambda z: abs(z - g)))
    if len(w_packet) == 2:
        K = contour_circle(center, radius, M * 2)
        polyv = [(z - w_packet[0]) * (z - w_packet[1]) for z in K]
        denv = [Xi_completed(z) for z in K]
        logs = [mp.log(polyv[i] / denv[i]) for i in range(len(K))]
        mean_log = sum(l.real for l in logs) / len(K)
        mean_phase = mp.atan2(sum(mp.sin(l.imag) for l in logs) / len(K),
                              sum(mp.cos(l.imag) for l in logs) / len(K))
        C0 = mp.exp(mp.mpc(mean_log, mean_phase))
        ratios = [polyv[i] / (C0 * denv[i]) for i in range(len(K))]
        marginB = max(abs(r - 1) for r in ratios)
        spreadB = max(abs(r) for r in ratios) / min(abs(r) for r in ratios)
        recB = {"margin": float(marginB), "ratio_spread": float(spreadB),
                "C_abs": float(abs(C0)),
                "w_packet": [mp.nstr(w, 20) for w in w_packet],
                "center": float(center), "radius": float(radius)}
        result["contours"]["packet_gamma12_polyB"] = recB
        print(f"  Paket-PolyB-Test (z-w1)(z-w2) vs C*Xi: margin={recB['margin']:.4g}, "
              f"C-spread={recB['ratio_spread']:.4g}")
    else:
        print("  Paket-PolyB-Test uebersprungen (Paket-Wurzeln nicht gefunden)")

    result["runtime_s"] = round(time.time() - t0, 1)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--c", type=float, default=14.0)
    ap.add_argument("--N-list", type=str, default="16,24,32")
    ap.add_argument("--T", type=int, default=200)
    ap.add_argument("--dps", type=int, default=60)
    ap.add_argument("--dps-scale", type=float, default=0.0,
                    help="wenn >0: dps = max(--dps, dps_scale*N) pro N (Floor-Hygiene, Proof-Note §4.1)")
    ap.add_argument("--M", type=int, default=48)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    N_list = [int(x) for x in args.N_list.split(",")]
    all_results = []
    for N in N_list:
        dps_N = max(args.dps, int(args.dps_scale * N)) if args.dps_scale > 0 else args.dps
        all_results.append(run_for_N(args.c, N, args.T, dps_N, M=args.M))

    stamp = time.strftime("%Y-%m-%d")
    out = args.out or os.path.join(
        RESULTS_DIR, f"MS1_DET_ZERO_COUNTING_c{int(args.c)}_{stamp}.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump({"script": os.path.basename(__file__),
                   "purpose": "MS1 det-zero-counting Rouché-Erstkartierung (Programmschritt 3)",
                   "results": all_results}, f, indent=2)
    print(f"\nRohdaten: {out}")

    print("\n=== N-Trend (Rouché-Margin, pol-adjustiert) ===")
    print(f"{'N':>4} {'g1 adj':>12} {'g2 adj':>12} {'g3 adj':>12} {'Paket g1+g2':>14}")
    for r in all_results:
        cs = r["contours"]
        print(f"{r['N']:>4} {cs['single_gamma1']['pole_adjusted']['margin']:>12.4g} "
              f"{cs['single_gamma2']['pole_adjusted']['margin']:>12.4g} "
              f"{cs['single_gamma3']['pole_adjusted']['margin']:>12.4g} "
              f"{cs['packet_gamma12']['margin']:>14.4g}")


if __name__ == "__main__":
    main()
