#!/usr/bin/env python3
# coding: utf-8
"""
anti_diagonal_route -- Antest der Anti-Diagonal-Route (Thm 4.1) gegen die
Auslöschungsstruktur des Sektor-Gaps.

BRANCH: explore/anti-diagonal-route (Explorationszweig, nicht main)
STATUS: Antest. Kein Beweis, kein Claim-Upgrade. Claim-Ceiling unverändert
(diagnostischer Atlas, keine prädiktive GRH-Behauptung, kein C2-Upgrade).

-------------------------------------------------------------------------------
AUSGANGSLAGE
-------------------------------------------------------------------------------
Der Lauf vom 2026-08-01 (_proof-notes/ETA_SOURCE_PROJECTION_2026-08-01.md) hat
gezeigt: der Gap entsteht als Differenz dreier Terme, die jeder ein Vielfaches
von ihm sind (Median-Toleranz 0.53, Extremfall chi_60 0.007). Jede Näherung, die
einen Term approximiert, verfehlt deshalb das Vorzeichen -- das erklärt v1..v4.

ANALYTIC_PIPELINE.md Thm 4.1 verspricht einen Ausweg: der gesamte Gap-Beitrag
sitze im anti-diagonalen Prim-Reflexionskern kappa_chi(x+y), also in einem
Objekt, aus dem der charakter-unabhängige Bulk bereits herausgekürzt ist.

-------------------------------------------------------------------------------
DIE ROUTE IN MATRIXFORM
-------------------------------------------------------------------------------
In build_W gilt (Zeile "ovl = I1 + I2 if sector == 'cos' else I1 - I2"):
    W_cos enthält w*(I1 + I2),   W_sin enthält w*(I1 - I2)
mit I1 = I1(n-m) ("diagonal") und I2 = I2(n+m, Phase in m) ("anti-diagonal").

Die Basen sind gegeneinander verschoben: cos-Matrixindex i entspricht Mode i,
sin-Matrixindex j entspricht Mode j+1. Auf den GEMEINSAMEN Moden 1..N-1 gilt
daher mit A := W_sin[:-1,:-1] und B := W_cos[1:,1:]

    D := B - A = 2 * sum_p w_p * I2 * NM      (rein anti-diagonal)

Archimedische Diagonale und I1-Anteil heben sich exakt weg, weil sie bei
gleicher Mode identisch sind. Der Bulk ist damit ANALYTISCH eliminiert, nicht
numerisch subtrahiert -- genau das, was die Auslöschung vermeiden soll.

Der Gap ist gap = lambda_1(W_sin) - lambda_1(W_cos), auf gemeinsamen Moden also
lambda_1(A) - lambda_1(A+D), und in erster Ordnung

    gap ~ - <phi_A | D | phi_A>.

-------------------------------------------------------------------------------
WAS GEPRUEFT WIRD
-------------------------------------------------------------------------------
V1 STRUKTUR: Ist D exakt der reine anti-diagonale I2-Anteil? Geprüft durch
             direkten Vergleich mit einer nur aus I2 gebauten Matrix.
             Das ist eine unabhängige numerische Prüfung von Thm 4.1.
             ACHTUNG: Ein Hankel-Test (D[i,j] = f(i+j)) wäre hier FALSCH und
             würde Thm 4.1 scheinbar widerlegen -- siehe build_I2_only().
V2 BULK-ELIMINATION: Ist ||D|| klein gegen ||W_cos||?
V3 MODE-0: Wie stark verfälscht das Weglassen der cos-Mode 0 den Gap?

P1 gap_ref     = lambda_1(W_sin) - lambda_1(W_cos)          [Referenz, exakt]
P2 P_D_exact   = -<phi_A|D|phi_A>, phi_A exakter GZ von A   [1. Ordnung, 1 Eigenvektor]
P3 P_D_crude   = -<phi_c|D|phi_c>, phi_c GZ der DIAGONALE   [kein Eigenvektor -> nicht-taut.]
P4 P_D_mode    = -D[n*,n*] an der Resonanzmode              [Ein-Moden, nicht-taut.]

ENTSCHEIDUNGSKRITERIUM (Antest):
  aussichtsreich, wenn P2 die Vorzeichen über der Basisrate trifft UND P3/P4 dem
  Vorzeichen von P2 weitgehend folgen (= robust gegen grobe Näherung von phi);
  negativ, wenn P2 schon scheitert oder P3/P4 beliebig vom Vorzeichen abweichen
  (= die Auslöschung ist nur verschoben, nicht beseitigt).

EXITCODE 0 = Rechnung samt Selbstprüfungen sauber (auch bei negativem Befund).
EXITCODE 1 = Selbstprüfung verletzt, Ergebnis nicht verwertbar.
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import sympy
from scipy.special import digamma

sys.stdout.reconfigure(line_buffering=True)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

from analytic_gap_formula_test_v2 import build_W  # noqa: E402
from analytic_gap_formula_test import make_chi_D, CHARS  # noqa: E402


def results_dir() -> Path:
    """Robust: der Klon nutzt results/, der OneDrive-Arbeitsordner _results/."""
    for name in ("_results", "results"):
        p = ROOT / name
        if p.is_dir():
            return p
    raise SystemExit("FEHLER: weder _results/ noch results/ gefunden.")


RES = results_dir()
LAM = 20000


def load_empirical_gaps() -> dict[str, tuple[float, int]]:
    """Empirische Gaps bei LAM, jeweils hoechstes verfuegbares N."""
    out: dict[str, tuple[float, int]] = {}
    for fn in ("oscillating6_N600_results.json", "all10_high_N_results.json",
               "chi12_high_N_results.json"):
        path = RES / fn
        if not path.exists():
            continue
        for row in json.loads(path.read_text(encoding="utf-8")):
            if int(row.get("lambda", 0)) != LAM:
                continue
            if fn == "chi12_high_N_results.json" and int(row.get("N", 0)) != 400:
                continue
            name, n_used = row["chi"], int(row["N"])
            prev = out.get(name)
            if prev is None or n_used > prev[1]:
                out[name] = (float(row["gap"]), n_used)
    return out


def sector_W(chi_fn, N: int, sector: str) -> np.ndarray:
    L = math.log(LAM)
    primes = [p for p in sympy.primerange(2, LAM + 1) if chi_fn(p) != 0]
    return build_W(sector, N, L, primes, [chi_fn(p) for p in primes])


def build_I2_only(N: int, L: float, primes, chi_vals, sector: str) -> np.ndarray:
    """Baut NUR den anti-diagonalen I2-Anteil (k = n+m) von build_W.

    Dient der direkten Verifikation von Thm 4.1 in Matrixform: die Differenz der
    Sektormatrizen auf gemeinsamen Moden muss exakt 2*I2 sein.

    HINWEIS: Ein Test auf Hankel-Struktur (D[i,j] = f(i+j)) waere hier FALSCH.
    I2 haengt ueber die Phase ph2 = -m*pi*ds/L separat von m ab und ist in der
    Modenbasis daher nicht Toeplitz-in-(n+m), obwohl der zugrunde liegende Kern
    kappa_chi(x+y) im Ortsraum anti-diagonal ist.
    """
    bs = 0 if sector == "cos" else 1
    idx = np.arange(N) + bs
    nn, mm = idx[:, None], idx[None, :]
    W = np.zeros((N, N))
    norm = (np.where(idx > 0, 1.0 / np.sqrt(L), 1.0 / np.sqrt(2 * L))
            if sector == "cos" else np.full(N, 1.0 / np.sqrt(L)))
    NM = np.outer(norm, norm)

    def A(k, phi, t, Lv):
        r = np.zeros_like(k, dtype=float)
        mk = (k != 0)
        r[mk] = (Lv / (k[mk] * np.pi)) * np.sin(k[mk] * np.pi * t / Lv + phi[mk])
        r[~mk] = np.cos(phi[~mk]) * t
        return r

    for p, cp in zip(primes, chi_vals):
        if cp == 0:
            continue
        lp = np.log(p)
        for me in range(1, int(2 * L / lp) + 1):
            d = me * lp
            if d >= 2 * L:
                break
            w = (cp ** me) * lp / (p ** (me / 2.0))
            for ds in (d, -d):
                a_, b_ = max(-L, -L + ds), min(L, L + ds)
                if a_ >= b_:
                    continue
                k2 = nn + mm
                P2 = np.broadcast_to(-mm * np.pi * ds / L, (N, N)).copy()
                I2 = 0.5 * (A(k2, P2, b_, L) - A(k2, P2, a_, L))
                W += w * I2 * NM
    return 0.5 * (W + W.T)


def sgn(x: float) -> int:
    return (x > 0) - (x < 0)


def analyse(name: str, D_char: int, N: int) -> dict:
    chi = make_chi_D(D_char)
    L = math.log(LAM)
    primes = [p for p in sympy.primerange(2, LAM + 1) if chi(p) != 0]
    chi_vals = [chi(p) for p in primes]
    Wc = build_W("cos", N, L, primes, chi_vals)
    Ws = build_W("sin", N, L, primes, chi_vals)

    ev_c = float(np.linalg.eigvalsh(Wc)[0])
    ev_s = float(np.linalg.eigvalsh(Ws)[0])
    gap_ref = ev_s - ev_c

    # Gemeinsame Moden 1..N-1
    A = Ws[:-1, :-1]
    B = Wc[1:, 1:]
    D = B - A

    evA, vecA = np.linalg.eigh(A)
    phi_A = vecA[:, 0]

    # P2: erste Ordnung mit exaktem Grundzustand von A
    P_D_exact = -float(phi_A @ D @ phi_A)

    # P3: Grundzustand der reinen DIAGONALE von A (kein Eigenvektor noetig)
    n_crude = int(np.argmin(np.diag(A)))
    P_D_crude = -float(D[n_crude, n_crude])

    # P4: Resonanzmode aus der vollen cos-Matrix (wie im Vorlauf definiert)
    n_res = int(np.argmin(np.diag(B)))
    P_D_mode = -float(D[n_res, n_res])

    # V1: Ist D exakt der reine anti-diagonale I2-Anteil? (Thm 4.1 in Matrixform)
    I2c = build_I2_only(N, L, primes, chi_vals, "cos")
    D_pred = 2.0 * I2c[1:, 1:]
    i2_residual = float(np.linalg.norm(D - D_pred) / max(np.linalg.norm(D), 1e-30))

    # V3: Mode-0-Effekt der Restriktion
    gap_restricted = float(np.linalg.eigvalsh(A)[0]) - float(np.linalg.eigvalsh(B)[0])

    return {
        "chi": name, "D": D_char, "N": N,
        "gap_ref": gap_ref, "gap_restricted": gap_restricted,
        "ev_cos": ev_c, "ev_sin": ev_s,
        "P_D_exact": P_D_exact, "P_D_crude": P_D_crude, "P_D_mode": P_D_mode,
        "norm_D": float(np.linalg.norm(D)), "norm_Wcos": float(np.linalg.norm(Wc)),
        "bulk_ratio": float(np.linalg.norm(D) / np.linalg.norm(Wc)),
        "i2_residual": i2_residual,
        "n_crude": n_crude, "n_res": n_res,
        "phi_A_ipr": float(1.0 / np.sum(phi_A ** 4)),
    }


def main() -> int:
    quick = "--quick" in sys.argv
    N = 80 if quick else 200
    subset = {"chi_8", "chi_12", "chi_13", "chi_33"} if quick else None

    print("=" * 78)
    print("ANTI-DIAGONAL-ROUTE -- Antest (Branch explore/anti-diagonal-route)")
    print("=" * 78)
    print(f"lambda={LAM}  N={N}  Modus={'QUICK' if quick else 'VOLL'}")
    print("Status: Antest. Kein Beweis, kein Claim-Upgrade.\n")

    emp = load_empirical_gaps()
    rows, failures = [], []

    for name, Dc in CHARS:
        if name not in emp:
            continue
        if subset and name not in subset:
            continue
        t0 = time.time()
        r = analyse(name, Dc, N)
        r["gap_emp"] = emp[name][0]
        r["N_emp"] = emp[name][1]
        rows.append(r)
        print(f"[{name:>7}] gap_emp={r['gap_emp']:+8.5f} gap_ref={r['gap_ref']:+8.5f} | "
              f"P2={r['P_D_exact']:+9.5f} P3={r['P_D_crude']:+9.5f} "
              f"P4={r['P_D_mode']:+9.5f} | ||D||/||W||={r['bulk_ratio']:.4f} "
              f"I2_res={r['i2_residual']:.1e} [{time.time()-t0:.0f}s]")

    if not rows:
        print("FEHLER: keine Zeilen.")
        return 1

    # --- Selbstpruefung: Restriktion darf den Gap nicht umdrehen
    for r in rows:
        if sgn(r["gap_restricted"]) != sgn(r["gap_ref"]):
            failures.append(f"{r['chi']}: Restriktion auf gemeinsame Moden dreht das "
                            f"Gap-Vorzeichen ({r['gap_restricted']:+.5f} vs "
                            f"{r['gap_ref']:+.5f}) -- Mode-0-Effekt nicht vernachlaessigbar")

    emps = [r["gap_emp"] for r in rows]
    n = len(rows)
    n_pos = sum(1 for g in emps if g > 0)
    base = max(n_pos, n - n_pos)

    print("\n" + "=" * 78)
    print("STRUKTURPRUEFUNG")
    print("=" * 78)
    hd = [r["i2_residual"] for r in rows]
    br = [r["bulk_ratio"] for r in rows]
    print(f"  V1 Residuum ||D - 2*I2|| / ||D||: max={max(hd):.2e}")
    print(f"     -> Thm 4.1 in Matrixform {'EXAKT BESTAETIGT' if max(hd) < 1e-10 else 'NICHT exakt'}")
    print(f"  V2 ||D||/||W_cos||: min={min(br):.4f} max={max(br):.4f}  "
          f"(klein = Bulk eliminiert)")
    print(f"  V3 Mode-0-Effekt: Vorzeichen erhalten bei "
          f"{sum(1 for r in rows if sgn(r['gap_restricted']) == sgn(r['gap_ref']))}/{n}")

    print("\n" + "=" * 78)
    print(f"VORZEICHEN-BILANZ (Basisrate {base}/{n})")
    print("=" * 78)
    summary = {}
    for key, label in [("gap_ref", "Referenz: exakter Sektor-Gap"),
                       ("P_D_exact", "P2 erste Ordnung, exakter GZ von A"),
                       ("P_D_crude", "P3 Diagonal-GZ (nicht-tautologisch)"),
                       ("P_D_mode", "P4 Resonanzmode (nicht-tautologisch)")]:
        vals = [r[key] for r in rows]
        ok = sum(1 for v, g in zip(vals, emps) if sgn(v) == sgn(g))
        changes = 0 < sum(1 for v in vals if v > 0) < n
        informative = changes and ok > base
        summary[key] = {"label": label, "sign_ok": ok, "n": n,
                        "sign_changes": changes, "informative": informative}
        print(f"  {key:>11} {ok:>2}/{n}  wechselt={str(changes):>5}  "
              f"informativ={'JA' if informative else 'nein':>4}   {label}")

    # Robustheit: folgen die groben Praediktoren dem Vorzeichen von P2?
    agree_crude = sum(1 for r in rows if sgn(r["P_D_crude"]) == sgn(r["P_D_exact"]))
    agree_mode = sum(1 for r in rows if sgn(r["P_D_mode"]) == sgn(r["P_D_exact"]))
    print(f"\n  ROBUSTHEIT: P3 folgt P2 in {agree_crude}/{n}, P4 folgt P2 in {agree_mode}/{n}")
    print("  (hohe Uebereinstimmung = Vorzeichen haengt nicht am genauen Eigenvektor)")

    verdict = ("AUSSICHTSREICH" if summary["P_D_exact"]["informative"]
               and agree_crude >= n - 1 else "NEGATIV")
    print(f"\n  ANTEST-VERDIKT: {verdict}")

    out = RES / f"ANTI_DIAGONAL_ROUTE_{'quick' if quick else 'full'}_2026-08-02.json"
    out.write_text(json.dumps({
        "meta": {"branch": "explore/anti-diagonal-route", "lambda": LAM, "N": N,
                 "mode": "quick" if quick else "full",
                 "status": "Antest; kein Beweis, kein Claim-Upgrade",
                 "self_checks_passed": not failures, "failures": failures},
        "rows": rows, "summary": summary,
        "structure": {"i2_residual_max": max(hd), "bulk_ratio_min": min(br),
                      "bulk_ratio_max": max(br)},
        "robustness": {"crude_follows_exact": agree_crude,
                       "mode_follows_exact": agree_mode, "n": n},
        "verdict": verdict,
    }, indent=2), encoding="utf-8")
    print(f"\n[out] {out}")

    if failures:
        print("\nSELBSTPRUEFUNG FEHLGESCHLAGEN:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("\nSelbstpruefung bestanden.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
