#!/usr/bin/env python3
# coding: utf-8
"""
anti_diagonal_rank -- Rang- und Supportfeld der anti-diagonalen Primquelle.

BRANCH: explore/anti-diagonal-route
STATUS: Diagnose. Kein Beweis, kein Claim-Upgrade. Claim-Ceiling unverändert
(diagnostischer Atlas, keine prädiktive GRH-Behauptung, kein C2-Upgrade).

Bearbeitet: BEWEISNOTIZ.md, offene proof-grade Aufgabe 3 `anti_diagonal_rank`
("Support-/Rangfeld für die anti-diagonale Primquelle formulieren").

-------------------------------------------------------------------------------
AUSGANGSLAGE
-------------------------------------------------------------------------------
Der Antest vom 2026-08-02 (_proof-notes/ANTI_DIAGONAL_ROUTE_2026-08-02.md) hat
belegt: auf den gemeinsamen Moden ist

    D := W_cos[1:,1:] - W_sin[:-1,:-1]

exakt die reine anti-diagonale Primquelle (Residuum 1.2e-13), und der Bulk ist
damit eliminiert (||D||/||W_cos|| = 0.068). Der Vorhersageteil scheiterte aber:
<phi|D|phi> haengt weiter empfindlich vom Eigenvektor ab, weil D klein, aber
stark oszillierend ist. Die Ausloeschung ist verschoben, nicht beseitigt.

LEITFRAGE DIESES LAUFS
Hat D eine ausnutzbare Struktur? Konkret: Wenn D niedrigen effektiven Rang hat,
koennte eine Rangreduktion D_k den oszillierenden Anteil abschneiden und das
Vorzeichen gegen Eigenvektor-Fehler STABILISIEREN. Dann waere die Empfindlichkeit
nicht nur verschoben, sondern beherrschbar.

-------------------------------------------------------------------------------
WAS GERECHNET WIRD
-------------------------------------------------------------------------------
R1 RANG      Singulaerwertspektrum von D; effektiver Rang bei 90% und 99% Energie.
R2 SUPPORT   Verteilung der Zeilennormen ueber die Moden: auf wie viele Moden
             konzentriert sich die Masse (IPR, Anzahl Moden fuer 90%)?
R3 PHASE     Anzahl Vorzeichenwechsel des fuehrenden Singulaervektors -- Mass
             fuer die Oszillation, die die Mittelung zerstoert.
R4 STABIL    Praediktoren P_k = -<phi|D_k|phi> fuer k = 1,2,3,5,10,voll, jeweils
             mit exaktem Grundzustand phi_A UND mit dem nicht-tautologischen
             Diagonal-Grundzustand phi_c. Entscheidend ist, ob es ein k gibt, bei
             dem die NICHT-TAUTOLOGISCHE Variante ueber die Basisrate kommt.

-------------------------------------------------------------------------------
BEWERTUNGSREGELN (beide verbindlich aus dem Lauf vom 2026-08-01)
-------------------------------------------------------------------------------
BASISRATEN-GATE: Eine Vorzeichenquote ist nur dann ein Signal, wenn der
Praediktor das Vorzeichen WECHSELT und die Basisrate UEBERSCHREITET. Die
Basisrate wird hier aus der jeweils bewerteten Teilmenge berechnet.

VORAB-FEHLERSCHRANKE: chi_21 und chi_60 haben tol < 0.02, d.h. ihr Galerkin-Gap
weicht bei N=200 um mehr als eine Groessenordnung vom empirischen Wert ab. Sie
werden aus der BEWERTUNG ausgeschlossen (die Rechnung laeuft mit, weil die
Rangstruktur dort weiterhin aussagekraeftig ist). Beide Bilanzen werden
ausgewiesen, damit der Ausschluss nachpruefbar bleibt.

EXITCODE 0 = Rechnung und Selbstpruefungen sauber (auch bei negativem Befund).
EXITCODE 1 = Selbstpruefung verletzt.
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import sympy

sys.stdout.reconfigure(line_buffering=True)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

from analytic_gap_formula_test_v2 import build_W  # noqa: E402
from analytic_gap_formula_test import make_chi_D, CHARS  # noqa: E402

LAM = 20000
N_GAL = 200
STAMP = "2026-08-02"
K_VALUES = [1, 2, 3, 5, 10, 20]
# tol < 0.02 aus ETA_SOURCE_PROJECTION_EVAL_2026-08-01.json
EXCLUDED = {"chi_21", "chi_60"}


def results_dir() -> Path:
    for name in ("_results", "results"):
        p = ROOT / name
        if p.is_dir():
            return p
    raise SystemExit("FEHLER: weder _results/ noch results/ gefunden.")


RES = results_dir()


def load_empirical_gaps() -> dict[str, tuple[float, int]]:
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


def sgn(x: float) -> int:
    return (x > 0) - (x < 0)


def sign_changes(v: np.ndarray, tol: float = 1e-12) -> int:
    s = np.sign(v[np.abs(v) > tol])
    return int(np.sum(s[1:] != s[:-1])) if s.size > 1 else 0


def analyse(name: str, D_char: int, N: int) -> dict:
    chi = make_chi_D(D_char)
    L = math.log(LAM)
    primes = [p for p in sympy.primerange(2, LAM + 1) if chi(p) != 0]
    cv = [chi(p) for p in primes]
    Wc = build_W("cos", N, L, primes, cv)
    Ws = build_W("sin", N, L, primes, cv)

    gap_ref = float(np.linalg.eigvalsh(Ws)[0] - np.linalg.eigvalsh(Wc)[0])

    A = Ws[:-1, :-1]
    D = Wc[1:, 1:] - A

    # R1: Singulaerwerte (D ist symmetrisch -> |Eigenwerte|, aber SVD ist robust)
    U, S, Vt = np.linalg.svd(D)
    energy = np.cumsum(S ** 2) / max(float(np.sum(S ** 2)), 1e-300)
    rank90 = int(np.searchsorted(energy, 0.90) + 1)
    rank99 = int(np.searchsorted(energy, 0.99) + 1)

    # R2: Support ueber Moden (Zeilennormen)
    rn = np.linalg.norm(D, axis=1)
    w = rn ** 2 / max(float(np.sum(rn ** 2)), 1e-300)
    ipr = float(1.0 / np.sum(w ** 2))
    order = np.argsort(-w)
    modes90 = int(np.searchsorted(np.cumsum(w[order]), 0.90) + 1)
    peak_mode = int(order[0])

    # R3: Oszillation des fuehrenden Singulaervektors
    sc_u1 = sign_changes(U[:, 0])

    # R4: Praediktoren mit rangreduziertem D_k
    evA, vecA = np.linalg.eigh(A)
    phi_exact = vecA[:, 0]
    n_crude = int(np.argmin(np.diag(A)))
    phi_crude = np.zeros_like(phi_exact)
    phi_crude[n_crude] = 1.0

    preds: dict[str, float] = {}
    for k in K_VALUES:
        kk = min(k, len(S))
        Dk = (U[:, :kk] * S[:kk]) @ Vt[:kk]
        preds[f"P_exact_k{k}"] = -float(phi_exact @ Dk @ phi_exact)
        preds[f"P_crude_k{k}"] = -float(phi_crude @ Dk @ phi_crude)
    preds["P_exact_full"] = -float(phi_exact @ D @ phi_exact)
    preds["P_crude_full"] = -float(phi_crude @ D @ phi_crude)

    # Selbstpruefung: volle Rangrekonstruktion muss D treffen
    D_rec = (U * S) @ Vt
    svd_residual = float(np.linalg.norm(D - D_rec) / max(np.linalg.norm(D), 1e-300))

    return {
        "chi": name, "D": D_char, "N": N,
        "gap_ref": gap_ref,
        "sv_top10": [float(x) for x in S[:10]],
        "rank90": rank90, "rank99": rank99, "sv_max": float(S[0]),
        "sv_decay_ratio": float(S[1] / S[0]) if S[0] > 0 else float("nan"),
        "support_ipr": ipr, "modes_for_90pct": modes90, "peak_mode": peak_mode,
        "sign_changes_u1": sc_u1,
        "norm_D": float(np.linalg.norm(D)),
        "svd_residual": svd_residual,
        "n_crude": n_crude,
        **preds,
    }


def random_sign_control(D_char: int, N: int, n_runs: int = 4, seed: int = 20260802) -> dict:
    """ZUFALLSVORZEICHEN-KONTROLLE -- das Struktur-Analogon zum Basisraten-Gate.

    Ersetzt chi(p) durch i.i.d. Zufallsvorzeichen +-1 bei identischen Primzahlen
    und Gewichten. Ein Strukturmerkmal (Rang, Oszillation), das dabei UNVERAENDERT
    bleibt, ist eine Eigenschaft der Galerkin-Overlap-Konstruktion und KEIN
    charakterspezifischer Befund.

    Motivation: Genau wie eine Vorzeichenquote gegen die Basisrate zu messen ist,
    muss ein Strukturbefund gegen diese Kontrolle gemessen werden. Ohne sie wird
    ein generisches Konstruktionsmerkmal als zahlentheoretisches Signal gelesen.
    """
    chi = make_chi_D(D_char)
    L = math.log(LAM)
    primes = [p for p in sympy.primerange(2, LAM + 1) if chi(p) != 0]
    rng = np.random.default_rng(seed)

    def ranks(vals):
        Wc = build_W("cos", N, L, primes, vals)
        Ws = build_W("sin", N, L, primes, vals)
        Dm = Wc[1:, 1:] - Ws[:-1, :-1]
        S = np.linalg.svd(Dm, compute_uv=False)
        e = np.cumsum(S ** 2) / max(float(np.sum(S ** 2)), 1e-300)
        return (int(np.searchsorted(e, 0.90) + 1), int(np.searchsorted(e, 0.99) + 1),
                float(np.linalg.norm(Dm)))

    r90, r99, nrm = ranks([chi(p) for p in primes])
    controls = [ranks(list(rng.choice([-1, 1], size=len(primes)))) for _ in range(n_runs)]
    ctrl_r90 = [c[0] for c in controls]
    ctrl_r99 = [c[1] for c in controls]
    ctrl_nrm = [c[2] for c in controls]
    return {
        "chi": f"chi_{D_char}", "N": N, "n_runs": n_runs,
        "real": {"rank90": r90, "rank99": r99, "norm_D": nrm},
        "control": [{"rank90": a, "rank99": b, "norm_D": c} for a, b, c in controls],
        # Getrennt ausgewiesen -- ein Unterschied in EINER Kennzahl macht die
        # anderen nicht signifikant.
        "rank90_differs": bool(r90 not in ctrl_r90),
        "rank99_differs": bool(r99 not in ctrl_r99),
        "norm_below_all_controls": bool(nrm < min(ctrl_nrm)),
        "norm_ratio_to_control_mean": float(nrm / (sum(ctrl_nrm) / len(ctrl_nrm))),
        "note": "Bleibt eine Kennzahl unter Zufallsvorzeichen gleich, ist sie ein "
                "Konstruktionsmerkmal und KEIN charakterspezifischer Befund. "
                "rank90 ist erfahrungsgemaess identisch (kein Signal); ein "
                "Unterschied in rank99 allein traegt wenig. Aussagekraeftiger ist "
                "die Norm: eine systematisch kleinere anti-diagonale Masse als bei "
                "Zufallsvorzeichen waere echte Charakter-Cancellation. Wenige Laeufe "
                "-- als Hinweis lesen, nicht als Befund.",
    }


def bilanz(rows: list[dict], keys: list[str], label: str) -> dict:
    emps = [r["gap_emp"] for r in rows]
    n = len(rows)
    n_pos = sum(1 for g in emps if g > 0)
    base = max(n_pos, n - n_pos)
    print(f"\n  --- {label}  (n={n}, Basisrate {base}/{n}) ---")
    out = {"n": n, "base_rate": base, "predictors": {}}
    for key in keys:
        vals = [r[key] for r in rows]
        ok = sum(1 for v, g in zip(vals, emps) if sgn(v) == sgn(g))
        changes = 0 < sum(1 for v in vals if v > 0) < n
        informative = changes and ok > base
        out["predictors"][key] = {"sign_ok": ok, "sign_changes": changes,
                                  "informative": informative}
        mark = "  <== INFORMATIV" if informative else ""
        print(f"    {key:>16} {ok:>2}/{n}  wechselt={str(changes):>5}{mark}")
    return out


def main() -> int:
    if "--control" in sys.argv:
        print("=" * 78)
        print("ZUFALLSVORZEICHEN-KONTROLLE (Struktur-Analogon zum Basisraten-Gate)")
        print("=" * 78)
        out = []
        for Dc in (12, 33):
            c = random_sign_control(Dc, N=80)
            out.append(c)
            print(f"chi_{Dc}: real rank90={c['real']['rank90']} rank99={c['real']['rank99']} "
                  f"||D||={c['real']['norm_D']:.3f}")
            for i, ct in enumerate(c["control"], 1):
                print(f"   Zufall {i}: rank90={ct['rank90']} rank99={ct['rank99']} "
                      f"||D||={ct['norm_D']:.3f}")
            print(f"   -> rank90 unterscheidet sich: {c['rank90_differs']} | "
                  f"rank99: {c['rank99_differs']} | ||D|| unter allen Kontrollen: "
                  f"{c['norm_below_all_controls']} "
                  f"(Verhaeltnis {c['norm_ratio_to_control_mean']:.2f})")
        path = RES / f"ANTI_DIAGONAL_RANK_CONTROL_{STAMP}.json"
        path.write_text(json.dumps({"meta": {"stamp": STAMP, "mode": "random-sign control"},
                                    "results": out}, indent=2), encoding="utf-8")
        print(f"[out] {path}")
        return 0

    print("=" * 78)
    print("ANTI-DIAGONAL-RANK -- Rang- und Supportfeld der anti-diagonalen Quelle")
    print("=" * 78)
    print(f"lambda={LAM}  N={N_GAL}  Aufgabe 3 aus BEWEISNOTIZ.md")
    print("Status: Diagnose. Kein Beweis, kein Claim-Upgrade.")
    print(f"Aus der Bewertung ausgeschlossen (tol<0.02): {sorted(EXCLUDED)}\n")

    emp = load_empirical_gaps()
    rows, failures = [], []

    for name, Dc in CHARS:
        if name not in emp:
            continue
        t0 = time.time()
        r = analyse(name, Dc, N_GAL)
        r["gap_emp"], r["N_emp"] = emp[name]
        r["excluded_from_scoring"] = name in EXCLUDED
        rows.append(r)
        if r["svd_residual"] > 1e-10:
            failures.append(f"{name}: SVD-Rekonstruktion ungenau ({r['svd_residual']:.2e})")
        flag = " [aus Bewertung ausgeschl.]" if name in EXCLUDED else ""
        print(f"[{name:>7}] gap_emp={r['gap_emp']:+8.5f} | rank90={r['rank90']:>3} "
              f"rank99={r['rank99']:>3} sv2/sv1={r['sv_decay_ratio']:.3f} "
              f"IPR={r['support_ipr']:.1f} Moden90%={r['modes_for_90pct']:>3} "
              f"VZW(u1)={r['sign_changes_u1']:>3}{flag} [{time.time()-t0:.0f}s]")

    if not rows:
        print("FEHLER: keine Zeilen.")
        return 1

    print("\n" + "=" * 78)
    print("R1-R3  STRUKTUR DER ANTI-DIAGONALEN QUELLE")
    print("=" * 78)
    r90 = [r["rank90"] for r in rows]
    r99 = [r["rank99"] for r in rows]
    ipr = [r["support_ipr"] for r in rows]
    vzw = [r["sign_changes_u1"] for r in rows]
    print(f"  Effektiver Rang (90% Energie): min={min(r90)} max={max(r90)} "
          f"Median={int(np.median(r90))}   von {N_GAL-1} moeglichen")
    print(f"  Effektiver Rang (99% Energie): min={min(r99)} max={max(r99)} "
          f"Median={int(np.median(r99))}")
    print(f"  Support-IPR ueber Moden:       min={min(ipr):.1f} max={max(ipr):.1f}")
    print(f"  Vorzeichenwechsel in u1:       min={min(vzw)} max={max(vzw)}")
    low_rank = int(np.median(r90)) <= 5
    print(f"  -> D ist {'NIEDRIGRANGIG' if low_rank else 'NICHT niedrigrangig'} "
          f"(Kriterium: Median rank90 <= 5)")

    print("\n" + "=" * 78)
    print("R4  STABILISIERT EINE RANGREDUKTION DAS VORZEICHEN?")
    print("=" * 78)
    keys = ([f"P_exact_k{k}" for k in K_VALUES] + ["P_exact_full"]
            + [f"P_crude_k{k}" for k in K_VALUES] + ["P_crude_full"])
    scored = [r for r in rows if not r["excluded_from_scoring"]]
    b_all = bilanz(rows, keys, "ALLE Charaktere")
    b_scored = bilanz(scored, keys, "OHNE chi_21/chi_60 (Vorab-Fehlerschranke)")

    crude_inf = [k for k in keys if k.startswith("P_crude")
                 and b_scored["predictors"][k]["informative"]]
    print("\n" + "=" * 78)
    print("VERDIKT")
    print("=" * 78)
    if crude_inf:
        print(f"  AUSSICHTSREICH: nicht-tautologische Praediktoren ueber Basisrate: {crude_inf}")
    else:
        print("  NEGATIV: kein nicht-tautologischer Praediktor kommt ueber die Basisrate.")
        print("  Eine Rangreduktion stabilisiert das Vorzeichen nicht.")

    out = RES / f"ANTI_DIAGONAL_RANK_{STAMP}.json"
    out.write_text(json.dumps({
        "meta": {"branch": "explore/anti-diagonal-route", "lambda": LAM, "N": N_GAL,
                 "k_values": K_VALUES, "excluded_from_scoring": sorted(EXCLUDED),
                 "status": "Diagnose; kein Beweis, kein Claim-Upgrade",
                 "self_checks_passed": not failures, "failures": failures},
        "rows": rows,
        "structure": {"rank90_median": int(np.median(r90)),
                      "rank99_median": int(np.median(r99)),
                      "ipr_min": min(ipr), "ipr_max": max(ipr),
                      "sign_changes_u1_min": min(vzw), "sign_changes_u1_max": max(vzw),
                      "low_rank": low_rank},
        "bilanz_all": b_all, "bilanz_scored": b_scored,
        "verdict": "AUSSICHTSREICH" if crude_inf else "NEGATIV",
        "informative_crude": crude_inf,
    }, indent=2), encoding="utf-8")
    print(f"\n[out] {out}")

    if failures:
        print("\nSELBSTPRUEFUNG FEHLGESCHLAGEN:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("\nSelbstpruefung bestanden (SVD-Rekonstruktion exakt).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
