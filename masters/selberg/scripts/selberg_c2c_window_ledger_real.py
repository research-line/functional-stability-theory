#!/usr/bin/env python3
"""Real-data ledger for the Selberg C2c minimizer-identification gate (rev 2).

Successor of the toy run `selberg_c2c_window_ledger.py` (2026-05-27).
Rev 2 (2026-07-03) incorporates the Codex adversarial review
(`_results/CODEX_REVIEW_C2C_REAL_LEDGER_2026-07-03.md`):

  * The ledger is declared as what it is: a FINITE-SPECTRUM PROXY of the
    spectral-theorem certificate, augmented by an explicit tail budget.
    `bound_holds` on the truncated list is near-tautological and is kept
    only as an implementation sanity check; candidates are judged by
    residual/gap and off-window mass, tail-adjusted for multipliers.
  * Explicit tail budgets for multiplier candidates: with a conservative
    counting bound  dN/dlambda <= rho  (twice the Weyl main term
    Area/4pi; documented per surface) we integrate
        T0 >= sum_{lambda > Lambda} mult |h(lambda)|^2
        T2 >= sum_{lambda > Lambda} mult |h(lambda)|^2 (lambda-nu)^2
    numerically (Simpson) and report worst-case corrected off-window
    mass, residual and ratio.  Pure eigenvector candidates carry no
    unknown tail by construction.
  * Modular surface rows are labelled Delta_cusp / sigma(Delta_cusp) /
    P_I^cusp throughout: the surface is non-compact and only the
    self-adjoint restriction to L^2_cusp is addressed; the continuous
    Eisenstein spectrum [1/4, oo) is outside this ledger's scope.
  * e_comm = 0 is a CONDITIONAL statement about the modelled candidate
    class (functions h(Delta) via functional calculus); the ledger does
    not prove that an external Selberg/Connes minimizer belongs to this
    class.  The restriction of exp(-t Delta) to L^2_0 is legitimate
    because L^2_0 is Delta-invariant.
  * Controls are classified: `positive`, `radial_miscandidate`
    (genuine h(Delta), wrongly centred/tuned), `mathematical_negative`
    (adversarial spectral distributions), `bookkeeping_ablation`
    (structure breakers: not functions of Delta / wrong bookkeeping).
    Additional genuine radial fail-candidates (wide window, miscentred
    bandpass) were added so that failure is exhibited INSIDE the
    theorem-relevant candidate class as well.
  * pass/ambiguous thresholds are audit conventions, not part of the
    spectral theorem.

Verified spectra (web-verified 2026-07-03):

  * Bolza surface (compact, genus 2): first ten positive Laplace
    eigenvalues with multiplicities, rigorous certified computation by
    Strohmaier & Uski 2013 (Comm. Math. Phys. 317; arXiv:1110.2150).
    Their method finds ALL eigenvalues in a given interval, so the list
    is complete below the cutoff (bracketing/completeness assumption
    documented in the report).
  * Modular surface PSL(2,Z)\\H, level 1: first eight cuspidal spectral
    parameters r from the rigorous LMFDB Maass-form database
    (Booker/Seymour-Howell line, arXiv:2201.08760, arXiv:2502.01442);
    lambda = 1/4 + r^2.

No Selberg-RH claim.  This fills the C2c window ledger with real
spectral windows; it does not identify the Connes minimizer.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

DATE = "2026-07-03"
PASS_RATIO = 0.25       # audit convention, not part of the spectral theorem
AMBIGUOUS_RATIO = 0.50  # audit convention

# ---------------------------------------------------------------------------
# Verified spectra (web-verified 2026-07-03 against the cited sources).
# ---------------------------------------------------------------------------

# Bolza surface: (eigenvalue, multiplicity).  lambda_0 = 0 (constants)
# belongs to the spectrum and matters for outer gaps.
BOLZA_SPECTRUM: List[Tuple[float, int]] = [
    (0.0, 1),
    (3.838887258842199518586622450435464597, 3),
    (5.353601341189050410918048311031446376, 4),
    (8.249554815200658121890106450682456568, 2),
    (14.726216787788832041289318442184835983, 4),
    (15.048916133267048746181584340258811275, 3),
    (18.658819627260193806296234661340993631, 3),
    (20.519859734142002001149771260642099824, 4),
    (23.078558481381635155075206299574552996, 1),
    (28.079605737677729081562207945001124964, 3),
    (30.833042737932549674243957560470189329, 4),
]

# Modular surface, level 1, cuspidal spectral parameters r (LMFDB,
# rigorous database).  lambda = 1/4 + r^2, all multiplicity 1.
MODULAR_R: List[float] = [
    9.53369526,
    12.1730083,
    13.7797513,
    14.3585095,
    16.1380731,
    16.6442592,
    17.7385633,
    18.1809178,
]
MODULAR_SPECTRUM: List[Tuple[float, int]] = [
    (0.25 + r * r, 1) for r in MODULAR_R
]


@dataclass(frozen=True)
class Surface:
    name: str
    operator_label: str  # "Delta" (compact) or "Delta_cusp" (cuspidal restriction)
    spectrum: List[Tuple[float, int]]  # (eigenvalue, multiplicity)
    spectral_scope: str
    advice_source: str
    completeness_note: str
    # Conservative counting-density bound rho >= dN/dlambda for
    # lambda > cutoff (documented assumption: 2x Weyl main term Area/4pi).
    tail_density_bound: float


@dataclass(frozen=True)
class Window:
    surface: Surface
    interval: Tuple[float, float]
    nu: float
    label: str


@dataclass(frozen=True)
class Candidate:
    name: str
    role: str
    control_class: str  # positive | radial_miscandidate | mathematical_negative | bookkeeping_ablation
    # Either a multiplier h(lambda) (function of Delta via functional
    # calculus; commutes with Delta exactly, e_comm = 0 BY DEFINITION of
    # this candidate class) ...
    multiplier: Optional[Callable[[float], float]] = None
    # ... or explicit squared masses per spectrum index (adversarial
    # spectral distributions / structure breakers; finite by construction).
    masses: Optional[Dict[int, float]] = None
    # Optional permutation of the spectrum the multiplier is evaluated on
    # (bookkeeping ablation: randomized lengths / wrong spectral records).
    permutation: Optional[Sequence[int]] = None
    # Optional indices to exclude (restriction to an invariant subspace
    # such as L^2_0; legitimate because the subspace is Delta-invariant).
    exclude: Optional[Sequence[int]] = None


# Gauss-Bonnet: Bolza has genus 2, Area = 4pi  -> Weyl main term Area/4pi = 1.
BOLZA = Surface(
    name="Bolza (compact, genus 2)",
    operator_label="Delta",
    spectrum=BOLZA_SPECTRUM,
    spectral_scope="full discrete spectrum (compact surface; lemma assumption satisfied)",
    advice_source="Strohmaier-Uski 2013 certified eigenvalues (arXiv:1110.2150)",
    completeness_note=(
        "Strohmaier-Uski compute ALL eigenvalues in a prescribed interval with "
        "rigorous error bounds; the list is complete below the cutoff "
        "lambda_max = 30.833..."
    ),
    tail_density_bound=2.0,  # 2 x (Area/4pi) = 2 x 1.0
)

# PSL(2,Z): Area = pi/3 -> cuspidal Weyl main term (pi/3)/(4pi) = 1/12.
MODULAR = Surface(
    name="Modular surface PSL(2,Z), level 1 (cuspidal restriction)",
    operator_label="Delta_cusp",
    spectrum=MODULAR_SPECTRUM,
    spectral_scope=(
        "sigma(Delta_cusp) on L^2_cusp ONLY (non-compact surface; the continuous "
        "Eisenstein spectrum [1/4, oo) of the full Delta is outside ledger scope; "
        "all gap statements refer to Delta_cusp)"
    ),
    advice_source="LMFDB rigorous Maass database / Seymour-Howell (arXiv:2201.08760)",
    completeness_note=(
        "LMFDB rigorous Maass-form database lists the cuspidal spectrum for level 1 "
        "completely in the covered r-range; list complete below lambda_max = 330.8..."
    ),
    tail_density_bound=2.0 / 12.0,  # 2 x (Area/4pi) = 2/12
)


def gaussian_multiplier(center: float, sigma: float) -> Callable[[float], float]:
    def h(lam: float) -> float:
        arg = ((lam - center) / sigma) ** 2
        return math.exp(-arg) if arg < 700.0 else 0.0

    return h


def heat_multiplier(t: float) -> Callable[[float], float]:
    def h(lam: float) -> float:
        arg = t * lam
        return math.exp(-arg) if arg < 700.0 else 0.0

    return h


def two_mode_rayleigh_masses(
    spectrum: List[Tuple[float, int]], i_low: int, i_high: int, nu: float
) -> Dict[int, float]:
    """Squared masses on two outside modes with Rayleigh quotient exactly nu."""
    lam_low = spectrum[i_low][0]
    lam_high = spectrum[i_high][0]
    w = (lam_high - nu) / (lam_high - lam_low)
    return {i_low: w, i_high: 1.0 - w}


def simpson_tail(
    f: Callable[[float], float], a: float, b: float, n: int = 20001
) -> float:
    """Plain Simpson rule on [a, b] with odd n nodes."""
    if n % 2 == 0:
        n += 1
    step = (b - a) / (n - 1)
    total = f(a) + f(b)
    for i in range(1, n - 1):
        x = a + i * step
        total += f(x) * (4.0 if i % 2 == 1 else 2.0)
    return total * step / 3.0


def tail_budgets(
    h: Callable[[float], float], cutoff: float, nu: float, rho: float
) -> Tuple[float, float]:
    """Conservative tail budgets T0, T2 for a multiplier candidate.

    Uses the documented counting-density bound rho >= dN/dlambda above the
    cutoff.  Integration range [cutoff, cutoff+600] suffices: all used
    multipliers underflow to 0 long before the upper end.
    """
    upper = cutoff + 600.0
    t0 = simpson_tail(lambda x: rho * h(x) ** 2, cutoff, upper)
    t2 = simpson_tail(lambda x: rho * (h(x) ** 2) * (x - nu) ** 2, cutoff, upper)
    return t0, t2


def build_windows_and_candidates() -> List[Tuple[Window, List[Candidate]]]:
    out: List[Tuple[Window, List[Candidate]]] = []

    # ---- Bolza, window around lambda_1 (multiplicity 3) ------------------
    nu1 = BOLZA_SPECTRUM[1][0]
    win1 = Window(BOLZA, (3.5, 4.2), nu1, "bolza_lambda1_mult3")
    out.append(
        (
            win1,
            [
                Candidate(
                    "pure_window_eigenspace",
                    "positive control: aggregated spectral mass on the full lambda_1 eigenspace (multiplicity 3)",
                    "positive",
                    masses={1: 1.0},
                ),
                Candidate(
                    "gaussian_window_multiplier",
                    "radial multiplier h(Delta), Gaussian window sigma=0.6",
                    "positive",
                    multiplier=gaussian_multiplier(nu1, 0.6),
                ),
                Candidate(
                    "heat_kernel_multiplier_t1",
                    "genuine h(Delta), untuned heat kernel t=1 (concentrates on the constants)",
                    "radial_miscandidate",
                    multiplier=heat_multiplier(1.0),
                ),
                Candidate(
                    "heat_kernel_t1_on_L2_0",
                    "heat multiplier exp(-t Delta) restricted to the Delta-invariant subspace L^2_0, t=1",
                    "positive",
                    multiplier=heat_multiplier(1.0),
                    exclude=[0],
                ),
                Candidate(
                    "wide_gaussian_sigma5",
                    "genuine h(Delta), window-shaped but far too wide (sigma=5): mass spills over neighbouring eigenvalues",
                    "radial_miscandidate",
                    multiplier=gaussian_multiplier(nu1, 5.0),
                ),
                Candidate(
                    "miscentered_bandpass_at_10",
                    "genuine h(Delta), Gaussian bandpass centred at 10 (no eigenvalue near centre, window missed)",
                    "radial_miscandidate",
                    multiplier=gaussian_multiplier(10.0, 1.5),
                ),
                Candidate(
                    "same_rayleigh_wrong_tail",
                    "adversarial spectral distribution: Rayleigh = nu from lambda_0 and lambda_2 only",
                    "mathematical_negative",
                    masses=two_mode_rayleigh_masses(BOLZA_SPECTRUM, 0, 2, nu1),
                ),
                Candidate(
                    "false_center_control",
                    "adversarial control: clean lambda_2 eigenspace against the lambda_1 window",
                    "mathematical_negative",
                    masses={2: 1.0},
                ),
                Candidate(
                    "broken_radial_commutation",
                    "bookkeeping ablation: Gaussian weights swapped between lambda_1 and lambda_4 (NOT a function of Delta)",
                    "bookkeeping_ablation",
                    multiplier=gaussian_multiplier(nu1, 0.6),
                    permutation=[0, 4, 2, 3, 1, 5, 6, 7, 8, 9, 10],
                ),
                Candidate(
                    "permuted_spectrum_control",
                    "bookkeeping ablation: multiplier evaluated on reversed spectrum (randomized-lengths analogue)",
                    "bookkeeping_ablation",
                    multiplier=gaussian_multiplier(nu1, 0.6),
                    permutation=list(range(len(BOLZA_SPECTRUM)))[::-1],
                ),
            ],
        )
    )

    # ---- Bolza, isolated window around lambda_8 (multiplicity 1) ---------
    nu8 = BOLZA_SPECTRUM[8][0]
    win8 = Window(BOLZA, (22.0, 24.0), nu8, "bolza_lambda8_isolated")
    out.append(
        (
            win8,
            [
                Candidate(
                    "pure_window_eigenspace",
                    "positive control: exact lambda_8 eigenvector (multiplicity 1)",
                    "positive",
                    masses={8: 1.0},
                ),
                Candidate(
                    "gaussian_window_multiplier",
                    "radial multiplier h(Delta), Gaussian window sigma=1.0",
                    "positive",
                    multiplier=gaussian_multiplier(nu8, 1.0),
                ),
                Candidate(
                    "same_rayleigh_wrong_tail",
                    "adversarial spectral distribution: Rayleigh = nu from lambda_7 and lambda_9 only",
                    "mathematical_negative",
                    masses=two_mode_rayleigh_masses(BOLZA_SPECTRUM, 7, 9, nu8),
                ),
            ],
        )
    )

    # ---- Modular surface (Delta_cusp), window around second cusp form ----
    # The second eigenvalue is bracketed by cuspidal modes on both sides,
    # so the same-Rayleigh negative control is genuinely constructible.
    nu_m = MODULAR_SPECTRUM[1][0]
    win_m = Window(MODULAR, (144.0, 153.0), nu_m, "modular_r2_cuspidal")
    out.append(
        (
            win_m,
            [
                Candidate(
                    "pure_window_eigenvector",
                    "positive control: exact second cusp form (Delta_cusp eigenvector)",
                    "positive",
                    masses={1: 1.0},
                ),
                Candidate(
                    "gaussian_window_multiplier",
                    "radial multiplier h(Delta_cusp), Gaussian window sigma=12",
                    "positive",
                    multiplier=gaussian_multiplier(nu_m, 12.0),
                ),
                Candidate(
                    "same_rayleigh_wrong_tail",
                    "adversarial spectral distribution: Rayleigh = nu from first and third cusp form only",
                    "mathematical_negative",
                    masses=two_mode_rayleigh_masses(MODULAR_SPECTRUM, 0, 2, nu_m),
                ),
            ],
        )
    )
    return out


def evaluate(window: Window, cand: Candidate) -> Dict[str, object]:
    spec = window.surface.spectrum
    lo, hi = window.interval

    def in_window(lam: float) -> bool:
        return lo <= lam <= hi

    # unnormalized squared mass per spectrum index (aggregated per
    # eigenvalue, multiplicity-weighted for multiplier candidates)
    raw: Dict[int, float] = {}
    tail_t0: Optional[float] = None
    tail_t2: Optional[float] = None
    if cand.multiplier is not None:
        perm = list(cand.permutation) if cand.permutation is not None else list(range(len(spec)))
        excluded = set(cand.exclude or [])
        for idx, (lam, mult) in enumerate(spec):
            if idx in excluded:
                continue
            lam_eval = spec[perm[idx]][0]
            h = cand.multiplier(lam_eval)
            raw[idx] = mult * h * h
        # Tail budgets only make sense for genuine (unpermuted) multipliers;
        # for bookkeeping ablations the construction itself is invalid and
        # the finite verdict already fails.
        if cand.permutation is None:
            cutoff = spec[-1][0]
            tail_t0, tail_t2 = tail_budgets(
                cand.multiplier, cutoff, window.nu, window.surface.tail_density_bound
            )
    else:
        raw = dict(cand.masses or {})

    known_total = sum(raw.values())
    if known_total <= 0:
        raise ValueError(f"zero mass for {cand.name}")

    nu = window.nu
    out_known = sum(m for i, m in raw.items() if not in_window(spec[i][0]))
    res2_known = sum(m * (spec[i][0] - nu) ** 2 for i, m in raw.items())
    rayleigh = sum(m * spec[i][0] for i, m in raw.items()) / known_total

    residual = math.sqrt(res2_known / known_total)
    off_mass = math.sqrt(out_known / known_total)
    outside = [lam for lam, _ in spec if not in_window(lam)]
    gap = min(abs(lam - nu) for lam in outside)

    # Bracketing check: unknown eigenvalues lie above the cutoff; they can
    # only shrink the gap if cutoff - nu < gap (never the case here, but
    # verified explicitly).
    cutoff = spec[-1][0]
    tail_gap_ok = (cutoff - nu) > gap

    ratio = residual / gap
    bound_holds = off_mass <= ratio + 1e-12  # finite-list sanity check only

    # Worst-case tail-adjusted metrics for genuine multipliers: the unknown
    # tail mass (<= T0) lies entirely outside the window and contributes
    # <= T2 to the squared residual.
    if tail_t0 is not None:
        adj_off = math.sqrt((out_known + tail_t0) / (known_total + tail_t0))
        adj_res = math.sqrt((res2_known + tail_t2) / known_total)
        adj_ratio = adj_res / gap
    else:
        adj_off = off_mass
        adj_res = residual
        adj_ratio = ratio

    def decide(r: float, om: float) -> str:
        if r <= PASS_RATIO and om <= PASS_RATIO:
            return "pass"
        if r <= AMBIGUOUS_RATIO and om <= AMBIGUOUS_RATIO:
            return "ambiguous"
        return "fail"

    return {
        "surface": window.surface.name,
        "operator": window.surface.operator_label,
        "spectral_scope": window.surface.spectral_scope,
        "completeness_note": window.surface.completeness_note,
        "window_label": window.label,
        "target_window": [lo, hi],
        "advice_source": window.surface.advice_source,
        "candidate": cand.name,
        "role": cand.role,
        "control_class": cand.control_class,
        "nu": nu,
        "rayleigh": rayleigh,
        "outer_gap": gap,
        "residual": residual,
        "certificate_ratio": ratio,
        "off_window_mass": off_mass,
        "tail_T0": tail_t0,
        "tail_T2": tail_t2,
        "tail_gap_ok": tail_gap_ok,
        "adjusted_residual": adj_res,
        "adjusted_ratio": adj_ratio,
        "adjusted_off_window_mass": adj_off,
        "bound_holds_finite_list": bound_holds,
        "decision_finite_proxy": decide(ratio, off_mass),
        "decision": decide(adj_ratio, adj_off),
    }


def fmt(value: object) -> str:
    if isinstance(value, float):
        return f"{value:.6f}" if (value == 0 or abs(value) >= 1e-6) else f"{value:.2e}"
    if value is None:
        return "-"
    return str(value)


def write_markdown(rows: List[Dict[str, object]], path: Path) -> None:
    lines = [
        "# Selberg C2c-Window-Ledger (Real Data, Rev. 2)",
        "",
        f"Datum: {DATE} (Rev. 2 nach Codex-Gegenreview, siehe",
        "`CODEX_REVIEW_C2C_REAL_LEDGER_2026-07-03.md`)",
        "",
        "Status: Real-Daten-Lauf des C2c-Fensterzertifikats als **Finite-Spectrum-Proxy",
        "mit explizitem Tail-Budget**. Kein Selberg-RH-Claim, keine Minimizer-",
        "Identifikation. Der Lauf ersetzt das Toy-Spektrum vom 2026-05-27 durch",
        "publizierte, unabhängig berechnete Spektren zweier Positivkontrollflächen.",
        "",
        "## Datenquellen (web-verifiziert 2026-07-03)",
        "",
        "- **Bolza-Fläche** (kompakt, Genus 2): erste zehn positive Laplace-Eigenwerte mit",
        "  Multiplizitäten, zertifizierte Berechnung Strohmaier–Uski 2013",
        "  (Comm. Math. Phys. 317, arXiv:1110.2150). Kompakt ⇒ Lemma-Annahme erfüllt.",
        "  Vollständigkeit: das Verfahren findet ALLE Eigenwerte im vorgegebenen Intervall;",
        "  die Liste ist unterhalb des Cutoffs λ_max ≈ 30.833 vollständig.",
        "- **Modulfläche PSL(2,Z), Level 1**: erste acht kuspidale Spektralparameter r aus der",
        "  rigorosen LMFDB-Maass-Datenbank (Seymour-Howell-Linie, arXiv:2201.08760,",
        "  arXiv:2502.01442); λ = 1/4 + r². **Operator-Scope:** Die Fläche ist NICHT kompakt;",
        "  alle Zeilen beziehen sich ausschließlich auf die selbstadjungierte Restriktion",
        "  **Δ_cusp auf L²_cusp** (Spektrum σ(Δ_cusp), Projektor P_I^cusp). Das kontinuierliche",
        "  Eisenstein-Spektrum [1/4, ∞) des vollen Δ liegt AUSSERHALB des Ledger-Scopes;",
        "  die Gap-Spalte darf NICHT als dist(ν, σ(Δ)\\I) der vollen Fläche gelesen werden.",
        "",
        "## Geprüfter Mechanismus (Finite-Spectrum-Proxy + Tail-Budget)",
        "",
        "`||(I-P_I) psi|| <= ||(D-nu) psi|| / dist(nu, sigma(D) \\ I)`  mit D = Delta bzw. Delta_cusp.",
        "",
        "Die Ungleichung wird auf der endlichen, laut Quelle unterhalb des Cutoffs",
        "vollständigen Eigenwertliste ausgewertet (`decision_finite_proxy`). Für",
        "Multiplikator-Kandidaten wird zusätzlich ein **konservatives Tail-Budget**",
        "eingepreist: mit der Zähldichte-Schranke ρ ≥ dN/dλ oberhalb des Cutoffs",
        "(dokumentierte Annahme: 2 × Weyl-Hauptterm Area/4π; Bolza ρ=2, Modulfläche",
        "kuspidal ρ=1/6) werden",
        "",
        "`T0 >= Σ_{λ>Λ} mult·|h(λ)|²`  und  `T2 >= Σ_{λ>Λ} mult·|h(λ)|²·(λ-ν)²`",
        "",
        "numerisch integriert und Worst-Case-korrigierte Werte berichtet",
        "(`adjusted_*`; Tail-Masse vollständig außerhalb des Fensters angesetzt).",
        "Die Spalte `bound_holds_finite_list` ist auf der endlichen Liste nahezu",
        "tautologisch und dient nur als Implementations-Sanity-Check; das Urteil",
        "fällt über Residual/Gap und Off-Window-Masse (tail-adjustiert).",
        "",
        f"**Audit-Schwellen** (Konvention, NICHT Teil des Spektralsatzes): pass ≤ {PASS_RATIO},",
        f"ambiguous ≤ {AMBIGUOUS_RATIO}.",
        "",
        "**e_comm = 0 (konditional):** Für die hier modellierten Kandidaten, die per",
        "Definition Funktionen h(D) im Funktionalkalkül sind, kommutiert h(D) exakt mit D;",
        "der Kommutator-Leakage-Term ist für DIESE Klasse strukturell 0. Das Ledger prüft",
        "NICHT, dass ein realer Selberg-/Trace-Formula-/Connes-Minimizer-Kandidat in diese",
        "Klasse fällt — genau das bleibt der offene C2c-Kern. Die Restriktion von",
        "exp(-tΔ) auf L²₀ ist zulässig, weil L²₀ Δ-invariant ist.",
        "",
        "## Kontrollklassen",
        "",
        "- `positive` — Positivkontrolle (echte h(D)-Fensterkandidaten bzw. Eigenvektoren).",
        "- `radial_miscandidate` — echte Funktionen h(D) (kommutierend!), aber falsch",
        "  zentriert/getunt: Scheitern INNERHALB der theoremrelevanten Kandidatenklasse.",
        "- `mathematical_negative` — adversarielle Spektralverteilungen (handgesetzte Massen),",
        "  z. B. die Pflichtklasse `same_rayleigh_wrong_tail`.",
        "- `bookkeeping_ablation` — Strukturbrecher (keine Funktionen von D / permutierte",
        "  Spektralbuchhaltung): testen Implementations-Robustheit, nicht die Kernaussage.",
        "",
        "## Ergebnisse",
        "",
        "| Fenster | Kandidat | Klasse | Rayleigh | Res/Gap | Off-Win | T0 | Res/Gap adj. | Off-Win adj. | Decision |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        lines.append(
            "| {w} | {c} | {cl} | {ray} | {ratio} | {off} | {t0} | {aratio} | {aoff} | **{d}** |".format(
                w=row["window_label"],
                c=row["candidate"],
                cl=row["control_class"],
                ray=fmt(row["rayleigh"]),
                ratio=fmt(row["certificate_ratio"]),
                off=fmt(row["off_window_mass"]),
                t0=fmt(row["tail_T0"]),
                aratio=fmt(row["adjusted_ratio"]),
                aoff=fmt(row["adjusted_off_window_mass"]),
                d=row["decision"],
            )
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- Die Positivkontrollen bestehen auch NACH Tail-Adjustierung: die Tail-Budgets",
            "  der getunten Gauss-Fenster sind numerisch vernachlässigbar klein, die des",
            "  Heat-Kernels auf L²₀ bleiben unterhalb der Audit-Schwelle. Der Casimir-",
            "  verankerte Kanal reproduziert den bekannten positiven Fall (Kalibrierungsziel).",
            "- `heat_kernel_multiplier_t1` (ungetunt) FÄLLT DURCH: der rohe Heat-Kernel",
            "  konzentriert auf den konstanten Grundzustand λ0=0 statt auf das Zielfenster.",
            "- Die neuen echten radialen Fehlkandidaten (`wide_gaussian_sigma5`,",
            "  `miscentered_bandpass_at_10`) zeigen Scheitern INNERHALB der h(D)-Klasse:",
            "  Kommutation allein rettet kein falsch gebautes Fenster.",
            "- `same_rayleigh_wrong_tail` scheitert auch mit realen Daten trotz exaktem",
            "  Rayleigh-Schwerpunkt: die Pflicht-Negativklasse aus dem Toy-Lauf trägt real.",
            "- `broken_radial_commutation` / `permuted_spectrum_control` sind als",
            "  Bookkeeping-Ablationen gekennzeichnet: sie falsifizieren nicht die Kernaussage",
            "  über echte Multiplikatoren, sondern belegen, dass das Ledger Strukturbruch",
            "  erkennt.",
            "- Multiplizitäten: aggregierte Spektralmasse pro Eigenwert, multiplizitäts-",
            "  gewichtet bei Multiplikatoren (Bolza λ1 hat Multiplizität 3; P_I projiziert",
            "  auf den ganzen Eigenraum; eine Verteilung über einzelne Eigenrichtungen wird",
            "  nicht geprüft und ist für die Normgrößen unerheblich).",
            "- Bracketing: für jedes Fenster gilt cutoff − ν > gap (`tail_gap_ok`), d. h.",
            "  unbekannte Eigenwerte oberhalb des Cutoffs können den Gap nicht verkleinern.",
            "",
            "## Grenzen (aus dem Codex-Review übernommen)",
            "",
            "1. Das Ledger ist ein Finite-Spectrum-Proxy mit konservativem Tail-Budget,",
            "   kein abgeschlossenes Zertifikat: die Zähldichte-Schranke (2×Weyl) ist eine",
            "   dokumentierte Annahme, kein bewiesener Restterm-Bound.",
            "2. Modulflächen-Zeilen gelten nur für Δ_cusp auf L²_cusp.",
            "3. e_comm=0 ist eine Klassenaussage (Funktionalkalkül), keine geprüfte",
            "   Eigenschaft eines externen Minimizer-Kandidaten.",
            "4. pass/ambiguous sind Audit-Konventionen.",
            "",
            "## Konsequenz für C2c",
            "",
            "Das Fensterzertifikat ist mit realen Spektraldaten und Tail-Budget",
            "operationalisiert (TODO-Punkt 'Toy-Spektrum durch echte Laplace-Fenster",
            "ersetzen': erledigt). OFFEN bleibt der eigentliche C2c-Kern: die Identifikation",
            "des *Connes-Minimizers* mit dem Selberg-Spektralfaktor — dafür müsste der",
            "Minimizer selbst als Vektor konstruiert, als h(Δ)-Klasse nachgewiesen und durch",
            "dieses Ledger geschickt werden. Bis dahin bleibt C2c Companion-/Audit-Status;",
            "kein Claim-Upgrade.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    project_root = Path(__file__).resolve().parents[1]
    result_dir = project_root / "_results"
    result_dir.mkdir(exist_ok=True)
    rows: List[Dict[str, object]] = []
    for window, cands in build_windows_and_candidates():
        for cand in cands:
            rows.append(evaluate(window, cand))

    json_path = result_dir / f"selberg_c2c_window_ledger_real_{DATE}.json"
    md_path = result_dir / f"SELBERG_C2C_WINDOW_LEDGER_REAL_{DATE}.md"
    json_path.write_text(json.dumps(rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(rows, md_path)

    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")
    for row in rows:
        print(
            f"{row['window_label']:>24} | {row['candidate']:>28} | {row['control_class']:>22} | "
            f"adj_ratio={row['adjusted_ratio']:.4f} adj_off={row['adjusted_off_window_mass']:.4f} "
            f"-> {row['decision']}"
        )
    problems = [
        f"{row['window_label']}/{row['candidate']}"
        for row in rows
        if not row["bound_holds_finite_list"] or not row["tail_gap_ok"]
    ]
    if problems:
        raise SystemExit(f"sanity check failed for: {', '.join(problems)}")


if __name__ == "__main__":
    main()
