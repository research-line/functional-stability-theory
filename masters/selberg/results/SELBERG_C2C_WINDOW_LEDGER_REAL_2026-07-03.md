# Selberg C2c-Window-Ledger (Real Data, Rev. 2)

Datum: 2026-07-03 (Rev. 2 nach Codex-Gegenreview, siehe
`CODEX_REVIEW_C2C_REAL_LEDGER_2026-07-03.md`)

Status: Real-Daten-Lauf des C2c-Fensterzertifikats als **Finite-Spectrum-Proxy
mit explizitem Tail-Budget**. Kein Selberg-RH-Claim, keine Minimizer-
Identifikation. Der Lauf ersetzt das Toy-Spektrum vom 2026-05-27 durch
publizierte, unabhängig berechnete Spektren zweier Positivkontrollflächen.

## Datenquellen (web-verifiziert 2026-07-03)

- **Bolza-Fläche** (kompakt, Genus 2): erste zehn positive Laplace-Eigenwerte mit
  Multiplizitäten, zertifizierte Berechnung Strohmaier–Uski 2013
  (Comm. Math. Phys. 317, arXiv:1110.2150). Kompakt ⇒ Lemma-Annahme erfüllt.
  Vollständigkeit: das Verfahren findet ALLE Eigenwerte im vorgegebenen Intervall;
  die Liste ist unterhalb des Cutoffs λ_max ≈ 30.833 vollständig.
- **Modulfläche PSL(2,Z), Level 1**: erste acht kuspidale Spektralparameter r aus der
  rigorosen LMFDB-Maass-Datenbank (Seymour-Howell-Linie, arXiv:2201.08760,
  arXiv:2502.01442); λ = 1/4 + r². **Operator-Scope:** Die Fläche ist NICHT kompakt;
  alle Zeilen beziehen sich ausschließlich auf die selbstadjungierte Restriktion
  **Δ_cusp auf L²_cusp** (Spektrum σ(Δ_cusp), Projektor P_I^cusp). Das kontinuierliche
  Eisenstein-Spektrum [1/4, ∞) des vollen Δ liegt AUSSERHALB des Ledger-Scopes;
  die Gap-Spalte darf NICHT als dist(ν, σ(Δ)\I) der vollen Fläche gelesen werden.

## Geprüfter Mechanismus (Finite-Spectrum-Proxy + Tail-Budget)

`||(I-P_I) psi|| <= ||(D-nu) psi|| / dist(nu, sigma(D) \ I)`  mit D = Delta bzw. Delta_cusp.

Die Ungleichung wird auf der endlichen, laut Quelle unterhalb des Cutoffs
vollständigen Eigenwertliste ausgewertet (`decision_finite_proxy`). Für
Multiplikator-Kandidaten wird zusätzlich ein **konservatives Tail-Budget**
eingepreist: mit der Zähldichte-Schranke ρ ≥ dN/dλ oberhalb des Cutoffs
(dokumentierte Annahme: 2 × Weyl-Hauptterm Area/4π; Bolza ρ=2, Modulfläche
kuspidal ρ=1/6) werden

`T0 >= Σ_{λ>Λ} mult·|h(λ)|²`  und  `T2 >= Σ_{λ>Λ} mult·|h(λ)|²·(λ-ν)²`

numerisch integriert und Worst-Case-korrigierte Werte berichtet
(`adjusted_*`; Tail-Masse vollständig außerhalb des Fensters angesetzt).
Die Spalte `bound_holds_finite_list` ist auf der endlichen Liste nahezu
tautologisch und dient nur als Implementations-Sanity-Check; das Urteil
fällt über Residual/Gap und Off-Window-Masse (tail-adjustiert).

**Audit-Schwellen** (Konvention, NICHT Teil des Spektralsatzes): pass ≤ 0.25,
ambiguous ≤ 0.5.

**e_comm = 0 (konditional):** Für die hier modellierten Kandidaten, die per
Definition Funktionen h(D) im Funktionalkalkül sind, kommutiert h(D) exakt mit D;
der Kommutator-Leakage-Term ist für DIESE Klasse strukturell 0. Das Ledger prüft
NICHT, dass ein realer Selberg-/Trace-Formula-/Connes-Minimizer-Kandidat in diese
Klasse fällt — genau das bleibt der offene C2c-Kern. Die Restriktion von
exp(-tΔ) auf L²₀ ist zulässig, weil L²₀ Δ-invariant ist.

## Kontrollklassen

- `positive` — Positivkontrolle (echte h(D)-Fensterkandidaten bzw. Eigenvektoren).
- `radial_miscandidate` — echte Funktionen h(D) (kommutierend!), aber falsch
  zentriert/getunt: Scheitern INNERHALB der theoremrelevanten Kandidatenklasse.
- `mathematical_negative` — adversarielle Spektralverteilungen (handgesetzte Massen),
  z. B. die Pflichtklasse `same_rayleigh_wrong_tail`.
- `bookkeeping_ablation` — Strukturbrecher (keine Funktionen von D / permutierte
  Spektralbuchhaltung): testen Implementations-Robustheit, nicht die Kernaussage.

## Ergebnisse

| Fenster | Kandidat | Klasse | Rayleigh | Res/Gap | Off-Win | T0 | Res/Gap adj. | Off-Win adj. | Decision |
|---|---|---|---:|---:|---:|---:|---:|---:|---|
| bolza_lambda1_mult3 | pure_window_eigenspace | positive | 3.838887 | 0.000000 | 0.000000 | - | 0.000000 | 0.000000 | **pass** |
| bolza_lambda1_mult3 | gaussian_window_multiplier | positive | 3.838893 | 0.001971 | 0.001971 | 0.000000 | 0.001971 | 0.001971 | **pass** |
| bolza_lambda1_mult3 | heat_kernel_multiplier_t1 | radial_miscandidate | 0.005804 | 2.532543 | 0.999306 | 1.65e-27 | 2.532543 | 0.999306 | **fail** |
| bolza_lambda1_mult3 | heat_kernel_t1_on_L2_0 | positive | 3.931009 | 0.247654 | 0.246255 | 1.65e-27 | 0.247654 | 0.246255 | **pass** |
| bolza_lambda1_mult3 | wide_gaussian_sigma5 | radial_miscandidate | 4.650224 | 1.123128 | 0.758298 | 2.21e-26 | 1.123128 | 0.758298 | **fail** |
| bolza_lambda1_mult3 | miscentered_bandpass_at_10 | radial_miscandidate | 8.249555 | 2.911882 | 1.000000 | 1.54e-169 | 2.911882 | 1.000000 | **fail** |
| bolza_lambda1_mult3 | same_rayleigh_wrong_tail | mathematical_negative | 3.838887 | 1.591979 | 1.000000 | - | 1.591979 | 1.000000 | **fail** |
| bolza_lambda1_mult3 | false_center_control | mathematical_negative | 5.353601 | 1.000000 | 1.000000 | - | 1.000000 | 1.000000 | **fail** |
| bolza_lambda1_mult3 | broken_radial_commutation | bookkeeping_ablation | 14.726189 | 7.187702 | 1.000000 | - | 7.187702 | 1.000000 | **fail** |
| bolza_lambda1_mult3 | permuted_spectrum_control | bookkeeping_ablation | 28.079601 | 16.003492 | 1.000000 | - | 16.003492 | 1.000000 | **fail** |
| bolza_lambda8_isolated | pure_window_eigenspace | positive | 23.078558 | 0.000000 | 0.000000 | - | 0.000000 | 0.000000 | **pass** |
| bolza_lambda8_isolated | gaussian_window_multiplier | positive | 23.078537 | 0.002869 | 0.002869 | 3.80e-54 | 0.002869 | 0.002869 | **pass** |
| bolza_lambda8_isolated | same_rayleigh_wrong_tail | mathematical_negative | 23.078558 | 1.398044 | 1.000000 | - | 1.398044 | 1.000000 | **fail** |
| modular_r2_cuspidal | pure_window_eigenvector | positive | 148.432131 | 0.000000 | 0.000000 | - | 0.000000 | 0.000000 | **pass** |
| modular_r2_cuspidal | gaussian_window_multiplier | positive | 148.432131 | 0.000006 | 0.000006 | 8.28e-203 | 0.000006 | 0.000006 | **pass** |
| modular_r2_cuspidal | same_rayleigh_wrong_tail | mathematical_negative | 148.432131 | 1.172134 | 1.000000 | - | 1.172134 | 1.000000 | **fail** |

## Interpretation

- Die Positivkontrollen bestehen auch NACH Tail-Adjustierung: die Tail-Budgets
  der getunten Gauss-Fenster sind numerisch vernachlässigbar klein, die des
  Heat-Kernels auf L²₀ bleiben unterhalb der Audit-Schwelle. Der Casimir-
  verankerte Kanal reproduziert den bekannten positiven Fall (Kalibrierungsziel).
- `heat_kernel_multiplier_t1` (ungetunt) FÄLLT DURCH: der rohe Heat-Kernel
  konzentriert auf den konstanten Grundzustand λ0=0 statt auf das Zielfenster.
- Die neuen echten radialen Fehlkandidaten (`wide_gaussian_sigma5`,
  `miscentered_bandpass_at_10`) zeigen Scheitern INNERHALB der h(D)-Klasse:
  Kommutation allein rettet kein falsch gebautes Fenster.
- `same_rayleigh_wrong_tail` scheitert auch mit realen Daten trotz exaktem
  Rayleigh-Schwerpunkt: die Pflicht-Negativklasse aus dem Toy-Lauf trägt real.
- `broken_radial_commutation` / `permuted_spectrum_control` sind als
  Bookkeeping-Ablationen gekennzeichnet: sie falsifizieren nicht die Kernaussage
  über echte Multiplikatoren, sondern belegen, dass das Ledger Strukturbruch
  erkennt.
- Multiplizitäten: aggregierte Spektralmasse pro Eigenwert, multiplizitäts-
  gewichtet bei Multiplikatoren (Bolza λ1 hat Multiplizität 3; P_I projiziert
  auf den ganzen Eigenraum; eine Verteilung über einzelne Eigenrichtungen wird
  nicht geprüft und ist für die Normgrößen unerheblich).
- Bracketing: für jedes Fenster gilt cutoff − ν > gap (`tail_gap_ok`), d. h.
  unbekannte Eigenwerte oberhalb des Cutoffs können den Gap nicht verkleinern.

## Grenzen (aus dem Codex-Review übernommen)

1. Das Ledger ist ein Finite-Spectrum-Proxy mit konservativem Tail-Budget,
   kein abgeschlossenes Zertifikat: die Zähldichte-Schranke (2×Weyl) ist eine
   dokumentierte Annahme, kein bewiesener Restterm-Bound.
2. Modulflächen-Zeilen gelten nur für Δ_cusp auf L²_cusp.
3. e_comm=0 ist eine Klassenaussage (Funktionalkalkül), keine geprüfte
   Eigenschaft eines externen Minimizer-Kandidaten.
4. pass/ambiguous sind Audit-Konventionen.

## Konsequenz für C2c

Das Fensterzertifikat ist mit realen Spektraldaten und Tail-Budget
operationalisiert (TODO-Punkt 'Toy-Spektrum durch echte Laplace-Fenster
ersetzen': erledigt). OFFEN bleibt der eigentliche C2c-Kern: die Identifikation
des *Connes-Minimizers* mit dem Selberg-Spektralfaktor — dafür müsste der
Minimizer selbst als Vektor konstruiert, als h(Δ)-Klasse nachgewiesen und durch
dieses Ledger geschickt werden. Bis dahin bleibt C2c Companion-/Audit-Status;
kein Claim-Upgrade.
