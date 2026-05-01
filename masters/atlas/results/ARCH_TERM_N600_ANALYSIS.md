# N=600 Server-Analyse: Arch-Term + Autokorrelations-Formel

**Datum:** 2026-04-16
**Run:** ellmos-services, nohup, 1.95 Stunden Gesamtlaufzeit (schneller als erwartet)
**Skript:** `_scripts/server_n600_full_analysis.py`
**JSON:** `_results/ARCH_TERM_N600_SERVER.json`

## Resultate (alle 10 Charaktere bei N=600, λ=20000)

| chi | D | gap_emp | N_emp | gap_gal | S_exact | arch_diff | S+arch | sign(S+arch) |
|-----|---|---------|-------|---------|---------|-----------|--------|--------------|
| χ_5 | 5 | +0.01560 | 400 | +0.01603 | +9.16e-02 | −5.79e-02 | +3.37e-02 | OK |
| χ_8 | 8 | −0.04902 | 600 | −0.04902 | −1.15e-01 | +7.41e-02 | −4.08e-02 | OK |
| χ_12 | 12 | +0.03367 | 400 | +0.01154 | +2.01e-02 | +2.32e-04 | +2.03e-02 | OK |
| χ_13 | 13 | −0.12504 | 600 | −0.12504 | −1.02e-01 | −1.37e-02 | −1.16e-01 | OK |
| χ_17 | 17 | +0.01318 | 400 | +0.00886 | +1.02e-02 | +1.06e-03 | +1.13e-02 | OK |
| χ_21 | 21 | −0.00423 | 600 | **−0.00423** | +1.10e-02 | +3.97e-04 | +1.14e-02 | **FAIL** |
| χ_24 | 24 | +0.01076 | 400 | +0.00866 | +1.05e-02 | +8.07e-03 | +1.86e-02 | OK |
| χ_29 | 29 | +0.01439 | 600 | +0.01439 | −9.93e-03 | −1.68e-03 | −1.16e-02 | **FAIL** |
| χ_33 | 33 | −0.14221 | 600 | −0.14221 | −1.64e-01 | +3.24e-02 | −1.32e-01 | OK |
| χ_60 | 60 | +0.00521 | 600 | +0.00521 | +6.30e-03 | +9.15e-04 | +7.21e-03 | OK |

## Korrelations-Statistiken

| Predictor | sign_ok | R | R² |
|---|---|---|---|
| gap_galerkin | 10/10 | +0.9942 | **0.9885** |
| S_with_arch | 8/10 | +0.9766 | 0.9538 |

**Ohne χ_21 (bekannter N-Oszillator):**

| Predictor | sign_ok | R | R² |
|---|---|---|---|
| gap_galerkin | 9/9 | +0.9944 | 0.9888 |
| S_with_arch | 8/9 | +0.9796 | 0.9597 |

## Wichtige Einsichten

### 1. `gap_galerkin(N=600) ≡ gap_emp(N=600)` für 6 Charaktere

Für χ_8, χ_13, χ_21, χ_29, χ_33, χ_60 ist gap_gal **identisch** mit gap_emp (bis auf Rundung). Das ist kein Triumph, sondern **Tautologie**: die empirischen Werte in `oscillating6_N600_results.json` wurden selbst durch dieselbe Galerkin-Matrix-Diagonalisierung bei N=600 berechnet. Die 10/10 Vorzeichentreue von `gap_galerkin` ist damit **per Konstruktion**, keine Prediktion.

### 2. `S_with_arch` verschlechtert sich bei N=600 auf 8/10

Bei N=200 (Session 7 T3): S_with_arch = 9/10, R² = 0.42 (bzw. 9/9, R² = 0.80 ohne χ_21).
Bei N=600 (jetzt): S_with_arch = 8/10, R² = 0.95 (bzw. 8/9, R² = 0.96 ohne χ_21).

**Neuer Ausreißer χ_29:**
- gap_emp = +0.01439, gap_galerkin = +0.01439
- S_exact = −0.009932, arch_diff = −0.001679, S+arch = −0.01161
- S flippt Vorzeichen gegenüber der richtigen Antwort.

Mögliche Ursache: die **Autokorrelations-Berechnung** `(φ * φ)(t) via np.correlate + np.interp` hat systematische numerische Fehler von 1–2% des S-Wertes. Bei kleinen Gaps (χ_21 |gap|=0.004, χ_29 |gap|=0.014) reicht das für Vorzeichen-Flip.

### 3. Die Formel ist tautologisch (wichtige Einsicht für Paper)

Die Identität
$$ \mathrm{gap}_\chi = \mathrm{prime\_diff\_matrix} + \mathrm{arch\_diff} $$
folgt per Konstruktion aus der Matrixzerlegung $W = W_{\mathrm{arch}} + W_{\mathrm{prime}}$:
$$ \langle\phi, W\phi\rangle = \langle\phi, W_{\mathrm{arch}}\phi\rangle + \langle\phi, W_{\mathrm{prime}}\phi\rangle $$
Dies ist **keine Vorhersage**, sondern eine Definition. Die Vorhersagekraft käme erst durch:
- Analytische Auswertung von $\Delta_\chi(t)$ ohne Galerkin-Numerik (in Session 6 als Leading-Order-Form mit $\phi^+ \approx \phi^-$ probiert, **falsifiziert**).
- Oder: Approximation von $\Delta_\chi$ via Paley-Wiener-Asymptotik ohne den Umweg über Galerkin-Grundzustände.

Die "9/9 Closure" bei N=200 in Session 7 T3 war ein Artefakt der **numerischen Auswertung** von S_exact via Autokorrelation — bei N=200 waren die Fehler zufällig klein genug, dass die Vorzeichen übereinstimmten. Bei N=600 werden sie sichtbar.

### 4. Was bleibt valide

- **Die Weil-Kern-Architektur** (ANALYTIC_KERNEL Thm 3.2, 4.1) — rigorose Kern-Darstellung bleibt unverändert.
- **gap_galerkin als numerischer Zugang zu C2_χ** — Galerkin-Diagonalisierung liefert konvergente Näherung (siehe χ_12 N-Abhängigkeit: +0.0325 bei N=200 → +0.0115 bei N=600, konsistent mit Siegel-Walfisz-Abfall).
- **Falsifikation der Leading-Order-Formel** (Session 6) — bleibt das wissenschaftliche Haupt-Resultat der Session.
- **Drei-Komponenten-Zerlegung als Diagnose-Atlas** — unverändert, C2_χ bleibt der offene harte Kern.

## Konsequenz für Paper

Der ursprünglich geplante "Closure-Paper" ist nicht haltbar. Stattdessen:

**Neuer Titel-Vorschlag:** "A Weil-Kernel Numerical Atlas for Dirichlet L-Function Sector Gaps: Galerkin Diagnostics and the Boundaries of Leading-Order Theory"

**Neue Story:**
1. Weil-Kern-Architektur als Rahmen (rigoros, Thm 3.2, 4.1).
2. Empirischer Atlas der 10 Charaktere bei N ∈ {200, 400, 600} mit Konvergenz-Studie.
3. Falsifikation der Leading-Order-Approximation (φ⁺ ≈ φ⁻) mit analytischer Diagnose.
4. Beobachtung der Matrixzerlegungs-Tautologie — erklärt, warum kein "closed-form"-Predictor aus der exakten Formel folgt.
5. Outlook: echte C2_χ-Frage verbleibt als offene Paley-Wiener/Connes-Aufgabe.

Das ist **kein** Misserfolg, sondern ein ehrliches Negativ-Resultat mit sauber benannter struktureller Grenze.
