# Woche-4-Validierung der Analytic-Gap-Formel — FALL C: Theorie scheitert

**Autor:** LG (Claude Opus 4.6 [1M], Session 6 Teil 5)
**Datum:** 2026-04-17
**Skripte:** `_scripts/analytic_gap_formula_test.py` (v1), `_scripts/analytic_gap_formula_test_v2.py` (v2)
**Daten-Outputs:** `_results/ANALYTIC_GAP_FORMULA_TEST.json`, `_results/ANALYTIC_GAP_FORMULA_TEST_V2.json`
**Kontext:** Validierung der Leading-Order-Gap-Formel aus ANALYTIC_KERNEL.md / ANALYTIC_GROUNDSTATE.md / ANALYTIC_GAP_FORMULA.md.

---

## 0. Zusammenfassung

**Verdikt: Fall C (Plan B aktivieren).** Die Leading-Order-Formel $\mathrm{gap}_\chi \approx 2\sum \chi(p)^m \log p/p^{m/2} \Phi_\infty(m\log p)$ **kann die empirischen Gap-Vorzeichen nicht vorhersagen** — unter keinem der getesteten Fenster-Modelle, inklusive des **empirisch berechneten Grundzustands**.

**Bester Fall:** 6/10 Vorzeichen korrekt, $R^2 = 0.15$ (Gauss um $t=0$, $\sigma = 0.5$). Erfolgskriterium des Handoffs §4.2 (sign_ok ≥ 9/10, $R^2 \geq 0.8$) **nicht erreicht**. Die getestete Form ist **nicht besser als Zufall** ($p = 0.5^{10} \cdot \binom{10}{6} = 0.21$ für reine Zufälligkeit, d.h. 6/10 liegt im normalen Schwankungsbereich).

**Konsequenz:** Die theoretische Architektur aus Woche 1–3 (Weil-Kern-Geometrie, Sektor-Differenz, Leading-Order-Prim-Autokorrelation) **erfasst nicht die dominante Struktur** des Dirichlet-Gaps. Wechsel zu Plan B (deskriptives Paper) ist indiziert.

---

## 1. Test-Aufbau

Empirische Gap-Werte bei $\lambda = 20000$ (primär N=600 wo verfügbar, sonst N=400):

| Charakter | $D$ | empirischer Gap | $N$ |
|---|:-:|---:|:-:|
| chi_5 | 5 | +0.01560 | 400 |
| chi_8 | 8 | -0.04902 | 600 |
| chi_12 | 12 | +0.03367 | 400 |
| chi_13 | 13 | -0.12504 | 600 |
| chi_17 | 17 | +0.01318 | 400 |
| chi_21 | 21 | -0.00423 | 600 |
| chi_24 | 24 | +0.01076 | 400 |
| chi_29 | 29 | +0.01439 | 600 |
| chi_33 | 33 | **-0.14221** | 600 |
| chi_60 | 60 | +0.00521 | 600 |

**Erwartung laut Theorie:** Vorhersage-Summe $\mathcal{S}_\chi(L)$ = Prim-Autokorrelation im Fenster $\Phi_\infty$ sollte in Vorzeichen und Relativ-Magnitude mit empirischen Gaps korrelieren.

---

## 2. Getestete Fenster-Modelle

### 2.1 Zusammenfassung

| # | Modell | Parameter | sign_ok | $R^2$ |
|---|---|---|:-:|:-:|
| 1 | Raw chi-Summe (keine Fenster-Gewichtung) | $p \leq 20000$ | 4/10 | 0.214 |
| 2 | Gauss zentriert $t=0$ | $\sigma = 0.5$ (v1 Best) | **6/10** | **0.146** |
| 3 | Gauss zentriert $t=L$ (Frontier) | $\sigma = 0.5$ | 4/10 | 0.000 |
| 4 | Gauss zentriert $t=L$ | $\sigma = 1.0$ | 6/10 | 0.047 |
| 5 | Gauss zentriert $t=L$ | $\sigma = 2.0$ | 6/10 | 0.070 |
| 6 | Gauss zentriert $t=L$ | $\sigma = 3.0$ | 6/10 | 0.070 |
| 7 | Rechteck $[0, T]$ | $T = L/4$ | 6/10 | 0.070 |
| 8 | Rechteck $[0, T]$ | $T = L/2$ | 4/10 | 0.008 |
| 9 | Rechteck $[0, T]$ | $T = L$ | 4/10 | 0.024 |
| 10 | Rechteck $[0, T]$ | $T = 3L/2$ | 6/10 | 0.023 |
| 11 | Rechteck $[0, T]$ | $T = 2L$ | 6/10 | 0.023 |
| 12 | Cosine² auf $[-T, T]$ | $T = L/2$ | 5/10 | 0.004 |
| 13 | Cosine² auf $[-T, T]$ | $T = L$ | 4/10 | 0.005 |
| 14 | Cosine² auf $[-T, T]$ | $T = 3L/2$ | 4/10 | 0.000 |
| 15 | Cosine² auf $[-T, T]$ | $T = 2L$ | 5/10 | 0.005 |
| **16** | **Empirischer Grundzustand (Galerkin N=200)** | Autokorrelation des echten $\phi_\chi^+$ | **4/10** | **0.022** |

### 2.2 Der kritische Befund: Empirischer Grundzustand

Test 16 ist das **schärfste** Argument: wir berechnen den **echten** Grundzustand $\phi_\chi^+$ durch Galerkin-Eigenwert-Diagonalisierung (N=200, analog zum Verfahren in `all10_high_N_server.py`), bilden seine Autokorrelation $(\phi * \phi)(t)$ numerisch, und werten damit die Leading-Order-Formel aus.

**Resultat:**
- Alle 10 Charaktere geben $\mathcal{S}_\chi \in [-7.7, -6.5]$ — **alle negativ, nahezu identisch**.
- Sign match mit empirischem Gap: nur 4/10 (die 4 Charaktere mit negativem empirischen Gap: chi_8, chi_13, chi_21, chi_33).
- $R^2 = 0.022$: praktisch **keine Korrelation**.

**Interpretation.** Der echte Grundzustand $\phi_\chi^+$ variiert extrem wenig zwischen Charakteren (konsistent mit Prop. 5.1 aus Woche 2: "same archimedean ground state at leading order"). Aber das heißt: die Sektor-Differenz-Formel reduziert auf eine nahezu charakter-**universelle** Summe — sie unterscheidet die Gap-Vorzeichen **nicht**.

Das ist ein **definitiver Hinweis**, dass die Leading-Order-Näherung $\phi_\chi^+ \approx \phi_\chi^- \approx \phi_\infty$ falsch ist (oder zu grob): **die Differenz** der Grundzustände zwischen Sektoren ist essentiell und nicht vernachlässigbar.

---

## 3. Was schief lief

### 3.1 Prop. 5.1 (Woche 2) war zu optimistisch

Die zentrale Approximation in Woche 2 war:

> *$\phi_\chi^+$ und $\phi_\chi^-$ haben zu führender Ordnung die gleiche archimedisch-dominierte Form. Die Differenz beginnt erst bei subleading-Korrekturen aus dem Prim-Term.*

**Das stimmt nicht.** Die empirische Analyse zeigt:

- Die Grundzustände sind zwischen **Sektoren** (cos vs. sin) unterschiedlich genug, dass ihre **Differenz** das Gap-Vorzeichen treibt.
- Die charakter-spezifische Signatur des Gaps sitzt in **$\phi_\chi^+ - \phi_\chi^-$**, nicht in $\phi_\chi^+$ (oder $\phi_\chi^-$) allein.
- Die Leading-Order-Corollary 4.4 (ANALYTIC_KERNEL.md) war eine Näherung, die exakt diese Differenz eliminiert hat — damit hat sie die essentielle Information weggeworfen.

### 3.2 Der Weil-Kern-Ansatz ist strukturell nicht feinkörnig genug

Die zugrundeliegende Architektur von Woche 1 (Kern-Darstellung via Konvolution $\kappa_\chi$) ist **mathematisch korrekt**. Aber die Zerlegung in "gemeinsamen Teil" + "Sektor-Differenz" reduziert die Gap-Formel auf eine Form, bei der das Vorzeichen durch Interferenz zwischen verschiedenen Sektoren erst bei höherer Ordnung in $1/L$ verständlich wird.

Mit anderen Worten: $\mathrm{gap}_\chi = \langle\phi^-|K^-|\phi^-\rangle - \langle\phi^+|K^+|\phi^+\rangle$ lässt sich **nicht** sauber in einen führenden Term separieren, dessen Vorzeichen-Struktur durch die Leading-Order-Analyse erfasst ist. Die Näherungsformel Cor. 4.4 (die $\phi^+ \approx \phi^-$ ansetzt) wirft die interessante Struktur weg.

### 3.3 Die "Naive" Raw-Summe ist überraschend gut ($R^2 = 0.21$)

Interessant: die **ungewichtete** Summe $\sum_{p \leq \lambda} \chi(p) \log p/\sqrt p$ hat das höchste $R^2 = 0.214$ aller Tests (trotz nur 4/10 Vorzeichen). Das legt nahe, dass der Raw-Summen-Trend mit den Gaps einen Zusammenhang hat, aber dass weder Magnitude noch Vorzeichen durch sie allein erfasst werden.

Dies ist konsistent mit den **Session-5-Regressions-Ergebnissen** (Handoff §0 Punkt 1: multivariates $R^2 = 0.42$ mit mehreren Prädiktoren). Die Raw-Summe ist einer dieser Prädiktoren, aber auch sie reicht nicht.

---

## 4. Abbruchkriterium greift

Handoff §4.3 und ANALYTIC_GAP_FORMULA.md §7.4 sahen explizit vor:

> *Fall C (< 7/10 sign ok): Analytik-Pfad scheitert. Wechsel zu Plan B (deskriptives Paper, vgl. Handoff §4.3).*

**Wir sind in Fall C.** Der Abbruch wurde nach **allen vier Wochen des Analytic-Plans** durchgeführt (nicht nach Woche 2 wie ursprünglich im Handoff vorgeschlagen), weil die Kern-Darstellung (Woche 1) und Variationsrechnung (Woche 2) gut fundiert waren und erst die numerische Validierung in Woche 4 die Unzulänglichkeit der Leading-Order-Formel aufdeckte.

**Wissenschaftlich wertvolles negatives Resultat:**
1. Die Weil-Kern-Architektur selbst ist **korrekt**, liefert aber keine scharfe Leading-Order-Vorhersage.
2. Die Hauptschwäche ist die Approximation $\phi_\chi^+ \approx \phi_\chi^-$: Sektor-Differenz der Grundzustände ist substanziell, nicht sub-leading.
3. Alle einfachen Fenster-Funktionen (Gauss, Rect, Cos²) und sogar der empirische Grundzustand selbst produzieren unzureichende Vorzeichen-Vorhersagen.

---

## 5. Plan B — Deskriptives Paper

Handoff §4.3 sieht Plan B vor: **deskriptives Paper über die 10-Charakter-Studie**.

### 5.1 Was ins Paper kommt

- **Tabelle der 10 Charaktere** mit Gaps, Nullstellen, $L(1)$, $L'(1)/L(1)$, Modul-Struktur.
- **Dokumentation der $\chi_{33}$-Anomalie:** stabil negativer Gap $-0.142$ (N ∈ {200, 400, 600}), Charakter des Moduls 33 = 3 × 11 und seine CRT-Faktorisierung.
- **Dokumentation der vier falsifizierten Vorhersage-Ansätze:**
  1. Einzel-Prädiktor-Regression (Session 4): UBA, CRT, $L'(1)$-Proportionalität (alle $R^2 < 0.32$).
  2. Resonanz-$t_1$-Hypothese (Session 5 Teil 1): $R^2 < 0.11$.
  3. Explicit-Formula-Test-Funktionen (Session 5 Teil 2): alle sechs $R^2 < 0.03$.
  4. **Leading-Order-Weil-Kern (diese Woche 4):** 6/10 Vorzeichen, $R^2 = 0.15$.
- **Analytische Architektur (Woche 1–3) als deskriptiver Rahmen:** Die Weil-Kern-Herleitung bleibt wertvoll als strukturelle Einsicht, auch wenn sie keine prädiktive Leading-Order-Formel liefert. Speziell:
  - Kern-Darstellung $K_\chi^\pm$ (ANALYTIC_KERNEL.md Thm 3.2) ist rigoros.
  - Sektor-Differenz $K^- - K^+$ als Anti-Diagonal-Prim-Kern (ANALYTIC_KERNEL.md Thm 4.1) ist rigoros.
  - Explicit-Formula-Session-5-Scheitern ist analytisch erklärt (ANALYTIC_GAP_FORMULA.md §5.4).
- **Offene Forschungs-Frage:** "Was bestimmt $C_\chi$?" — als **Challenge** an die Community formuliert.

### 5.2 Was NICHT ins Paper kommt (als Behauptung)

- Die Leading-Order-Formel Cor. 4.1 (ANALYTIC_GAP_FORMULA.md) darf **nicht** als prädiktiv verkauft werden. Sie kann als "heuristischer Ansatz" präsentiert werden, der empirisch ungenau ist.
- Die Charakter-spezifischen Vorhersagen (ANALYTIC_GAP_FORMULA.md §4.4) müssen als **falsifiziert** dokumentiert werden.

### 5.3 Titel-Vorschlag

> *"Asymmetric Dirichlet Gaps: A Counter-Example-Driven Analysis of Character-Specific Signatures in the Weil Quadratic Form"*

oder neutraler:

> *"Empirical and Analytic Investigation of Sector Gaps for Dirichlet L-Functions in the v2.1 Weil Framework"*

### 5.4 Zenodo-Qualität

Das Paper ist **qualitativ Zenodo-reif** nach Woche 1–3 Analytik (rigorose Kern-Darstellung, Sektor-Differenz, Siegel-Walfisz-Auflösung) + Woche 4 Falsifikation. Empirische Daten aus 5 Sessions (140 Min Server-Laufzeit für N=600-Scan allein) bieten soliden Testkorpus. Ist geeignet als Arbeits-Preprint (PP-READY-Status), nicht als Journal-Einreichung.

---

## 6. Konsequenzen für das Meta-Paper

### 6.1 META_GALERKIN_BRIDGE.tex §5.5 Cartography

Die Dirichlet-Zeile in der Cartography-Tabelle muss **erneut** aktualisiert werden. Vorherige Versionen:

1. **v0.3 (2026-04-13):** "$C2^{\mathrm{parity}}$ expected closable via v2.1-twist" (optimistisch).
2. **v0.4 (nach Session 4 Teil 3, 2026-04-16):** "$\chi_{12}$ exceptional" (vorübergehend).
3. **v0.5 (nach Session 5 Teil 2, 2026-04-16):** "mixed behavior; 4/10 uniformly positive, 3/10 negative asymptotic, 3/10 near-zero" (deskriptiv).
4. **v0.6 (nach dieser Woche 4-Falsifikation):** "v2.1-twist on leading order does NOT capture Dirichlet gap signs; character-specific structure requires non-leading-order analysis that is currently open."

### 6.2 Drei-Komponenten-Zerlegung

Die Formel $RH_\chi = GBC_\chi + C2_\chi + \mathrm{Hurwitz}_\chi$ bleibt **gültig** — die Woche-4-Falsifikation zeigt nur, dass der $C2_\chi$-Anteil für non-triviale Dirichlet-Charaktere **nicht** durch die Leading-Order-Analyse geschlossen werden kann. $C2_\chi$ bleibt offen für $\chi \neq \chi_0$.

Das ist konsistent mit der Meta-Paper-Aussage (§5.3 Meta-Theorem 1): $C2_X$ ist das "Hilbert-Pólya-Problem" des spezifischen $X$, das für trivial-character-Riemann via Connes-Route geschlossen ist, für nicht-triviale Dirichlet-L **offen bleibt**.

---

## 7. Lessons Learned

### 7.1 Methodologisch

- **Analytische Heuristik ist nicht ausreichend.** Die 4-Wochen-Roadmap (Kern → Grundzustand → Gap-Formel → Validierung) hat alle drei theoretischen Wochen erfolgreich durchgezogen, aber die empirische Validierung war das kritische Filter.
- **Die Heuristik des "charakter-universellen Grundzustands" (Prop. 5.1) wurde zu schnell gemacht.** Ein früher empirischer Check des Grundzustands (z.B. nach Woche 2) hätte die Schwäche der Leading-Order-Formel schon früher aufgezeigt.
- **Rigorosität-pending war ehrlich markiert.** Alle Working Hypotheses wurden als solche gekennzeichnet (ANALYTIC_GROUNDSTATE.md §7.2), sodass die Falsifikation strukturell vorbereitet war.

### 7.2 Mathematisch

- **Die Weil-Kern-Architektur bleibt ein Fortschritt.** Rigoros hergeleitet sind: sektor-projizierte Kerne, Anti-Diagonal-Struktur der Sektor-Differenz, Spektral-Darstellung mit $L'/L$-Polen, Siegel-Walfisz-Entkopplung. Diese Resultate überleben die Falsifikation der Leading-Order-Gap-Formel.
- **Das $C_\chi$-Mysterium bleibt.** Aber es ist jetzt klar, dass es NICHT in der Prim-Summe mit glattem Fenster sitzt; es liegt in der **feinen Struktur der Grundzustand-Differenz** $\phi_\chi^+ - \phi_\chi^-$. Das ist eine subtile Aussage über die **Interaktion** der Primzahl-Phasen mit den **modalen** Strukturen der Galerkin-Matrix, nicht mit dem kontinuierlichen Grundzustand.

### 7.3 Strategisch

- **Plan B ist kein Abstieg, sondern eine Neu-Positionierung.** Ein Paper, das 4 widerlegte Hypothesen und eine rigorose Kern-Theorie (ohne prädiktive Leading-Order) präsentiert, ist **wissenschaftlich substantiell**. Es ist das erste systematische Dokument, das die Grenzen der v2.1-Übertragung auf Dirichlet konkretisiert.
- **Die Meta-Paper-Aussage bleibt robust.** Die Drei-Komponenten-Zerlegung $RH_X = GBC_X + C2_X + \mathrm{Hurwitz}_X$ reformuliert das Dirichlet-Problem: $GBC_\chi$ und $\mathrm{Hurwitz}_\chi$ sind strukturell analog zum Riemann-Fall; $C2_\chi$ bleibt offen. Das ist weniger als ein direkter Beweis, aber es ist ein **klassifizierendes Resultat**.

---

## 8. Aktualisierungs-Log

| Datum | Autor | Änderung |
|---|---|---|
| 2026-04-17 | LG (Opus 4.6, Session 6 Teil 5) | Erste Version: Dokumentation aller 16 getesteten Modelle (Gauss, Rect, Cos², Raw, Empirical), Verdikt Fall C mit Best-Case 6/10 sign_ok und $R^2 = 0.15$ (Gauss $\sigma=0.5$), Empfehlung Plan B. |

---

**Ende von WOCHE4_VALIDATION_RESULT.md.** Die Analytic-Pfad-Exploration ist abgeschlossen. Nächster Schritt: Umstellung auf Plan B (deskriptives Paper).
