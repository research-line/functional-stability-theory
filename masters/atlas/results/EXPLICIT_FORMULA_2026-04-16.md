# Explicit-Formula-Experiment: naive Null-Stellen-Resonanz ist falsifiziert

**Datum:** 2026-04-16 (Session 5 Teil 4)
**Script:** `_scripts/zeros_and_testfunctions.py` (9 Min lokal, 10 Charaktere)
**Status:** **NEGATIVER BEFUND** — keine der sechs getesteten Explicit-Formula-Test-Funktionen sagt das Vorzeichen oder die Magnitude des gemessenen Gaps vorher.

---

## 1. Was getestet wurde

Für alle 10 Charaktere wurden die ersten 6–13 nicht-trivialen Null-Stellen $\{t_1, \ldots, t_n\}$ via mpmath berechnet (Minima von $|L(1/2+it, \chi)|$ bis $t_{\max} = 20$).

Für jede der folgenden Test-Funktionen $F_\lambda(t)$ wurde die Summe $S_\chi = 2 \sum_k F_\lambda(t_k)$ (Paar $\pm t_k$) berechnet und gegen den gemessenen asymptotischen Gap bei $\lambda = 20000$ korreliert:

| Name | $F_\lambda(t)$ | Motivation |
|---|---|---|
| F_weil | $\mathrm{Re}(-1/\rho) = (1/2)/(1/4+t^2)$ | Explicit-Formula-Hauptterm |
| F_cos_weil | $\cos(Lt) \cdot F_\text{weil}$ | Weil-Kern mit Frontier-Oszillation |
| F_gauss | $\exp(-t^2/(2/L^2))$ | Gauss, Breite $1/L$ |
| F_fejer | $\mathrm{sinc}^2(Lt/2)$ | Fejer-Kern aus Cos-Basis-Projektion |
| F_inv_t2 | $1/t^2$ | Dominanz niedrigster Null-Stellen |
| F_sinc | $\sin(Lt)/(Lt)$ | sinc |

---

## 2. Ergebnisse

| Test-Funktion | $r$ | $R^2$ | sign_ok |
|---|---:|---:|---:|
| F_weil | $+0.10$ | $0.009$ | 4/10 |
| F_cos_weil | $-0.07$ | $0.005$ | 4/10 |
| F_gauss | $-0.10$ | $0.009$ | 4/10 |
| F_fejer | $+0.17$ | $0.029$ | 4/10 |
| F_inv_t2 | $+0.11$ | $0.011$ | 4/10 |
| F_sinc | $+0.06$ | $0.004$ | 4/10 |

**Multivariate Regression auf allen 6 Test-Funktionen:**
- $R^2 = 0.31$, sign_ok = 6/10
- Das ist schlechter als die frühere Multivariate-Regression auf $\{t_1, L(1), -L'/L, \log D\}$ ($R^2 = 0.42$).

**Konklusion:** Keine der naiven Null-Stellen-Summen-Modelle sagt das Gap-Verhalten vorher. **Die Resonanz-Interpretation (gap als gewichtete Summe über Null-Stellen) in ihrer einfachen Form ist falsifiziert.**

---

## 3. Was das bedeutet

### 3.1 Die einfache Resonanz-Hypothese ist tot

Die Session-4-Vermutung: *"Der Gap ist dominiert von den niedrigsten Null-Stellen"* ist empirisch falsch. Weder $1/t_1^2$-Gewichtung noch Fejer-Kern noch Gauss-Summen produzieren eine verwertbare Vorhersage. Alle univariaten $R^2 < 0.03$.

### 3.2 Vorzeichen-Vorhersage schlägt Zufall nicht

4/10 sign_ok ist **schlechter** als der Base-Rate-Anteil (6/10 positive Gaps, also 6/10 durch Konstante vorhersagbar). Das heißt: die Test-Funktionen tragen **aktiv negative Information** zum Vorzeichen bei.

### 3.3 Die echte Struktur ist komplexer

Drei mögliche Gründe:

1. **Die Test-Funktion ist nicht elementar.** Die wahre $F_\lambda$ aus Cos-Basis + Weil-Distribution hat **kompliziertere Struktur**, z.B. Interferenzen zwischen den Moden $n=1, 2, 3, \ldots$. Ein Explicit-Formula-Ansatz, der nur eine Test-Funktion nutzt, reicht nicht.

2. **Der Gap hat Primzahl-Komponente.** Der v2.1-Gap könnte eine Mischung aus Null-Stellen-Summe und Primzahl-Beitrag sein, die sich für verschiedene Charaktere unterschiedlich ausbalancieren. Das passt zur Chebyshev-Bias-Analyse (CHI33_ANOMALIE): die Primzahlen-Struktur ist teilweise aussagekräftig.

3. **Archimedischer Effekt.** Die $\Gamma(s/2)$-Faktor-Analyse bei $t = t_k$ ist nicht-trivial und könnte je nach Charakter unterschiedliche Beiträge liefern. Das ist bisher in den Test-Funktionen nicht berücksichtigt.

---

## 4. Die 10 Charaktere mit ihren Null-Stellen

| Charakter | $t_1$ | $t_2$ | $t_3$ | $t_4$ | $t_5$ | $t_6$ |
|---|---:|---:|---:|---:|---:|---:|
| $\chi_5$ | 6.648 | 9.831 | 11.959 | 16.034 | 17.567 | 19.541 |
| $\chi_8$ | 4.900 | 7.628 | 10.807 | 12.311 | 15.196 | 17.022 |
| $\chi_{12}$ | 3.805 | 6.692 | 8.891 | 11.188 | 12.966 | 15.181 |
| $\chi_{13}$ | 3.119 | 7.232 | 8.625 | 10.336 | 12.617 | 15.148 |
| $\chi_{17}$ | 3.728 | 5.636 | 7.283 | 10.617 | 11.978 | 13.210 |
| $\chi_{21}$ | 2.315 | 5.780 | 7.655 | 9.465 | 10.611 | 12.629 |
| $\chi_{24}$ | 2.689 | 5.292 | 6.972 | 9.225 | 10.446 | 12.631 |
| $\chi_{29}$ | 1.794 | 5.317 | 6.758 | 8.631 | 10.420 | 11.195 |
| $\chi_{33}$ | 2.997 | 4.459 | 6.362 | 7.803 | 10.169 | 11.871 |
| $\chi_{60}$ | 1.881 | 3.986 | 5.461 | 7.166 | 8.485 | 9.943 |

**Wichtige Beobachtung:** χ_33 und χ_17 haben ähnliche $t_1$ (3.00 vs 3.73) und ähnliche höhere Null-Stellen — aber entgegengesetzte Gap-Vorzeichen ($-0.142$ vs $+0.013$). Die Null-Stellen-Verteilung allein kann den Unterschied **nicht** erklären.

---

## 5. Revidierte Interpretation

Nach Session 4 (Resonanz-Hypothese) und Session 5 (multivariate Regression, χ_33-Anomalie, Explicit-Formula-Experiment) ist klar:

> **Der v2.1-Dirichlet-Gap hängt von einer Struktur ab, die nicht durch globale Charakter-Invarianten (Null-Stellen, $L(1)$, Regulator, Klassenzahl, Chebyshev-Bias, Primsum) erklärbar ist — zumindest nicht in den bisher getesteten linearen und multiplikativen Formen. Die Variabilität ist entweder (i) durch feinere, kollektiv-interferente Strukturen bestimmt, oder (ii) durch ein bislang nicht identifiziertes Funktional, das mehrere dieser Größen in nichttrivialer Weise kombiniert.**

---

## 6. Konsequenz für das Programm

### 6.1 Task #36 (analytische C_χ-Ableitung)

Die rein theoretische Ableitung aus Weil-Distribution + Cos-Basis bleibt **der einzige nicht-ausgeschlossene Weg**. Sie müsste folgende Struktur produzieren:

$$
\mathrm{gap}_\chi(\lambda) = \mathrm{Re} \int_{-\infty}^{\infty} G_\lambda(t) \cdot \frac{L'/L(1/2+it, \chi)}{1} dt + (\text{boundary terms})
$$

mit einer **nicht-trivialen Test-Funktion** $G_\lambda$, die aus der Matrixdarstellung von Cos/Sin-Basis auf $[-L, L]$ folgt. Diese Ableitung ist nicht-trivial und erfordert mehrere Wochen.

### 6.2 Dirichlet-Paper

Das Dirichlet-Paper wird **deskriptiv** sein müssen. Empirische Tabelle der gemessenen Gaps, Kommentar zu χ_33 als Gegenbeispiel, keine Vorhersage-Formel. Der Abschnitt "Was bestimmt $C_\chi$?" wird als **offenes Problem** geschrieben.

### 6.3 Meta-Paper

Bereits ausreichend ehrlich formuliert (Version nach REVISION_N400 + N600). Keine weitere Revision nötig.

---

## 7. Status

- ✅ Null-Stellen für alle 10 Charaktere berechnet (mpmath, 9 Min)
- ✅ 6 Test-Funktionen evaluiert — alle verwerfen sich
- ✅ Multivariate Regression über alle Test-Funktionen: $R^2 = 0.31$, schlechter als Primsum-Features
- ❌ **Resonanz-Interpretation falsifiziert** in einfacher Form
- ⏳ Task #36 (analytische Ableitung) bleibt als einzige Option offen

---

## 8. Artefakte

- `_scripts/zeros_and_testfunctions.py`
- `_results/zeros_all_chars.json` — Null-Stellen aller 10 Charaktere
- `_results/explicit_formula_experiment.json` — Summen + Fits
- Dieses Dokument

---

**Session 5 Teil 4 Fazit:** Das Explicit-Formula-Experiment ist durchgeführt und die naive Resonanz-Interpretation **empirisch falsifiziert**. Das war ein sauberer Test einer plausiblen Hypothese, und das negative Ergebnis hat Informationsgehalt: die Gap-Struktur ist **nicht** durch einfache Null-Stellen-Gewichtung erfassbar. Der Weg vorwärts ist (a) echte analytische Ableitung aus der Weil-Distribution mit korrekter Test-Funktion (mehrere Wochen, Task #36), oder (b) Akzeptanz der deskriptiven Dirichlet-Paper-Formulierung ohne prädiktive Theorie.
