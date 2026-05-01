# Multi-variate Regression — die Resonanz-Hypothese hält nicht

**Datum:** 2026-04-16 (Session 4 real-real-Ende)
**Script:** `_scripts/compute_t1_all.py`
**Status:** **REVISION** — die 3-Punkt-Korrelation aus LMFDB-Recherche war Zufall.

---

## 1. Alle 10 getesteten Charaktere

Berechne $t_1$ numerisch via mpmath (L(1/2+it, χ)-Minima-Suche) für alle in Session 4 getesteten primitiven reellen even Charaktere:

| Charakter | $D$ | $t_1$ | $L(1,\chi)$ | $-L'/L$ | mean gap (N=200) | pos/8 |
|---|---:|---:|---:|---:|---:|---:|
| $\chi_5$ | 5 | 6.648 | 0.430 | 1.833 | $+0.160$ | 5/8 |
| $\chi_8$ | 8 | 4.900 | 0.623 | 2.256 | $+0.055$ | 5/8 |
| $\chi_{12}$ | 12 | **3.805** | 0.760 | 2.629 | $+0.131$ | **8/8** |
| $\chi_{13}$ | 13 | 3.119 | 0.663 | 2.704 | $+0.271$ | 6/8 |
| $\chi_{17}$ | 17 | 3.728 | 1.016 | 2.954 | $+0.026$ | 7/8 |
| $\chi_{21}$ | 21 | 2.315 | 0.684 | 3.154 | $+0.094$ | 5/8 |
| $\chi_{24}$ | 24 | 2.689 | 0.936 | 3.280 | $+0.001$ | 5/8 |
| $\chi_{29}$ | 29 | 1.794 | 0.612 | 3.460 | $+0.132$ | 7/8 |
| $\chi_{33}$ | 33 | 2.997 | 1.333 | 3.584 | $-0.011$ | 4/8 |
| $\chi_{60}$ | 60 | 1.881 | 1.066 | 4.159 | $+0.157$ | 6/8 |

---

## 2. Univariate Regressionen

| Prädiktor | Slope | Intercept | $R^2$ |
|---|---:|---:|---:|
| $t_1$ | $+0.003$ | $+0.093$ | **0.002** (keine Korrelation) |
| $L(1,\chi)$ | $-0.180$ | $+0.248$ | **0.313** (beste Einzel-Korrelation) |
| $-L'/L$ | $-0.026$ | $+0.180$ | $0.042$ |
| $\log D$ | $-0.025$ | $+0.173$ | $0.043$ |

**Keine einzelne Variable erklärt mehr als 31% der Varianz.**

## 3. Multivariate Regression

Mit allen vier Prädiktoren als linearer Kombination:
$$
\mathrm{gap} = -0.35 \cdot t_1 + 0.52 \cdot L(1) + 21.9 \cdot (-L'/L) - 21.3 \cdot \log D - 3.57
$$
$R^2 = 0.70$ — aber bei 10 Datenpunkten und 4 Features ist das **Overfitting-verdächtig**.

## 4. Exponentieller Fit: `gap ~ exp(α · t_1)`

Der vorher gefundene Fit `gap ≈ 0.130 · exp(-0.9 · (t_1 − 3.805))` aus 3 Datenpunkten:
- Slope = $-0.058$ (statt $-0.9$!)
- $R^2 = 0.017$

**Die exponentielle Abhängigkeit ist mit 10 Datenpunkten nicht mehr vorhanden.**

---

## 5. Interpretation

### 5.1 Was falsch war

Der "Durchbruch" mit drei Punkten (χ_5, χ_8, χ_12) war ein **statistisches Artefakt**. Mit 10 Punkten zeigt sich: die tatsächliche Varianz der Gap-Werte wird nicht durch niedrige Null-Stellen erklärt.

**Gegenbeispiele:**
- χ_29 hat sehr niedriges $t_1 = 1.79$, aber gap=+0.13 (nicht besonders groß).
- χ_24 hat $t_1 = 2.69$, aber gap=+0.001 (fast null trotz niedrigem $t_1$).
- χ_13 hat $t_1 = 3.12$ und gap=+0.27 (größter gap in den Daten!) — das wäre laut Hypothese nicht zu erwarten.

### 5.2 Was die Daten wirklich zeigen

**Die Gap-Werte sind komplex verteilt**, nicht monoton in einem einzelnen Parameter. Die multivariate Regression verbessert sich auf $R^2 = 0.70$, aber das ist mit 4 Features und 10 Punkten wahrscheinlich Overfitting.

**L(1, χ) ist der beste Einzel-Prädiktor**, mit negativem Slope (höheres L(1) → kleinerer gap), aber $R^2 = 0.31$ ist schwach.

### 5.3 Mögliche Gründe für mangelnde Korrelation

1. **N=200-Daten sind Truncation-verzerrt.** Wir wissen aus Session 4 Teil 4, dass bei N=400 die Werte ganz anders aussehen (χ_5 fällt von 0.16 auf 0.01). Mit N=400-Daten für alle 10 Charaktere könnte die Struktur anders aussehen.

2. **Komplizierteres funktionales Muster.** Die Abhängigkeit könnte von mehreren Null-Stellen $t_1, t_2, \ldots$ gleichzeitig sein, nicht nur $t_1$.

3. **Chebyshev-Bias oder Klassen-Zahl-Effekte** könnten eine Rolle spielen, die nicht in den getesteten Prädiktoren erscheint.

---

## 6. Nächste Schritte

### 6.1 Unmittelbar

**Multi-N-Erweiterung.** Server-Lauf mit N=400 für alle 10 Charaktere. Dann Regression wiederholen. Das könnte klären, ob die bisherigen Daten Truncation-verzerrt waren.

### 6.2 Mittelfristig

Wenn N=400-Daten auch keine klare Korrelation mit t_1 zeigen:
- **L(1/2, χ) testen** als Prädiktor — der zentrale Wert könnte relevant sein.
- **Höhere Null-Stellen** $t_2, t_3$ einbeziehen.
- **Explicit-Formula-Integral** numerisch berechnen und direkt mit gap vergleichen.

---

## 7. Revidierte Meta-Paper-Formulierung

Die vorherige Formulierung "exponentieller Fit mit $t_1$" wird zurückgenommen. Die ehrliche Aussage:

> **For primitive real even Dirichlet characters, the v2.1 Galerkin gap varies over roughly one order of magnitude (from $\approx 0$ to $\approx 0.27$) across the tested 10 characters. The $N = 200$ data show no single-parameter correlation (all $R^2 < 0.32$). Understanding the character-specific amplitude $C_\chi$ requires either higher-$N$ data or more sophisticated structural analysis.**

---

## 8. Wissenschaftliche Einsicht

**Die Physik-Forschungs-Lektion:** 3-Punkt-"Korrelationen" sind statistisch nicht signifikant. Mit 10 Punkten wurde klar, dass die Struktur **reicher und widerspenstiger** ist als die Drei-Punkt-Extrapolation vermuten ließ.

χ_12 bleibt ein interessanter Einzelfall (8/8 perfekt positiv), aber die Annahme, dass seine Stabilität durch $t_1 = 3.805$ allein erklärt wird, ist falsch. **χ_12 ist wirklich besonders**, und die Ursache ist weiterhin offen.

---

**Session 4 FINALES Fazit:** Die scheinbare Lösung via Null-Stellen-Korrelation war ein statistisches Artefakt von zu wenigen Datenpunkten. Die ehrliche Situation ist: mit N=200-Daten kein einfacher Prädiktor gefunden. Die nächste Priorität für Session 5 ist der **N=400-Lauf für alle 10 Charaktere**, um die Truncation-Unsicherheit zu eliminieren und dann die Regression zu wiederholen.
