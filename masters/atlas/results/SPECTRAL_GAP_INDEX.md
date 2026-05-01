# Spectral-Gap-Index Probe

**Datum:** 2026-04-16
**Skript:** `_scripts/spectral_gap_index.py`
**Motivation:** Meta-Framework 'verborgene Schoenheit' #2 - Verhaeltnis gamma^(1)/sqrt(L) als charakter-spezifischer Gap-Index.

## Parameter
- lambda = 20000, L = 9.9035, sqrt(L) = 3.1470

## Daten pro Charakter

| chi | D | gamma^(1) | gap | C_chi | I_ratio=g/sqrtL | R_chi=1/g^3+1/g2^3 |
|---|---:|---:|---:|---:|---:|---:|
| chi_5 | 5 | 6.648 | +0.01603 | +0.00509 | 2.113 | 0.00446 |
| chi_8 | 8 | 4.900 | -0.04902 | -0.01558 | 1.557 | 0.01075 |
| chi_12 | 12 | 3.805 | +0.01154 | +0.00367 | 1.209 | 0.02149 |
| chi_13 | 13 | 3.119 | -0.12504 | -0.03973 | 0.991 | 0.03559 |
| chi_17 | 17 | 3.728 | +0.00886 | +0.00282 | 1.185 | 0.02489 |
| chi_21 | 21 | 2.315 | -0.00423 | -0.00134 | 0.736 | 0.08576 |
| chi_24 | 24 | 2.689 | +0.00866 | +0.00275 | 0.854 | 0.05820 |
| chi_29 | 29 | 1.794 | +0.01439 | +0.00457 | 0.570 | 0.17991 |
| chi_33 | 33 | 2.997 | -0.14221 | -0.04519 | 0.952 | 0.04843 |
| chi_60 | 60 | 1.881 | +0.00521 | +0.00166 | 0.598 | 0.16615 |

## Signed correlations (C_chi)

| Feature | Pearson R | R^2 |
|---|---:|---:|
| gamma_chi^(1) | +0.0532 | 0.0028 |
| 1/gamma^(1) | +0.1429 | 0.0204 |
| 1/gamma^(1)^2 | +0.2195 | 0.0482 |
| 1/gamma^(1)^3 (Session 8 R_chi single-zero) | +0.2702 | 0.0730 |
| R_chi = sum 1/gamma^(k)^3 (k<=2) | +0.2524 | 0.0637 |
| gamma^(1)/sqrt(L) | +0.0532 | 0.0028 |
| L(1,chi) | -0.3489 | 0.1217 |
| log D | -0.0520 | 0.0027 |

## Magnitude correlations (|C_chi|)

| Feature | Pearson R | R^2 |
|---|---:|---:|
| gamma_chi^(1) | +0.0054 | 0.0000 |
| 1/gamma^(1) | -0.1775 | 0.0315 |
| 1/gamma^(1)^2 | -0.2418 | 0.0585 |
| 1/gamma^(1)^3 (Session 8 R_chi single-zero) | -0.2817 | 0.0793 |
| R_chi = sum 1/gamma^(k)^3 (k<=2) | -0.2659 | 0.0707 |
| gamma^(1)/sqrt(L) | +0.0054 | 0.0000 |
| L(1,chi) | +0.3184 | 0.1014 |
| log D | +0.0063 | 0.0000 |

## Top-4 Predictors of |C_chi|

| Rank | Feature | R | R^2 |
|:---:|---|---:|---:|
| 1 | L(1,chi) | +0.3184 | 0.1014 |
| 2 | 1/gamma^(1)^3 (Session 8 R_chi single-zero) | -0.2817 | 0.0793 |
| 3 | R_chi = sum 1/gamma^(k)^3 (k<=2) | -0.2659 | 0.0707 |
| 4 | 1/gamma^(1)^2 | -0.2418 | 0.0585 |

## Fazit

**Hypothese widerlegt.** Kein einfacher Index auf $\gamma_\chi^{(1)}$ bzw. dessen Kombinationen erreicht substantielle Korrelation mit $C_\chi$. Alle $R^2 < 0.13$; Vorzeichen-Match bei Kombinationen max 7/10 (chance level).

**Entscheidendes Gegenbeispiel:**
- $\chi_{17}$: $\gamma^{(1)} = 3.73$, $C_\chi = +0.003$ (winzig positiv).
- $\chi_{33}$: $\gamma^{(1)} = 3.00$, $C_\chi = -0.045$ (stark negativ, Faktor 15).
- Ähnliche Spektrallücke, entgegengesetzte Vorzeichen, 15× unterschiedliche Magnitude.

**Strukturelle Einsicht:** Die niedrigste Nullstelle $\gamma^{(1)}$ allein kodiert die char-spezifische Gap-Signatur nicht. $L(1,\chi)$ ist der stärkste Einzel-Prädiktor, aber mit nur $R^2 \approx 0.10$ weit von prädiktiv entfernt. Konsistent mit Session 4/5 (multivariate Regression $R^2 = 0.42$, Session 7 T1 Sektor-Differenz $R^2 = 0.49$).

**Positiver Inhalt (Einschränkung des Suchraums):** Die char-Signatur hängt nicht-trivial von **mehreren** Nullstellen UND von Prim-Phasen UND von $L(1,\chi)$ ab. Eine Ein-Parameter-Reduktion existiert nicht. Das verstärkt die These, dass $C^{(2)}_\chi$ strukturell komplex (Connes-Niveau) ist, nicht durch niedrige Zero-Daten kollabierbar.

**Fünfter falsifizierter Predictor.** Zusammen mit den vier in Sessions 4–7 bereits dokumentierten (multivariate Regression, niedrig-Zero-Resonanz, Explicit-Formula-Testfunktionen, Leading-Order-Weil-Kern, Galerkin-Autokorrelation) ist dies die fünfte qualitativ verschiedene Ansatzklasse, die an $C^{(2)}_\chi$ scheitert. Der Suchraum für "einfache" Prädiktoren ist systematisch ausgeleuchtet und leer.
