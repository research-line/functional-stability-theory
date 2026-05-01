# χ_33-Anomalie: die Struktur bleibt widerspenstig

**Datum:** 2026-04-16 (Session 5 Teil 3)
**Script:** `_scripts/chi33_anomaly.py`
**Status:** **NEGATIVER BEFUND** — keine der getesteten strukturellen Merkmale erklärt, warum χ_33 robust negativ ist, während χ_17, χ_29 (primitive even Charaktere ähnlicher Größe) positiv sind.

---

## 1. Die Ausgangsfrage

Bei N=600 und λ=20000 (konvergierter asymptotischer Bereich):

| Gruppe | Charaktere | gap_asymp |
|---|---|---:|
| **Negativ** | χ_8, χ_13, χ_33 | $-0.049, -0.125, -0.142$ |
| **Positiv** | χ_5, χ_12, χ_17, χ_24, χ_29 | $+0.01$ bis $+0.03$ |
| **Nah-null** | χ_21, χ_60 | $-0.004, +0.005$ |

Warum ist χ_33 robust negativ? Getestet wurden vier strukturelle Hypothesen.

---

## 2. Hypothese 1: Globale Primsum $T_\chi(\lambda = 20000)$

$T_\chi(\lambda) = \sum_{p \leq \lambda} \chi(p) \log(p) / \sqrt{p}$

| Gruppe | $T_\chi(\lambda = 20000)$ |
|---|---:|
| Negativ-Gap (3) | mean $-4.78$, std $0.16$ |
| Positiv-Gap (5) | mean $-4.34$, std $0.33$ |
| Nah-null (2) | mean $-3.47$, std $0.32$ |

**Tendenz erkennbar**, aber nicht diskriminativ:
- Korrelation $T_\chi$ vs gap: $+0.46$ (schwach)
- Gegenbeispiele: χ_5 hat $T = -4.78$ (wie die negativen), aber gap $= +0.016$.
- χ_17 hat $T = -4.25$ (wie die positiven), aber kaum Unterschied zu Negativ-Gruppe.

**Verworfen als Einzelerklärung.**

---

## 3. Hypothese 2: Chebyshev-Bias (#negative vs #positive $\chi(p)$)

Alle getesteten Charaktere haben $\#\{\chi(p) = -1\} - \#\{\chi(p) = +1\} \in \{+5, +24\}$, mittelwertig $\sim +15$. 

**Differenzierung:** Die negative Gruppe hat **geringere Streuung** der Biases (std $0.47$ vs $6.77$ für positive), aber keine klare Schranke.

Korrelation bias vs gap: $-0.31$ (schwach, und falsches Vorzeichen für "mehr Bias → mehr Gap").

**Verworfen als Einzelerklärung.**

---

## 4. Hypothese 3: Ratio gap / $T_\chi$

Wenn gap $\propto T_\chi$ gilt (mit charakter-abhängigem Vorfaktor), sollte die Ratio uniform sein.

| Char | Ratio |
|---|---:|
| **χ_33** | **$+0.031$** (größtes Exzess) |
| **χ_13** | **$+0.025$** |
| **χ_8** | $+0.010$ |
| χ_21 | $+0.001$ |
| χ_60 | $-0.002$ |
| χ_24 | $-0.002$ |
| χ_17 | $-0.003$ |
| χ_5 | $-0.003$ |
| χ_29 | $-0.004$ |
| χ_12 | $-0.008$ |

**Das Vorzeichen der Ratio gap/$T_\chi$ trennt die Gruppen perfekt** (die negativen haben positive Ratio, weil beide Zähler und Nenner negativ sind mit unterschiedlichen Ordnungen). Aber das ist trivial, weil $T_\chi < 0$ für alle Charaktere.

**Interessant:** χ_33 und χ_13 haben **3× größere absolute Ratio** als die positiven Charaktere. Ihre Gap-Werte sind im Verhältnis zu $T_\chi$ überproportional groß.

**Aber:** Das ist eine deskriptive Quantifizierung der Anomalie, keine Erklärung.

---

## 5. Hypothese 4: χ_33-Faktorisierung

Bestätigt: $\chi_{33}(n) = \chi_{-3}(n) \cdot \chi_{-11}(n)$ (Produkt zweier **odd** primitive Charaktere).

| | Parität | $T_\chi$ |
|---|---|---:|
| $\chi_{-3}$ | odd | $-3.49$ |
| $\chi_{-11}$ | odd | $-3.12$ |
| $\chi_{33}$ | even | $-4.57$ |

Die Ratio $T_{\chi_{33}} / (T_{\chi_{-3}} + T_{\chi_{-11}}) \approx 0.69$. Kein klares additives Muster.

**Strukturelle Hypothese geprüft:** χ_21 = χ_{-3} · χ_{-7} hat **dieselbe** Form (even Produkt zweier odd), aber gap $= -0.004$ (nah null, nicht $-0.14$). Die Struktur "odd·odd → even" erklärt die Anomalie **nicht**.

**Verworfen.**

---

## 6. Hypothese 5: Niedrige-Primzahlen-Summe (signed-weighted, erste 50 Primzahlen)

$S_\chi^{(50)} = \sum_{p \in P_{50}} \chi(p) \log(p) / \sqrt{p}$

| Char | $S^{(50)}$ | gap |
|---|---:|---:|
| χ_17 | $-3.96$ | $+0.013$ |
| χ_8 | $-3.58$ | $-0.049$ |
| χ_13 | $-2.98$ | $-0.125$ |
| χ_5 | $-2.84$ | $+0.016$ |
| χ_24 | $-2.96$ | $+0.011$ |
| χ_12 | $-2.32$ | $+0.034$ |
| χ_29 | $-2.03$ | $+0.014$ |
| χ_21 | $-1.99$ | $-0.004$ |
| **χ_33** | **$-1.23$** | $-0.142$ |
| χ_60 | $-0.86$ | $+0.005$ |

**χ_33 hat die zweitkleinste (am wenigsten negative) Niedrig-Primzahlen-Summe**, aber den **stärksten negativen asymptotischen Gap**. Das ist umgekehrt zur Erwartung — wäre das Signal bei kleinen Primzahlen konzentriert, sollte χ_33 positive Gaps haben.

Interpretation: Der Gap für χ_33 wird durch **hohe Primzahlen** gespeist, nicht durch die ersten 50. Das passt zu einer **Resonanz mit höheren L-Funktion-Null-Stellen**, nicht mit der globalen Primsumme.

---

## 7. Was bleibt offen

Keine der vier globalen Metriken unterscheidet die drei negativen Charaktere (χ_8, χ_13, χ_33) **eindeutig** von den fünf positiven. Die schwachen Korrelationen ($R \approx 0.3-0.5$) zeigen: es gibt **Struktur**, aber sie ist **nicht durch einfache Skalare erfassbar**.

**Was könnte die Struktur sein?**

1. **Explicit-Formula-Resonanz:** Der v2.1-Gap hängt über die Weil-Distribution von einer **gewichteten Summe** über alle $L$-Funktion-Null-Stellen ab. Die Gewichte sind Test-Funktionen, die aus der Cos-Basis folgen. Das Vorzeichen-Gesamtergebnis hängt von der **Phase-Ausrichtung** aller Null-Stellen $t_k$ relativ zu diesen Gewichten ab.
   
   **Konkret:** $\mathrm{gap}_\chi \sim -\mathrm{Re} \sum_{\rho} F_\lambda(\rho)$ mit Test-Funktion $F_\lambda$. Für verschiedene Charaktere haben die $\{t_k\}$ unterschiedliche Dichten und Phasen, und das Vorzeichen ist **nicht-trivial**.

2. **$L$-Funktion-Höhen-Spektrum:** Die ersten 5-10 Null-Stellen könnten kollektiv die Phase bestimmen. Unser einziger Prädiktor $t_1$ ist zu grob; eine Analyse mit $\{t_1, t_2, \ldots, t_5\}$ könnte mehr Varianz erklären.

3. **Euler-Produkt-Parität:** Bei χ_33 ist $\chi(p) = -1$ für viele kleine Primzahlen (p=5, 7, 13 alle $-1$ → das Produkt ist stark "oszillatorisch"). Das könnte konstruktive Interferenz negativer Gap-Beiträge bewirken.

---

## 8. Befund

**Die χ_33-Anomalie ist nicht durch globale Einzelmetriken erklärbar.** Sie ist eine **strukturelle Eigenschaft**, die sich vermutlich nur durch:
- direkte Explicit-Formula-Rechnung mit mehreren Null-Stellen, oder
- eine globale Test-Funktions-Analyse der Weil-Distribution

auflösen lässt. Das ist dasselbe Problem wie in der Multivariate-Regression: Einzel- und Mehrfach-Parameter-Modelle mit $\{t_1, L(1), -L'/L, \log D\}$ versagen, weil der Gap eine **Funktion unendlich vieler Null-Stellen** ist.

---

## 9. Konsequenzen

### 9.1 Für das Meta-Paper

Die aktuelle deskriptive Formulierung ($C_\chi$ ist charakter-spezifisch, keine Einzelprädiktor-Erklärung) ist **bestätigt**. Keine Revision nötig.

### 9.2 Für das Dirichlet-Paper

Wenn geschrieben, wird das Paper **deskriptiv** sein müssen: Tabelle der gemessenen Gaps, Analyse der Ratio gap/$T_\chi$, Hinweis auf Resonanz-Interpretation der Anomalien (χ_8, χ_13, χ_33), aber **keine prädiktive Formel**.

### 9.3 Nächster theoretischer Schritt

**Task #36 (Analytische C_χ-Ableitung)** bleibt die richtige strategische Priorität. Die Weil-Distribution mit Cos-Basis-Test-Funktion ist der einzige Weg, der **potentiell** das Vorzeichen vorhersagen könnte — falls das in geschlossener Form überhaupt möglich ist.

**Alternative:** Ein Explicit-Formula-Numerik-Experiment, das die Summe $\sum_\rho F_\lambda(\rho)$ **direkt** für χ_33 vs χ_17 berechnet (mit gemeldeten Null-Stellen aus LMFDB). Das wäre ein numerischer Test der Resonanz-Interpretation, ohne vollständige Analytik.

---

## 10. Artefakte

- `_scripts/chi33_anomaly.py` — Analyse-Script
- Dieses Dokument

---

**Session 5 Teil 3 Fazit:** Die χ_33-Anomalie ist **bestätigt als echt** (nicht Truncation-Artefakt) und **bestätigt als nicht durch globale Metriken erklärbar**. Die Antwort liegt in der Feinstruktur der Weil-Distribution mit allen Null-Stellen gleichzeitig — eine Aufgabe für Task #36 (analytische Ableitung) oder ein gezieltes Explicit-Formula-Experiment.
