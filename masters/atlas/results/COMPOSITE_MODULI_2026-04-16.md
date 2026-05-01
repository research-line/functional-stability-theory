# Composite-Moduli-Test: CRT-Hypothese falsifiziert — χ₁₂ ist einzigartig

**Datum:** 2026-04-16 (Session 4 fortgesetzt)
**Script:** `_scripts/composite_moduli_server.py`
**Server:** ellmos-services (18.2 Min, 48 Datenpunkte)
**Status:** **HYPOTHESE FALSIFIZIERT** — χ₁₂-Stabilität ist nicht durch zusammengesetzte teilerfremde Modul-Struktur erklärbar.

---

## 1. Versuchsaufbau

Test der Hypothese aus `THEORIE_CRT_FAKTORISIERUNG.md`: "Primitive reelle Charaktere mit zusammengesetztem teilerfremden Modul zeigen Gap-Stabilität".

6 Charaktere getestet bei $N=200$, $\lambda \in [100, 20000]$:

| D | Faktorisierung | Typ |
|---:|---|---|
| 17 | prim | Kontrolle (erwartet: weak) |
| 21 | $3 \cdot 7$, teilerfremd | Test (erwartet: strong) |
| 24 | $8 \cdot 3$, teilerfremd | Test (erwartet: strong) |
| 29 | prim | Kontrolle (erwartet: weak) |
| 33 | $3 \cdot 11$, teilerfremd | Test (erwartet: strong) |
| 60 | $4 \cdot 3 \cdot 5$, drei Faktoren | Test (erwartet: sehr strong) |

---

## 2. Konsolidierte Ergebnisse

| $D$ | Typ | pos/n | mean | std | slope | slope_err |
|---:|---|---:|---:|---:|---:|---:|
| 17 | prim | 7/8 | $+0.026$ | $0.040$ | $+0.002$ | $0.009$ |
| 21 | komp. tf | 5/8 | $+0.094$ | $0.179$ | $-0.010$ | $0.042$ |
| 24 | komp. tf | 5/8 | $+0.001$ | $0.030$ | $+0.005$ | $0.007$ |
| 29 | prim | 7/8 | $+0.132$ | $0.235$ | $-0.046$ | $0.051$ |
| 33 | komp. tf | **4/8** | $-0.011$ | $0.051$ | $-0.020$ | $0.009$ |
| 60 | komp. tf | 6/8 | $+0.157$ | $0.188$ | $-0.029$ | $0.042$ |

**Vergleich mit Referenz-Charakteren aus vorherigem Lauf:**

| $D$ | Typ | pos/n | mean | slope |
|---:|---|---:|---:|---:|
| 5 | prim | 5/8 | $+0.160$ | $-0.166$ |
| 8 | $2^3$ | 5/8 | $+0.055$ | $-0.027$ |
| **12** | **komp. tf** | **8/8** | $+0.131$ | $+0.002$ |
| 13 | prim | 6/8 | $+0.271$ | $-0.218$ |

**Nur χ₁₂ erfüllt weiterhin 8/8 Positiv + Slope ≈ 0.**

---

## 3. Hypothesen-Auswertung

### 3.1 CRT-Hypothese: FALSIFIZIERT

Die Vorhersage "komposite teilerfremde Moduln → Stabilität" wird nicht bestätigt:
- **D=21** (3·7): 5/8 positiv, aber mit großer Streuung (std 0.18) — **oszillatorisch**.
- **D=24** (8·3): 5/8 positiv, aber Mean $\approx 0$ — praktisch keine Dominance.
- **D=33** (3·11): **4/8 positiv** — weniger als Hälfte, Mean negativ.
- **D=60** (3·4·5): 6/8 positiv, aber slope fallend.

Kein zusammengesetzter Charakter reproduziert die χ₁₂-Stabilität.

### 3.2 Primzahl-Kontrolle: ÜBERRASCHEND STABIL

- **D=17**: 7/8 positiv, aber alle Werte sehr klein ($|\mathrm{gap}| < 0.13$). Slope $+0.002 \pm 0.009$ — konsistent mit 0!
- **D=29**: 7/8 positiv, aber starke Fluktuation (std 0.24).

**D=17 zeigt ein *ähnliches* Verhalten wie χ₁₂ (Slope ≈ 0)** — aber mit **viel kleinerem Mean**. Das ist interessant: es gibt mehr als nur einen "stabilen" Charakter, aber die Stärke variiert.

### 3.3 Neue Klassifikation (differenziert)

| Klasse | Charakteristik | Beispiele |
|---|---|---|
| **Strong stabil** | Alle $> 0$, slope $\approx 0$, Mean $\gtrsim 0.1$ | $\chi_{12}$ einzigartig |
| **Weak stabil** | Fast alle $> 0$, slope $\approx 0$, aber Mean klein | $\chi_{17}$ |
| **Oszillatorisch** | Gemischte Vorzeichen, größerer slope_err | $\chi_5, \chi_8, \chi_{13}, \chi_{21}, \chi_{29}, \chi_{33}, \chi_{60}$ |

Die **Mehrheit** der getesteten Charaktere ist oszillatorisch. **Stabilität (slope ≈ 0)** tritt in zwei Varianten auf (stark und schwach). χ₁₂ ist unter den getesteten der einzige "strong stabil".

---

## 4. Warum ist χ₁₂ so besonders?

### 4.1 Strukturelle Eigenheiten

$\chi_{12}$ = $\chi_{-3} \cdot \chi_{-4}$ — Produkt zweier odd characters mit **minimalen** Leitern 3 und 4. Das ist die **kleinste** fundamentale Diskriminante, die als Produkt zweier odd characters mit teilerfremden Leitern geschrieben werden kann. Andere Produkte diesen Typs:

| Form | Leiter | D | stabil? |
|---|---|---:|---|
| $\chi_{-3} \cdot \chi_{-4}$ | (3, 4) | 12 | **JA** |
| $\chi_{-3} \cdot \chi_{-7}$ | (3, 7) | 21 | NEIN |
| $\chi_{-3} \cdot \chi_{-11}$ | (3, 11) | 33 | NEIN |
| $\chi_{-4} \cdot \chi_{-7}$ | (4, 7) | 28 | (nicht getestet) |
| $\chi_{-7} \cdot \chi_{-8}$ | (7, 8) | 56 | (nicht getestet) |

**Mögliche Erklärung:** Die **absolute Größe** der Leiter beeinflusst die Kovarianz-Struktur der zugrunde liegenden L-Funktion. Kleinste Leiter (3, 4) → niedrigste Konditionalität → stabilstes Gap.

Das wäre **quantitativ** eine Form der Analytic-Conductor-Theory: kleinere analytischer Leiter → weniger niedrigliegende Null-Stellen → weniger Oszillation.

### 4.2 Alternative: Numerisches Artefakt?

Eine wichtige Möglichkeit: der "konstante Gap bei χ₁₂" könnte ein **N-dependentes Artefakt** sein, das bei $N \to \infty$ verschwindet. Tests mit $N = 400, 600$ würden das klären.

### 4.3 Metaphysische Interpretation

Die Beobachtung ist: **Asymptotische Gap-Konstanz ist selten und von spezifischer Modul-Struktur abhängig.** Das ist ein wertvoller negativer Befund: die Meta-Paper-Vorhersage "Even-Dominance für Dirichlet even" muss als **generischer Trend** formuliert werden, nicht als **quantitative** Aussage über Gap-Asymptotik.

---

## 5. Revidierte Meta-Paper-Formulierung (Vorschlag)

Die §5.5 Dirichlet-Zeile sollte die quantitative Aussage abschwächen:

**Alt** (Session 4 Teil 3):
> χ₁₂ perfect 9/9 up to λ=50,000 with asymptotically constant gap ≈ 0.13; fitted slope −0.002±0.012

**Revidiert** (nach Composite-Test):
> Numerical evidence varies widely across characters. Most even characters show oscillatory gaps with amplitude $\lesssim 0.15$. χ₁₂ is exceptional: 9/9 positive, slope consistent with 0, but this stability is not reproduced by other composite-coprime-modulus characters ($D = 21, 24, 33, 60$). The asymptotic structure of $C2^{\mathrm{parity}}$ for Dirichlet even characters remains a delicate open question.

Das ist ehrlicher und offener für weitere Forschung.

---

## 6. Implikationen für Dirichlet-Paper

- χ₁₂ bleibt ein **starkes Einzelbeispiel**.
- Die **Theorie der Einzigartigkeit** muss entwickelt werden (analytic conductor + Klassenzahl-Analyse).
- Das **geplante Dirichlet-Paper** könnte als zentrale Behauptung haben: "Die v2.1-Methode zeigt Even-Dominance für alle primitiven reellen Dirichlet-Charaktere, aber mit charakter-spezifischen Raten; χ₁₂ ist ein außergewöhnliches Beispiel stabiler Asymptotik."

---

## 7. Nächste Schritte

### 7.1 Unmittelbar (Session 5)

1. **χ₁₂ bei höherem N prüfen** ($N = 400$ oder $600$): ist die Konstanz wirklich echt oder Truncation-Artefakt?
2. **D=28, D=56** testen: weitere Produkte zweier odd characters, um Leiter-Hypothese zu festigen.
3. Meta-Paper §5.5 revidieren mit ehrlicher Formulierung.

### 7.2 Mittelfristig

1. **Analytic-Conductor-Analyse** für Klassifikation (Iwaniec-Kowalski Ch 5).
2. **Null-Stellen-Verteilung** von $L(s, \chi_{12})$ vs. $L(s, \chi_5)$ bei $s=1/2$.
3. **Theorie der Gap-Fluktuation**: was genau bestimmt Amplitude und Vorzeichen?

---

**Session 4 Teil 6 Fazit:** Die CRT-Hypothese ist falsifiziert. χ₁₂ ist **einzigartig** unter den getesteten Charakteren — sowohl in der 4-Char-Gruppe als auch im 6-D-Composite-Test. Die Ursache bleibt unklar und ist ein zentrales Thema für Session 5. Die Meta-Paper-Vorhersage muss abgeschwächt werden, um die echte Variabilität zu reflektieren.
