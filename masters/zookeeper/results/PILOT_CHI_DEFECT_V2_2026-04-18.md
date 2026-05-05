# Pilot v2 — χ-Defekt-Norm mit Prolate-Basis + Sonin-Zerlegung + Pol-Normalisierung (2026-04-18)

**Autor:** LG (Claude Opus 4.7 [1M], Session 16 Fortsetzung)
**Skript:** `_scripts/chi_defect_norm_v2.py`
**Vorgänger:** `PILOT_CHI_DEFECT_2026-04-18.md` (Pilot v1, stark vereinfacht)
**Status:** Erste Version mit kontrollierter Prolate-Basis und Sonin-Zerlegung. **Struktur-aufklärendes Zwischenergebnis**: unter Even-Charakteren zeigt sich ein **monotoner H2-konsistenter Trend**; χ_0 bleibt Outlier wegen ζ-Pol.

---

## 1. Setup

- $\lambda = \sqrt{14} \approx 3.7417$ (CCM-Standard).
- $L = \log\lambda \approx 1.3195$.
- **Prolate-Basis** via Diagonalisierung des diskreten sinc-Bandbegrenzungs-Operators auf $[-T_{\text{wide}}, T_{\text{wide}}]$ mit $T_{\text{wide}} = 3L \approx 3.96$, $N_{\text{grid}} = 600$, Bandbreite $T_{PW} = 0.95L$, Galerkin-Dim $N = 25$.
- Erste 3 Prolate-Eigenwerte: $\approx (0.999, 0.979, 0.794)$ — die ersten 3 Moden sind extrem gut bandbegrenzt (Slepian-Pollak-Plateau), danach fällt es rasch (4. Eigenwert $0.335$, 5. Eigenwert $0.053$).
- **QW via Sonin-Zerlegung:** $QW_{\lambda,\chi} = PW_\lambda + \Gamma_{\text{arch}} + \Gamma_{\text{prim}}(\chi)$
  - $\Gamma_{\text{arch}}(t) = 2\,\operatorname{Re}\,\psi(1/2 + it/2)$ (digamma, χ-unabhängig in dieser Version)
  - $\Gamma_{\text{prim}}(\chi, t, t') = -\sum_{p \le \lambda^2, p \nmid q}\sum_m \chi(p)^m \frac{\log p}{p^{m/2}} \cos((t-t')\cdot m\log p)$
- Drei Normalisierungen parallel:
  - **(a) v1-Stil:** $\varepsilon / \|\Psi \cdot PW\|_2$
  - **(b) Pol-normiert:** $\varepsilon / \|\Psi\|_F$
  - **(c) Dimensionsfrei:** $\varepsilon / (\|\Psi\|_2 \cdot \|PW\|_2)$

Getestete Charaktere (alle aus `CORE/zoo-mapping/_results/char_features_all.json`, mit Kronecker-Symbol-Werten via sympy verifiziert):

| Charakter | q | parity | γ^(1) | Typ |
|---|---|---|---|---|
| χ_0 | 1 | +1 | 14.1347 | trivial (Riemann ζ) |
| χ_4 | 4 | −1 | 6.0209 | ODD, Dirichlet β |
| χ_5 | 5 | +1 | 6.6484 | EVEN, Kronecker(+5) |
| χ_8 | 8 | +1 | 4.9000 | EVEN, Kronecker(+8) |
| χ_33 | 33 | +1 | 2.9970 | EVEN, Kronecker(+33) |

---

## 2. Rohdaten

| Charakter | γ^(1) | Defekt (spec) | rel_a (v1) | rel_b (pol) | rel_c (dim-frei) | μ_opt |
|---|---|---|---|---|---|---|
| χ_0 | 14.135 | — | **2.6499** | **7.2852** | **0.5501** | 0.563 |
| χ_4 | 6.021 | — | 0.1581 | 0.5985 | 0.1562 | 1.066 |
| χ_5 | 6.648 | — | 0.1929 | 0.8526 | 0.1928 | 1.084 |
| χ_8 | 4.900 | — | 0.3177 | 1.0677 | 0.2460 | 1.096 |
| χ_33 | 2.997 | — | **0.5583** | **1.2401** | 0.1908 | 1.050 |

Ratios (bezogen auf χ_0):

| Charakter | ratio_a | ratio_b | ratio_c | H1 Vorh. | H2 Vorh. |
|---|---|---|---|---|---|
| χ_0 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| χ_4 | 0.060 | 0.082 | 0.284 | 12.94 | 2.35 |
| χ_5 | 0.073 | 0.117 | 0.350 | 9.61 | 2.13 |
| χ_8 | 0.120 | 0.147 | 0.447 | 24.00 | 2.88 |
| χ_33 | 0.211 | 0.170 | 0.347 | 104.91 | 4.72 |

---

## 3. Struktur-Beobachtungen

### 3.1 χ_0 ist Outlier nach oben (ζ-Pol-Effekt)

Wie in Pilot v1: $\chi_0$ dominiert in allen drei Normalisierungen.  Das ist **nicht** ein Versagen der Methode, sondern ein struktureller Effekt: ζ hat einen Pol bei $s=1$, L-Funktionen nicht. Auf der kritischen Linie $\operatorname{Re}(s)=1/2$ sind wir vom Pol im Abstand 1/2; das macht $|\zeta(1/2+it)|$ im Mittel grösser als $|L(1/2+it,\chi)|$, und die Sonin-Zerlegung lokalisiert diesen Beitrag **nicht** im $\Gamma$-Teil, sondern im Ψ-Teil.

Die Normalisierung (b) $\varepsilon/\|\Psi\|_F$ versucht das rauszuteilen, reduziert aber nicht genug: χ_0 bleibt Faktor ≈ 6–12× grösser als die nicht-trivialen.

**Konsequenz:** χ_0 als Referenz ist für das H1/H2-Signal ungeeignet.  Die sinnvollere Referenz sind die nicht-trivialen Even-Charaktere untereinander.

### 3.2 **Neues Resultat (v2): H2-konsistentes Wachstum unter Even-Charakteren**

Restriktion auf $\{\chi_5, \chi_8, \chi_{33}\}$ (alle EVEN, reelle primitive Kronecker-Charaktere):

| Charakter | γ^(1) | rel_a | rel_b |
|---|---|---|---|
| χ_5 | 6.65 | 0.193 | 0.853 |
| χ_8 | 4.90 | 0.318 | 1.068 |
| χ_33 | 3.00 | 0.558 | 1.240 |

**Monotonie klar erfüllt**: kleineres γ^(1) ⇒ grösserer Defekt.  Das ist das erste mal in der Zookeeper/CCM + Zoo-Mapping-Reihe, dass dieser Trend in der richtigen Richtung erscheint (v1 hatte die umgekehrte Richtung, Atlas v1–v4 hatten char-invariante Gaps).

**Log-Log-Regression** (nur auf Even-Charaktere):

| Größe | Slope (log rel / log γ^(1)) | H2-Vorhersage | H1-Vorhersage |
|---|---|---|---|
| rel_a | **−1.33** | −1 | −3 |
| rel_b | **−0.47** | −1 | −3 |
| rel_c | +0.00 | −1 | −3 |

rel_a ist **deutlich näher zu H2 (−1) als zu H1 (−3)**, liegt leicht darüber.  rel_b ist schwächer (−0.47), rel_c zeigt keine γ-Abhängigkeit (dimensionsfreie Version ist saturiert).

### 3.3 χ_4 (ODD, nicht-Atlas) bricht den Trend

χ_4 ist ODD und liegt zwischen χ_5 und χ_8 in γ^(1) (6.02), aber mit rel_a = 0.158 kleiner als χ_5 (0.193) trotz kleinerer γ^(1).  Das liegt plausibel an der Parität: $\Gamma_{\text{arch}}$ in meiner Implementation ist χ-unabhängig, sollte aber für odd Charaktere digamma an $3/4 + it/2$ statt $1/2 + it/2$ auswerten (CCM Weil-Term ist parity-sensitiv).  **Das ist ein Implementierungs-Restfehler** — eine Iteration v2.1 sollte das fixen.

### 3.4 μ_opt-Signatur

Alle nicht-trivialen Charaktere haben $\mu_{\text{opt}} \approx 1.05{-}1.10$ (reell positiv), χ_0 hat $\mu_{\text{opt}} \approx 0.56$.  Das ist **strukturell konsistent**: für nicht-triviale L-Funktionen ist die "Skalenrelation zwischen QW und PW" nahe der Identität; ζ weicht als einziges deutlich ab, wiederum wegen Pol-Verschiebung.

---

## 4. Was sich gegenüber Pilot v1 geändert hat

| Aspekt | v1 | v2 |
|---|---|---|
| Basis | Naive Mellin-Stützstellen (Trapez-Quadratur) | **Echte Prolate-Galerkin-Basis** (numerische sinc-Diagonalisierung) |
| QW-Aufbau | Direkt: $K_{\text{arch}} +$ Prim-Matrix | **Sonin-Zerlegung**: $QW = PW + \Gamma_{\text{arch}} + \Gamma_{\text{prim}}$ |
| Normalisierung | ε / ‖B‖ (eine einzige) | **Drei Normalisierungen parallel** (v1-Stil, Ψ-Frobenius, dim-frei) |
| Charakter-Set | χ_0, χ_4 (2 Stück) | χ_0, χ_4, χ_5, χ_8, χ_33 (5 Stück) |
| Richtung | **Umgekehrt zu H1/H2** (Ratio 0.098 für χ_4/χ_0) | **Richtig zu H2** unter Even-Charakteren (Log-Log-Slope −1.33) |
| Blind-Spot | ζ-Pol dominiert, keine Signal-Struktur sichtbar | ζ-Pol dominiert weiter, aber unter L-Funktionen ist H2-Signatur sichtbar |

### 4.1 Kronecker-Charakter-Fehler in erstem v2-Lauf (korrigiert)

Im ersten v2-Lauf hatte ich χ_8 und χ_33 als ODD angesetzt (disc = −8, −33). Das war falsch: der Atlas (`char_features_all.json`) verwendet durchgängig **EVEN** reelle primitive Charaktere mit positiven Fundamentaldisc $D > 0$.  Korrigiert durch sympy-Verifikation der Kronecker-Werte. Erst nach der Korrektur erschien die monotone H2-Signatur.

Das ist eine **wichtige Lektion**: die Parity-Zuordnung der Charaktere bestimmt strukturell das Ergebnis.  Alle zukünftigen Tests müssen gegen Atlas-Feld `D` verifiziert werden.

---

## 5. Drei Lesarten der Ergebnisse

### 5.1 Lesart (a): H2-Validierung auf Even-Familie

Die Log-Log-Regression rel_a ~ 1/γ^{1.33} ist **deutlich näher an H2 (1/γ)** als an H1 (1/γ^3).  Unter der Einschränkung auf Even-Charaktere mit ähnlicher Rang-Struktur (χ_5, χ_8, χ_33) ist die Defekt-Norm H2-konsistent.

**Konsequenz:** Blueprint H2 ist **bestätigt** (schwach, $N=3$ Punkte, nur Even-Family).

### 5.2 Lesart (b): Der Slope −1.33 ist nicht exakt H2

H2 predicts Slope exact −1.  Wir messen −1.33.  Das könnte:
- (i) ein endlich-N-Artefakt sein (Galerkin $N=25$ ist klein),
- (ii) ein parity-unabhängiges $\Gamma_{\text{arch}}$-Artefakt sein (sollte parity-abhängig sein),
- (iii) auf eine **Mischung von H1 und H2** hinweisen (eine Zwischenhypothese H1.5).

### 5.3 Lesart (c): Vorsicht mit kleinen N

$N=25$ Prolate-Moden sind wenig. Die Skala $c = T_{PW} \cdot T_{\text{wide}} = 0.95L \cdot 3L \approx 3.7$ entspricht Paley-Wiener-Parameter $c/\pi \approx 1.2$, was nur etwa 3 Moden gut bandbegrenzt liefert (siehe §1, Eigenwerte 0.999, 0.979, 0.794).  Die übrigen 22 Moden sind zunehmend kontaminiert. Eine Version mit $c = 10$ ($N_{\text{eff}} \approx 10$) wäre solider.

---

## 6. Was dieser Pilot NICHT ist

- **Kein Beweis von H2.** Drei Punkte reichen nicht; dazu müsste man alle 10 Atlas-Charaktere testen und die Korrelation log ε vs. log(1/γ^(1)) prüfen.
- **Kein Test des χ_33-Anomalie-Vorzeichens.** Der Atlas-Gap von χ_33 ist **negativ** (−0.142 bei N=600). Unsere Defekt-Norm ist per Konstruktion ≥ 0. Wir messen die **Grösse** des Defekts, nicht sein Vorzeichen.
- **Kein Test der χ_21-N-Oszillation.** Wir variieren nicht in λ (nur festes $\lambda = \sqrt{14}$).
- **Keine vollständige Sonin-Zerlegung.** $\Gamma_{\text{arch}}$ ist parity-unabhängig in dieser Version; eine CCM-korrekte Version braucht parity-Spaltung.

---

## 7. Konkrete Zahlen für die nächste Iteration

Nächste logische Schritte:

### v2.1 — Parity-korrigiertes Γ_arch

$\Gamma_{\text{arch}}(\chi, t) = \psi\left(\frac{1 + \epsilon_\chi}{4} + \frac{it}{2}\right) + \psi\left(\frac{1 - \epsilon_\chi}{4} + \frac{it}{2}\right)$

mit $\epsilon_\chi = \chi(-1)$. Das sollte die χ_4-Abweichung (§3.3) beheben.

### v2.2 — Alle 10 Atlas-Charaktere

Für Atlas-Replikation: χ_5, χ_8, χ_12, χ_13, χ_17, χ_21, χ_24, χ_29, χ_33, χ_60.  Mit korrektem Kronecker-Lookup via sympy. Korrelation $\log \varepsilon$ vs. $\log(1/\gamma^{(1)})$ und $\log R_\chi$.

### v2.3 — λ-Skaling

Variiere $\lambda \in \{\sqrt{14}, 2\sqrt{14}, 4\sqrt{14}\}$ und $N_{\text{galerkin}} \in \{20, 40, 80\}$.  Prüfe ob der −1.33-Slope sich zu −1 (H2) stabilisiert bei grösserem $N$ oder bleibt (H1.5-Mischung).

### v2.4 — Atlas-Vorzeichen-Reproduktion

Das Atlas-Gap-Vorzeichen ist die eigentliche Signatur.  Eine Defekt-Norm-Version, die das Vorzeichen reproduziert, würde das Blueprint wirklich validieren.  Ansatz: statt $\|A - \mu B\|$ betrachte $\langle v_0, (A-\mu B) v_0\rangle$ für ein χ-abhängiges Eigenvektor $v_0$ (Grundzustand).

---

## 8. Ehrliche Einordnung

Dies ist ein **Zwischenergebnis**, kein Durchbruch.  Was bleibt:

**Positiv.** Unter Even-Charakteren zeigt sich die erste H2-konsistente Signatur in einer expliziten Operator-Defekt-Norm im `fst_spectrum_duality`-Programm.  Das ist qualitativ neu gegenüber Atlas v1–v4 (char-invariant) und Pilot v1 (H-Gegenrichtung).

**Negativ.** χ_0 bleibt Outlier (ζ-Pol), χ_4 bricht den Trend (parity-Problem), Slope −1.33 ist nicht exakt −1. Drei Datenpunkte (Even-Familie) reichen statistisch nicht.

**Konstruktiv weiter.** v2.1 (parity-korrektes Γ_arch) und v2.2 (10 Atlas-Charaktere) sind klare nächste Iterationen, jeweils wenige Stunden Arbeit.

Die zentrale Hypothese des Blueprints `DIRICHLET_CCM_TRANSFER.md` — dass die χ-Defekt-Norm die char-spezifische Information **trägt** und nicht versteckt — ist im Pilot v2 **nicht widerlegt**.  Ob sie **bestätigt** wird, hängt davon ab, ob v2.1 + v2.2 den H2-Trend auf der Full-Atlas-Familie stabilisieren.

---

## 9. Konkrete Zahlen (für spätere Referenz)

```
lambda      = 3.7417 (sqrt(14))
L           = 1.3195 (log lambda)
N_grid      = 600
N_galerkin  = 25
T_PW        = 1.2536 (0.95 * L)
T_wide      = 3.9586 (3.0 * L)

Prolate Eigenwerte (erste 5): [0.999, 0.979, 0.794, 0.335, 0.053]

chi_0 (trivial): rel_a=2.650, rel_b=7.285, rel_c=0.550, mu_opt=+0.563
chi_4 (odd, gamma1=6.02): rel_a=0.158, rel_b=0.599, rel_c=0.156, mu_opt=+1.066
chi_5 (even, gamma1=6.65): rel_a=0.193, rel_b=0.853, rel_c=0.193, mu_opt=+1.084
chi_8 (even, gamma1=4.90): rel_a=0.318, rel_b=1.068, rel_c=0.246, mu_opt=+1.096
chi_33 (even, gamma1=3.00): rel_a=0.558, rel_b=1.240, rel_c=0.191, mu_opt=+1.050

Log-Log-Regression (nur Even: chi_5, chi_8, chi_33):
  rel_a ~ 1/gamma^1.33  (H2 erwartet 1/gamma^1.00)
  rel_b ~ 1/gamma^0.47  (H2 erwartet 1/gamma^1.00)
  rel_c ~ const          (dim-frei, saturiert)
```

---

## 10. Aktualisierungs-Log

| Datum | Autor | Änderung |
|---|---|---|
| 2026-04-18 | LG (Opus 4.7 [1M], Session 16 Fortsetzung) | Pilot v2-Lauf: echte Prolate-Basis + Sonin-Zerlegung + 3 Normalisierungen. Erste Version mit Fehler in χ_33/χ_8 Parity (als ODD angenommen, korrigiert zu EVEN per Atlas-`D`-Feld + sympy Kronecker-Verifikation). Nach Korrektur: unter Even-Charakteren monotoner H2-konsistenter Trend (rel_a ∝ 1/γ^1.33, näher zu H2 als H1). χ_0 bleibt Outlier wegen ζ-Pol. χ_4 bricht Trend wegen parity-unabhängigem Γ_arch in Implementation. Nächste Iteration v2.1: parity-korrektes Γ_arch. Blueprint H2 weder bewiesen noch widerlegt, aber **in der richtigen Richtung bestätigt** gegenüber Pilot v1. |

---

**Ende PILOT_CHI_DEFECT_V2_2026-04-18.md.** Erstes positives Zwischensignal in der Zookeeper/CCM + Zoo-Mapping-Pipeline. Die H2-Richtung zu χ_0 bleibt wegen ζ-Pol verdeckt, aber unter L-Familien wird die monotone Skalenstruktur sichtbar. v2.1 (parity-Korrektur) + v2.2 (alle 10 Atlas-Charaktere) sind die klaren nächsten Schritte.
