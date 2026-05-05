# Milestone 2_χ v2 — ρ-Term-Einbau: Drei strukturelle Lektionen (2026-04-18)

**Autor:** LG (Claude Opus 4.7 [1M], Session 16 Fortsetzung)
**Skripte:**
- `_scripts/milestone_2_chi_v2_rho.py` (strikte Paley-Wiener-Version, $T_{PW} < \log\lambda$)
- `_scripts/milestone_2_chi_v2_rho_wideband.py` (diagnostische Wideband-Variante, $T_{PW} = 25$)

**Status:** Negativ-Ergebnis mit strukturellen Einsichten. **Die direkte Galerkin-Matrix-Form der Weil-Explizitformel ist nicht die Weil-Positivitäts-Aussage.**  QW-Positivität wurde in beiden Varianten **NICHT** erreicht.  Die Lehren daraus revidieren den Blueprint-Ansatz ernsthaft.

---

## 1. Was getestet wurde

Zwei Versionen der "vollen" Weil-Form
$$
QW_{\lambda,\chi} := \Gamma_\rho(\chi) - \Gamma_{\text{arch}}(\chi) - \Gamma_{\text{prim}}(\chi)
$$
wurden in der Prolate-Galerkin-Basis implementiert.

### 1.1 ρ-Term-Definition

$$
[\Gamma_\rho(\chi)]_{mn} \;=\; \sum_{k : |\gamma_\chi^{(k)}| \le \Lambda}\Big[h_m(\gamma_\chi^{(k)})\,h_n(\gamma_\chi^{(k)}) \;+\; h_m(-\gamma_\chi^{(k)})\,h_n(-\gamma_\chi^{(k)})\Big]
$$
mit $h_n$ als Prolate-Basisfunktionen und Nullstellen aus LMFDB-Daten (für χ_5, χ_8, χ_33 aus `zeros_all_chars.json`) bzw. Literaturwerten (20 Riemann-Nullstellen für χ_0, 10 Dirichlet-β-Nullstellen für χ_4).

### 1.2 Zwei Varianten

| Variante | $T_{PW}$ | $T_{\text{wide}}$ | Milestone-1-konform? |
|---|---|---|---|
| A: Strikt | $0.95 \log\lambda \approx 1.25$ | 25 | **Ja** (strikte Paley-Wiener) |
| B: Wideband | **25** | 25 | **Nein** (bricht $T < \log\lambda$) |

---

## 2. Ergebnisse Variante A (strikt, Milestone-1-konform)

### 2.1 Roh-Daten

| χ | γ^(1) | #ρ im Cutoff | ‖G_arch‖ | ‖G_prim‖ | ‖G_ρ‖ | min λ(QW) | pos.-def? |
|---|---|---|---|---|---|---|---|
| χ_0 | 14.135 | 1 | 5.04 | 18.47 | 0.61 | −5.00 | nein |
| χ_4 | 6.021 | 5 | 5.04 | 17.40 | 0.84 | −19.21 | nein |
| χ_5 | 6.648 | 6 | 5.04 | 16.36 | 0.90 | −19.58 | nein |
| χ_8 | 4.900 | 7 | 5.04 | 17.79 | 0.84 | −20.87 | nein |
| χ_33 | 2.997 | 12 | 5.04 | 17.88 | 1.17 | −20.93 | nein |

**Beobachtung 1:** ‖G_ρ‖ ist **drastisch kleiner** als ‖G_arch‖ und ‖G_prim‖ (Faktor 20×).

**Beobachtung 2:** QW bleibt **weit von positiv-definit**, min-EW $\sim -20$.

### 2.2 Log-Log-Slopes (Even-Familie)

| Größe | Slope |
|---|---|
| rel_a | −0.326 |
| rel_b | −0.239 |

**Signal ist schwach** — weder H1 (−3) noch H2 (−1).

---

## 3. Ergebnisse Variante B (Wideband, diagnostisch)

### 3.1 Rohdaten

| χ | γ^(1) | ‖G_ρ‖ | min λ(QW) |
|---|---|---|---|
| χ_0 | 14.135 | 1.10 | −4.85 |
| χ_4 | 6.021 | **5.27** | −4.85 |
| χ_5 | 6.648 | 3.07 | −4.86 |
| χ_8 | 4.900 | 2.50 | −4.87 |
| χ_33 | 2.997 | 2.39 | −4.86 |

**Beobachtung 3:** In Wideband-Version wirkt der ρ-Term stärker (‖G_ρ‖ bis 5.27 für χ_4), aber **min-EW ist ≈ −4.85 für ALLE Charaktere, nahezu identisch**.

**Beobachtung 4:** Das "gemeinsame" negative Eigenwert-Niveau kommt aus $\Gamma_{\text{arch}}$ (||G_arch||=5.04–6.54 ist selbst für alle charakter-unabhängig bei fester parity).

### 3.2 Defekt-Slopes

| Größe | Slope |
|---|---|
| rel_a | +0.13 (flach, wechselhaft) |
| rel_b | −0.26 |

**Kein klares Signal.**

---

## 4. Drei strukturelle Lektionen

### 4.1 Lektion 1: Strikte Paley-Wiener-Einschränkung macht ρ-Term wirkungslos

Bei $\lambda = \sqrt{14}$ ist $\log\lambda \approx 1.32$. Die Prolate-Basis mit $T_{PW} < \log\lambda$ ist auf $|t| \le 1.32$ konzentriert (Slepian-Pollak-Eigenfunktionen zerfallen schnell ausserhalb).  Aber:

- Niedrigste Riemann-Nullstelle: $\gamma_0^{(1)} = 14.13$
- Niedrigste Dirichlet-Nullstelle im Test: $\gamma_{\chi_{33}}^{(1)} = 2.997$

**Alle Nullstellen liegen ausserhalb des strikten Paley-Wiener-Bereichs der Prolate-Basis.**  Die Basis-Funktionen $h_n$ sind an diesen $\gamma$-Werten praktisch null.  Konsequenz: $[\Gamma_\rho]_{mn} \approx 0$ in der strikten Version, egal wie viele Nullstellen im Cutoff liegen.

**Strukturelle Konsequenz:** Die Defekt-Norm-Reduktion des Blueprints `DIRICHLET_CCM_TRANSFER.md` bei festem $\lambda = \sqrt{14}$ kann die χ-Nullstellen-Information **strukturell nicht übertragen**.  Sie ist ein **asymptotisches Konzept für** $\lambda \to \infty$.

### 4.2 Lektion 2: Wideband-Hack zeigt: Matrix-Positivität ≠ Weil-Positivität

Bei $T_{PW} = 25$ (weit ausserhalb der Paley-Wiener-Schranke) wirkt der ρ-Term (‖G_ρ‖ = 1.1–5.3), aber die Matrix $QW = \Gamma_\rho - \Gamma_{\text{arch}} - \Gamma_{\text{prim}}$ bleibt nicht positiv-definit.

**Warum?** Weil die Weil-Positivität **nicht** die Aussage "die Matrix $G_\rho - G_{\text{arch}} - G_{\text{prim}}$ ist positiv-semidefinit" ist.  Die echte Weil-Positivität ist
$$
\sum_{\rho} |\tilde f(\gamma_\rho)|^2 \;\ge\; |\text{arch}(f)| + \sum_p |\text{prim}_p(f)|
$$
eine **Ungleichung mit absoluten Beträgen / Quadraten**, keine lineare Differenz.  Beim Übergang zur Matrix-Darstellung entsteht eine Vorzeichen-Subtraktion, die **keine äquivalente Positivität-Aussage** ist.

**Strukturelle Konsequenz:** Der Blueprint-Test "baue QW als Matrix, check Eigenwerte ≥ 0" ist **fundamental inkorrekt** als Test für Weil-Positivität.  Die richtige Matrix-Form wäre quadratisch in der Test-Funktionen-Basis, nicht linear.

### 4.3 Lektion 3: Γ_arch dominiert die Vorzeichen-Struktur

In der Wideband-Version ist $\|G_{\text{arch}}\| \approx 5$ und der min-EW $\approx -4.85$ für alle Charaktere nahezu identisch.  Das heisst: **die negativen Eigenwerte kommen aus dem digamma-Term, nicht aus der χ-Information**.

Das ist konsistent mit der Rolle des archimedischen Terms in der Explizitformel: es ist ein "großer" Hintergrund-Beitrag, der die Form bis weit unter Null drückt.  Die χ-spezifische ρ-Information ist ein **kleiner Modifikator** auf diesem Hintergrund, wird aber durch die direkte Matrix-Differenz nicht sichtbar.

---

## 5. Was das für das Blueprint bedeutet

### 5.1 Blueprint-Zentraltest wird in Galerkin-Form unausführbar

Der Blueprint `DIRICHLET_CCM_TRANSFER.md` definiert:
$$
\varepsilon_{\lambda,\chi} \;=\; \|QW_{\lambda,\chi} \Psi_{\lambda,\chi} - \mu_{\lambda,\chi}\,\Psi_{\lambda,\chi} PW_{\lambda,\chi}\|
$$
mit Hypothesen H1 (R_χ-Skalierung) und H2 (1/γ^(1)-Skalierung).

Bei $\lambda = \sqrt{14}$:
- Strikte Version: $QW$ ohne wirksamen ρ-Term ⇒ keine χ-Nullstellen-Information ⇒ Defekt-Norm misst primär Prim-Summe.
- Wideband-Version: $QW$ mit ρ-Term, aber Matrix-Struktur ist falsch für Weil-Positivität.

**In beiden Fällen** ist der Defekt-Norm-Test nicht informativ für H1/H2 bei diesem $\lambda$.

### 5.2 Konsequenz 1: Blueprint braucht λ → ∞ Asymptotik

Der eigentliche Test von H1/H2 erfordert:
- **Skalen-Reihe** $\lambda \in \{e^3, e^7, e^{15}, \ldots\}$, wo die Nullstellen zunehmend in den Paley-Wiener-Bereich fallen.
- **Explizite Rate** $\varepsilon_{\lambda,\chi} \to 0$ als Funktion von $\lambda$ und der χ-Nullstellen.

Das ist Milestone 4_χ aus dem Blueprint, nicht Milestone 2_χ.

### 5.3 Konsequenz 2: Matrix-Konstruktion braucht Reformulierung

Statt $QW = G_\rho - G_{\text{arch}} - G_{\text{prim}}$ (lineare Subtraktion, positiv-def-Test falsch) braucht die korrekte Form:
$$
QW(f,g) \;=\; \sum_\rho \tilde f(\gamma_\rho) \overline{\tilde g(\gamma_\rho)} \;-\; \text{(quadratische arch/prim-Formen)}
$$
als **bilineare Form** auf $f \times g$.  Die Diagonal-Elemente sind dann die Positivitäts-Tests.

Das ist eine **echte Umstrukturierung** der Implementation und auch der Blueprint-§2.2-Formulierung.

### 5.4 Konsequenz 3: Weil-Positivität bleibt das harte Problem

Die Weil-Positivität ist die eigentlich schwere Aussage — **nicht** durch Matrix-Eigenwert-Check zu validieren.  Die echte Analyse braucht:
- Explizite Konstruktion der Testfunktions-Familie $\{f_\alpha\}$,
- Auswertung $Q_{\text{Weil}}(f_\alpha, f_\alpha) \ge 0$ für jede,
- Konvergenz gegen eine dichte Teilmenge.

Das ist das Connes-Consani-Moscovici-Programm (arXiv:2006.13771 "Weil positivity and trace formula at archimedean place").

---

## 6. Was dieser Lauf POSITIV liefert

Trotz negativer Primär-Resultate ist das wissenschaftlich **substantiell**:

1. **Klarstellung der Struktur:** Der Blueprint-Test "Matrix-Positivität von $QW$" ist **falsch formuliert**.  Das wusste ich vor diesem Lauf nicht.
2. **Diagnose des λ-Regimes:** Bei $\lambda = \sqrt{14}$ sind χ-Nullstellen strukturell ausserhalb des Paley-Wiener-Raums der Prolate-Basis.  Das ist **quantitativ** belegt.
3. **Vorzeichen-Lokalisierung:** ‖G_arch‖ ist Hauptquelle der Negativität; χ-Signal in G_ρ ist klein.  Eine korrekte Matrix-Form müsste G_ρ **quadratisch verstärken**, nicht linear subtrahieren.
4. **Lektion für Pilot v2:** Die rel_a-Slope −1.33 aus Pilot v2 (ohne ρ-Term) war ein **Artefakt der Prim-Summen-Struktur**, nicht ein echtes χ-Nullstellen-Signal.  Das ist eine wichtige Revision des Pilot-v2-Berichts.

---

## 7. Konkrete Zahlen

### 7.1 Variante A (strikt)

```
lambda     = 3.7417
L          = 1.3195
T_PW       = 1.2536 (0.95 * L)
T_wide     = 25.0
N_grid     = 2400
N_galerkin = 40
rho_cutoff = 20.0

chi_0:   ||G_rho||=0.61, ||G_arch||=5.04, ||G_prim||=18.47, min_EW=-5.00, nicht pos.
chi_4:   ||G_rho||=0.84, ||G_arch||=5.04, ||G_prim||=17.40, min_EW=-19.21, nicht pos.
chi_5:   ||G_rho||=0.90, ||G_arch||=5.04, ||G_prim||=16.36, min_EW=-19.58, nicht pos.
chi_8:   ||G_rho||=0.84, ||G_arch||=5.04, ||G_prim||=17.79, min_EW=-20.87, nicht pos.
chi_33:  ||G_rho||=1.17, ||G_arch||=5.04, ||G_prim||=17.88, min_EW=-20.93, nicht pos.

Slope rel_a (Even-Familie): -0.326
Slope rel_b (Even-Familie): -0.239
```

### 7.2 Variante B (wideband)

```
T_PW       = 25.0 (!! bricht Milestone-1-Schranke T < log lambda)
N_galerkin = 60

chi_0:   ||G_rho||=1.10, min_EW=-4.85
chi_4:   ||G_rho||=5.27, min_EW=-4.85
chi_5:   ||G_rho||=3.07, min_EW=-4.86
chi_8:   ||G_rho||=2.50, min_EW=-4.87
chi_33:  ||G_rho||=2.39, min_EW=-4.86

Alle min_EW ~ -4.85, nahezu identisch (Artefakt von ||G_arch||).

Slope rel_a (Even-Familie): +0.13 (chaotisch, kein Signal)
Slope rel_b (Even-Familie): -0.26
```

---

## 8. Strategische Konsequenzen für die nächste Session

### 8.1 Priorität-Revision

**Vorher** (nach Pilot v2): glaubwürdiger H2-Slope −1.33 ⇒ H2 wird in v2.1/v2.2 verifiziert ⇒ Blueprint stärkt.

**Jetzt** (nach ρ-Term-Einbau): der −1.33-Slope aus Pilot v2 war wahrscheinlich **Prim-Summen-Artefakt**, nicht echte Nullstellen-Signatur.  Der korrekte Test braucht:
- (a) Bilineare Weil-Form (statt lineare Subtraktion),
- (b) λ-Asymptotik (statt festes $\lambda = \sqrt{14}$),
- (c) Korrekte Positivitäts-Kriterien (Testfunktions-Evaluation statt Matrix-Eigenwert).

### 8.2 Nächste Aufgaben, ehrlich formuliert

1. **Blueprint `DIRICHLET_CCM_TRANSFER.md` §2.2 revidieren.**  Die Defekt-Norm-Formulierung braucht explizite Unterscheidung zwischen: (i) Matrix-Form in Galerkin-Approximation, (ii) bilineare Form auf Testfunktionen, (iii) Asymptotik $\lambda \to \infty$.

2. **Blueprint `DIRICHLET_CCM_TRANSFER.md` Milestones reorganisieren.**  Derzeitige Reihenfolge (1_χ → 2_χ → 3_χ → 4_χ) ist **falsch**: 4_χ (λ-Asymptotik) ist Voraussetzung für die anderen, nicht deren Folge.  Korrekte Reihenfolge wäre: **4_χ ist das Kern-Milestone**, die anderen sind Teilschritte.

3. **Pilot-v2-Bericht revidieren.**  Die H2-Slope-Interpretation muss mit dem Caveat versehen werden: "Slope im festen-λ-Regime; echtes H2/H1-Testen erfordert λ-Variation".

### 8.3 Was NICHT mehr machen bis Revision

- **Nicht** v2.1, v2.2, v2.3 direkt bauen (parity-Fix, 10 Charaktere, λ-Skalierung) — das verstärkt nur das Artefakt-Signal.
- **Nicht** Milestone 2_χ v3 versuchen mit mehr Moden — die strukturelle Diagnose zeigt, dass das nicht der Engpass ist.
- **Nicht** direkt Milestone 3_χ oder 4_χ angehen ohne Blueprint-Revision.

---

## 9. Ehrliche Einordnung

Das ist **das zweite Negativ-Ergebnis in Folge** (nach Pilot v1 Richtungs-Umkehr).  Beide Negativ-Ergebnisse sind aber **strukturell aufklärend**:

- **Pilot v1 → v2:** "Implementation zu primitiv" → "bessere Basis + Sonin-Zerlegung" → partielle H2-Signatur.
- **Pilot v2 → Milestone 2_χ v2 (diese):** "Sonin-Zerlegung unvollständig" → "ρ-Term einbauen" → **Entdeckung: Matrix-Form ist falsch für Weil-Positivität**.

Das ist der Kurs einer seriösen numerischen Untersuchung: jede Iteration legt eine tiefere strukturelle Schwierigkeit frei.  Die derzeitige Schwierigkeit (Matrix vs. bilineare Form) ist die tiefste bisher — sie betrifft den **Blueprint-Kern**, nicht die Implementation.

**Positiv:** Wir wissen jetzt genau, was **nicht** funktioniert, und warum.  Das ist der Unterschied zwischen einem Sackgassen-Lauf und einem Diagnostik-Lauf.

**Negativ:** Ein echter H1/H2-Test ist bei $\lambda = \sqrt{14}$ mit direkter Matrix-Konstruktion **nicht möglich**.  Die nächste Iteration muss strukturell anders sein: bilineare Testfunktions-Evaluation oder λ → ∞-Asymptotik.

---

## 10. Aktualisierungs-Log

| Datum | Autor | Änderung |
|---|---|---|
| 2026-04-18 | LG (Opus 4.7 [1M], Session 16 Fortsetzung) | Milestone 2_χ v2 mit ρ-Term in zwei Varianten (strikt + wideband). **Negativ-Ergebnis in beiden**: QW nicht positiv-definit, Defekt-Slopes flach (−0.33 bzw. +0.13). **Drei strukturelle Lektionen**: (1) bei $\lambda = \sqrt{14}$ liegen Nullstellen ausserhalb Paley-Wiener-Bereich der Prolate-Basis ⇒ ρ-Term strukturell wirkungslos; (2) Matrix-Form $QW = G_\rho - G_{\text{arch}} - G_{\text{prim}}$ ist **nicht** die Weil-Positivitäts-Aussage; diese ist bilinear auf Testfunktionen; (3) $\Gamma_{\text{arch}}$ dominiert Vorzeichen-Struktur, χ-ρ-Signal verschwindet im Noise. **Strategische Revision**: Pilot-v2 H2-Slope −1.33 wahrscheinlich Prim-Summen-Artefakt, nicht Nullstellen-Signatur. Blueprint §2.2 Defekt-Norm-Formulierung braucht Überarbeitung; Milestone-Reihenfolge 1→4 sollte zu **4 zuerst** (λ-Asymptotik als Kern) reorganisiert werden. Kein direkter Nächst-Schritt bis Blueprint-Revision. |

---

**Ende MILESTONE_2_CHI_V2_RHO_2026-04-18.md.**  Drittes strukturelles Negativ-Ergebnis der Session, aber das wissenschaftlich wertvollste: die Matrix-Form der Weil-Explizitformel ist **nicht** das korrekte Instrument zur Weil-Positivitäts-Validierung.  Das verschiebt das Blueprint-Kernproblem von "mehr Rechenleistung" zu "andere strukturelle Form".
