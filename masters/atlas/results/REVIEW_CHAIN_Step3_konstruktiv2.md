# REVIEW CHAIN — Step 3: Konstruktiver Reviewer (zweite Runde)
## Paper: Dirichlet Character Atlas — A Weil-Kernel Cartography of Low-Conductor Sector Gaps (v2, Final Draft)
## Datum: 2026-04-29

---

## Abschnitt A: Was hat sich verbessert?

**E.1 — Parität-Sign-Note nach Gl. (5).** Umgesetzt in Zeilen 471–477. Die Note erklärt die `np.correlate`-Konvention, den Vorzeichen-Flip im odd sector und die numerische Gleichheit von Konvention und mathematischer Differenz. **Problem teilweise gelöst, ein Rest:** Step 2 nannte drei Optionen; das Paper hat *Option 3* (Note + Disclaimer) gewählt. Step 2 selbst kennzeichnete Option 3 als „akzeptabel für Preprint, blocking für Journal" — dieser Status ist also unverändert: Gl. (5) definiert weiterhin formal eine *Differenz*, was numerisch eine *Summe* ist. Für ein Journal muss das nachgezogen werden.

Zusätzlich: Der Einschub Z. 471–477 unterbricht den syntaktischen Fluss. Z. 470 endet mit Komma, Z. 479 beginnt mit „and the archimedean sector difference", was an die Definition vor dem Einschub anschließen will. **Empfehlung:** Note nach unten verschieben (nach beiden Definitionen, vor „Substituting into …") oder als `\begin{remark}` ausweisen.

**E.2 — Beweis Proposition 6.1 + Tone-Korrektur.** Umgesetzt in Z. 913–922 (Beweis) und Z. 924–927 (methodologische Konsequenz). Vollständig gelöst. **Rest (optional):** Die Verallgemeinerung auf alle Ritz-Galerkin-Verfahren (Step-1-C1) fehlt weiterhin — wäre die wertvollste optionale Aufwertung (2–3 Zeilen).

**E.3 — Tabelle 2 post-hoc + Tabelle 1 Caption + Abstract zwei Sign-Failures.** Alle drei Teilschritte umgesetzt und vollständig gelöst.

**E.4 — §7.3 umgeschrieben.** Umgesetzt. Titel ist „Sign convention: from empirical fix to rigorous derivation". Session-7-Narrativ entfernt. Vollständig gelöst.

**E.5 — Vier-Achsen-Präzisierung.** In §4 Einleitung (Z. 720–725) vollständig umgesetzt. **Verbleibender Rest:** Abstract Z. 95–96 enthält noch die *alte* Formulierung mit Klammer-Inhalt „(parity, root number, archimedean factor, gap signature)". E.5 wurde im Abstract nicht bis zum Ende durchgezogen. **Konkrete Korrektur (Pflicht vor arXiv):** Z. 95–96 ersetzen durch „We catalog the ten characters along two filter constants (parity and root number, both $+1$ in the sample) and two informative axes (archimedean factor weight, gap signature)".

**Zwischenfazit Abschnitt A.** Vier von fünf E-Punkten sind kompromisslos abgeschlossen. E.1 ist in der schwachen Variante (Preprint-tauglich, nicht Journal-tauglich). E.5 ist im Abstract noch nicht vollständig durchgezogen.

---

## Abschnitt B: Verbleibende Vorschläge (neue Perspektive)

**B.1 (NEU) — Tonaler Bruch in F1–F4 nicht aufgelöst.**
F1–F3 sind diagnostisch-negativ; F4 ist strukturell-positiv. Ein externer Reviewer wird fragen: Ist das Paper diagnostisch oder lieferst Du strukturelle Theoreme? **Vorschlag:** F4 entweder als Setup-Statement umformulieren und aus dem Findings-Block nach §III verschieben, oder den Abstract-Anspruch auf „diagnostic in the predictive register; structural in the algebraic register" präzisieren.

**B.2 (NEU) — N=400-Daten fehlen in Tabelle 1, sind aber in Tabelle 4 vorhanden.**
§3.1 nennt drei Galerkin-Resolutionen $N\in\{200,400,600\}$; Tabelle 1 zeigt nur N=200 und N=600. Inkonsistenz. Entweder N=400-Spalte ergänzen oder explizit begründen, warum sie weggelassen wurde.

**B.3 (NEU) — Title-Präfix „FST-Mathematics:" ist internes Jargon.**
Für eine Journal-Submission empfehle ich das Präfix zu entfernen. **Vorschlag:** „The Dirichlet Character Atlas: A Weil-Kernel Cartography of Low-Conductor Sector Gaps". Der FST-Bezug bleibt im Abstract erhalten.

**B.4 (NEU) — §7.3 und Appendix B doppeln sich.**
Beide Stellen erklären die Connes-vs-Weil-Q-Konventionen mit denselben sechs Sign-Sources. §7.3 (17 Zeilen) und Appendix B (51 Zeilen) sind redundant. **Empfehlung:** §7.3 auf 4–5 Zeilen kürzen.

**B.5 (Aus Step 2 nicht umgesetzt) — Spectral-Convergence-Verweis fehlt.**
Babuška–Osborn-Verweis in §3 fehlt. Nach Lektüre des überarbeiteten Papers ist klar: §6 (χ_21) sagt „der Galerkin-Operator selbst hat nicht konvergiert" ohne den Konvergenz-Begriff präzise zu definieren. Ein Verweis (Babuška–Osborn 1991 oder Trefethen–Bau 1997) im Setup würde die χ_21-Diskussion absichern.

**B.6 (Aus Step 2 nicht umgesetzt) — Theoreme 3.2/3.3 ohne Beweisskizzen.**
Beide zentralen strukturellen Theoreme verweisen auf [AnalyticKernel] als interne Note. Ein Reviewer ohne Kenntnis von [AnalyticKernel] kann das Paper nicht selbsttragend lesen.

**B.7 (Aus Step 1 nicht umgesetzt) — Selektionskriterium der zehn Charaktere bleibt implizit.**
$D \in \{5,8,12,13,17,21,24,29,33,60\}$ ohne explizite Begründung. Ein Satz im Setup würde Reviewer-Rückfragen verhindern.

**B.8 (Aus Step 1 nicht umgesetzt) — λ=20000-Begründung fehlt.**
Reviewer-Standardfrage. Ein Satz: „chosen as the largest cut-off for which N=600 Galerkin diagonalization remains feasible within ~2 hours wall-clock per ten-character batch."

**B.9 (Aus Step 1 nicht umgesetzt) — Interne Supplementary-Bibitems unverifizierbar.**
AnalyticKernel, AnalyticKernelV2 usw. verweisen auf `.md`-Pfade. **Konkrete Optionen:** (a) Als Zenodo-Supplementary mit eigenen DOIs hochladen; (b) Essentielle Inhalte als Anhänge ins Atlas-Paper integrieren.

**B.10 (NEU) — χ_21 Quasi-Degeneracy bleibt narrativ, nicht spektral belegt.**
Die mechanistisch-plausible Erklärung ist ohne spektralen Beleg (Tabelle mit $\lambda_2 - \lambda_1$ und Eigenvektor-Überlappungen für N=200/400/600). Für ein Journal reale Reviewer-Angriffsfläche.

**B.11 (NEU, kosmetisch) — Zwei Ungenauigkeiten in §3.1.**
- Z. 568: „Fast Fourier transform of $\rho_{\chifn}$" ohne Rückverweis auf die Definition in Z. 313.
- Z. 575: Notation für Sektor-Grundzustand inkonsistent ($\phi^{(N)}$ vs. $\phi^{\pm,(N)}$ vs. $\phi^\pm_{\chifn,N}$).

---

## Abschnitt C: Journal-Empfehlung

**Erstwahl: Experimental Mathematics (Taylor & Francis).**
Scope: „Computational and experimental mathematical research" — explizit für numerische Exploration mit Diagnose-Charakter. Akzeptiert ehrliche Negativ-Resultate und Methoden-Kritik. Erwartet computational reproducibility.

**Zweitwahl: Research in Number Theory (Springer Nature).**
Open Access. Breiter Scope inklusive computational und diagnostischer Arbeiten.

**Drittwahl: Mathematics of Computation (AMS).**
Strengerer numerisch-rigoroser Anspruch; schwieriger mit einem Paper, dessen zentrale Aussage „die Methode liefert keine prädiktive Kraft" ist.

**Nicht empfohlen:** *Journal of Number Theory* (erwartet positive Theoreme), *Acta Arithmetica*, *Inventiones / Annals*.

**Strategische Empfehlung:** arXiv-Preprint zuerst (math.NT + Cross-listing math.SP), dann Submission zu *Experimental Mathematics*. Falls Reject, *Research in Number Theory*.

---

## Abschnitt D: Narrative Analyse

**Der Spannungsbogen** funktioniert jetzt nach den E-Fixes:

1. *Setup* (§II–§III): Rigoroser Weil-Kernel-Apparat, vier Theoreme.
2. *Versuch und Falsifikation* (§III.B): Die natürliche Reduktion $\phi^+\approx\phi^-$ wird empirisch zerstört.
3. *Atlas* (§IV–§V): Was bleibt, ist eine Klassifikation entlang der intakten Achsen.
4. *Selbst-Demaskierung* (§VII): Die scheinbare 9/9-Konvergenz bei N=200 wird als Tautologie entlarvt.
5. *Anomalie-Analyse* (§VIII): χ_21 als N-Oszillator.
6. *Diagnose* (§IX): Hindernis = fehlende analytische Kontrolle der Sektor-Grundzustände.
7. *Programm* (§X): Sechs konkrete Wege, Route D (CCM-Twist) als Hauptpfad.

**Wo bricht die Narrative ab:**
- *F4 versus F1–F3* (wie in B.1 beschrieben).
- *§7.2 (Slope-0.72-Diskussion)* liest sich wie eine nachgereichte Verteidigung. Entweder stärker ins Hauptteil-Narrativ einbinden oder kürzen.
- *Übergang §VII → §VIII* zu abrupt. Brücke: „Beyond the systematic tautology of §VII, one character demands a separate spectral analysis: $\chi_{21}$."

**Abstract-Schwachstelle (dringend):**
Z. 95–96 noch mit alter Vier-Achsen-Formulierung. Korrektur: „We catalog the ten characters along two filter constants (parity and root number, both $+1$ in the sample) and two informative axes (archimedean factor weight, gap signature)."

**Stärken der jetzigen Narrative:**
- Self-deprecating tone funktioniert und macht das Paper vertrauenswürdig.
- Drei Galerkin-Resolutionen sind ein echtes methodologisches Argument, das im Setup noch nicht explizit verkauft wird.
- Trichotomie prime-dominated/mixed/arch-dominated ist die solideste positive Lieferung.

---

## Abschnitt E: Readiness-Score

**arXiv-Preprint: 8/10**
Begründung: Alle fünf Pflicht-Fixes aus Step 2 sind umgesetzt. Drei Mikro-Fixes vor arXiv: (1) Title-Präfix entfernen, (2) Abstract Z. 95–96 E.5-Korrektur, (3) Optional: Note-Block-Position Z. 471–477 verbessern.

**Journal-Einreichung: 7/10**
Blocking für Journal: (a) E.1-Lösung in der schwachen Variante (Definition vs. Implementation), (b) Theoreme 3.2/3.3 ohne Beweisskizzen, (c) Interne Supplementary-Bibitems unverifizierbar. Ferner: B.5 (Spectral-Convergence), B.10 (χ_21 spektral belegen). Mit diesen Fixes wäre das Paper bei 8.5/10 Journal-fähig.

**Empfehlung an die Pipeline:**
1. **Drei Mikro-Fixes vor arXiv** (Title, Abstract, optional Note-Block).
2. arXiv-Submission an math.NT + Cross-listing math.SP.
3. **Revision 2:** Fokus auf B.5, B.6, B.7, B.8, B.9, B.10.
4. Dann *Experimental Mathematics* einreichen.
5. B.1, B.2, B.3, B.4 als POSTRELEASE-PR-Items festhalten.

---

## Belegstellen für die Pflicht-Fixes (Quick Reference)

| Fix | Paper-Stelle | Status |
|---|---|---|
| E.1 (Parität-Note) | Z. 471–477 | Umgesetzt (Option 3); Journal-Rest offen |
| E.2 (Beweis Prop 6.1) | Z. 913–922 | Umgesetzt; C1-Verallgemeinerung optional offen |
| E.2 (Tone-Korrektur) | Z. 924–927 | Umgesetzt |
| E.3 (Tab 1 Caption) | Z. 686–691 | Umgesetzt |
| E.3 (Tab 2 Caption) | Z. 624–627 | Umgesetzt |
| E.3 (Abstract zwei Sign-Failures) | Z. 197–200 | Umgesetzt |
| E.4 (§7.3 umgeschrieben) | Z. 1123–1139 | Umgesetzt |
| E.5 (§4 Vier-Achsen) | Z. 720–725 | Umgesetzt |
| E.5 (Abstract Vier-Achsen) | Z. 95–96 | **OFFEN: Klammer-Inhalt noch alte Formulierung** |
