# REVIEW CHAIN — Step 2: Experte
## Paper: Dirichlet Character Atlas — A Weil-Kernel Cartography of Low-Conductor Sector Gaps (v1/v2)
## Datum: 2026-04-29

---

## A: Bewertung Step-1-Vorschläge

Step 1 hat eine solide Außensicht geliefert und die diagnostische Postur des Papers korrekt erkannt. Die B-Punkte zerfallen aber bei sachlich-fachlicher Prüfung in zwei Klassen, die Step 1 nicht klar trennt: **substanzkritisch** (mathematisch-inhaltlich) vs. **kosmetisch** (Präsentation, Lesefluss). Die folgende Priorisierung ordnet sie aus Experten-Sicht.

### A.1 Blocking (vor Zenodo-Preprint zwingend)

- **B4 (Widerspruch §7.3 ↔ Appendix B)** — **kritisch**. §7.3 sagt in Z. 1098–1115 sowohl „empirically fixed" und nennt die rigorous re-derivation in App B; das ist mathematisch in sich widersprüchlich (das Paper kann nicht gleichzeitig sagen „open" und „rigorously derived"). Die jetzt vorhandene Zeile 1113–1115 löst das eigentlich auf, aber die Formulierung „empirically fixed during Session 7 ... the rigorous re-derivation reconciling the two conventions is given in Appendix B" liest sich für einen externen Leser noch widersprüchlich. **Nachziehen:** „The coefficient $-2$ was originally fixed empirically during Session 7 of the development; in the present version it is rigorously derived in Appendix B and \cite{AnalyticKernelV2}, Theorem 5.1." Dadurch verschwindet der historische Geruch und die Aussage wird linear.

- **B1 (SKELETON-Header)** — gemäß meiner Lektüre bereits behoben in Z. 7 („Final Draft"). Step 1 hat eine ältere Version gesehen oder die Datei wurde nach Step 1 geändert. **Verifizieren** vor Submission, sonst wirklich kosmetischer Schaden.

- **B3 (Abstract erklärt nur einen Sign-Flip)** — **kritisch**. Der Abstract verschweigt $\chi_{21}$. Tabelle 1 markiert beide Charaktere $\chi_{21}$ und $\chi_{29}$ mit (F). Ein sorgfältiger Reviewer wird fragen „warum F2 nur einen sign failure benennt, wenn die zentrale Tabelle zwei zeigt". Die Argumentation in §6 (Tautologie + $\chi_{21}$-Oszillator) ist sinnvoll — aber der Abstract muss beide Fälle aufzählen. Vorschlag wie in Step 1: „two sign failures: $\chi_{21}$ ($N$-oscillator, §6) and $\chi_{29}$ (numerical interpolation artifact at near-zero gap)".

- **B6 (R²=0.99 inflationär)** — **kritisch**. Tabelle 2 zeigt $R^2=0.99$ für die „$\gap_{\mathrm{gal}}$ (all 10) at N=600"-Zeile, was per Tautologie wahr ist (sechs der zehn Einträge sind identisch mit der Empirie-Spalte). Ohne explizite Fußnote wird das missgelesen. **Nachziehen:** Zeile 631 in Tabelle 2 mit Fußnote „identity for 6/10 characters by construction (cf. §VI)" markieren.

### A.2 Wichtig, aber nicht blocking (Revision-2 ausreichend)

- **B2 (v1/v2 Mismatch)** — Datei `v1_en.tex`, `\date` sagt „v2". Konsistent ziehen. Ich präferiere v2 (entspricht dem Inhalt nach N=600-Lauf).

- **B5 (Root-Number-Achse degeneriert)** — Step 1 ist hier richtig, **aber unvollständig**. Auch die Parity-Achse (§4.1) ist konstant im Sample (alle Charaktere sind even, $\chi(-1)=+1$). Das Paper hat realistisch *zwei informative Achsen* (archimedean layer, gap signature) und *zwei Filter-Konstanten* (parity, root number). Das ist ehrlich gesagt in §4.1 und §4.2, aber die Überschrift „Character atlas in four axes" suggeriert Vier-Achsen-Information. **Vorschlag:** §4 mit „four-axis taxonomy (two filter constants, two informative axes)" einleiten. Das ist präziser und reviewer-resistent.

- **B7 (Selektionskriterium der zehn Charaktere)** — **wichtig** für Außensicht. Begründung in einem Satz: „We select the ten smallest fundamental discriminants $D > 0$ that yield primitive real even characters with $\chi(-1)=+1$ up to $D \leq 60$, providing a small but structurally diverse sample."

- **B8 (λ=20000 Begründung)** — **wichtig**. Standard-Reviewer-Frage. Ein Satz reicht.

### A.3 Optional / nachgeordnet

- **B9 (Internal Supplementary)** — Step 1 ist hier zurecht skeptisch, aber unter den FST-Programm-Gegebenheiten zumutbar. Lösung: Bibitems mit „Internal supplementary note, on file with the author" oder als Zenodo-Companion-Upload. Nicht blocking für ein Preprint, aber blocking für ein Journal.

### A.4 Kreativ-Punkte (C1–C5)

- **C1 (Tautologie als allgemeines Resultat)** ist die wertvollste Anregung von Step 1 — als Remark unter Prop 6.1 in maximal drei Zeilen umsetzbar und steigert die Zitierbarkeit deutlich. **Einbauen.**
- **C3 (N=1000-Checkpunkt für $\chi_{21}$)** ist hochrelevant (siehe B.4 unten); ohne diesen Datenpunkt bleibt die Quasi-Degeneracy-Erklärung **plausibel, aber nicht gesichert**.
- C2, C4, C5 — schöne Anregungen, aber für die Erstversion verzichtbar.

---

## B: Neue Schwachstellen (Step 1 übersehen)

### B.1 [BLOCKING] Δ_χ-Definitionsdiskrepanz: numerische vs. mathematische Faltung

**Befund.** Eq. (5) (Z. 467–469, ebenso Eq. (15)) definiert formal:
$$\Delta_\chi(t) = (\phi^+_\chi * \phi^+_\chi)(t) - (\phi^-_\chi * \phi^-_\chi)(t)$$

Appendix B (Z. 1281–1297) räumt ein:
- Das Skript verwendet `np.correlate`, also $R_{\phi^\pm}(u) = \int \phi^\pm(s)\phi^\pm(s-u)\,ds$
- Für $\phi^+$ (gerade): $R_{\phi^+} = (\phi^+ * \phi^+)$
- Für $\phi^-$ (ungerade): $R_{\phi^-} = -(\phi^- * \phi^-)_{\mathrm{math}}$
- Die Skript-Größe $\Delta^{\mathrm{num}}_\chi := R_{\phi^+} - R_{\phi^-}$ entspricht mathematisch $(\phi^+ * \phi^+) + (\phi^- * \phi^-)$ — eine **Summe**, nicht eine Differenz

**Kritik.** Das Paper definiert in den Theoremen ein Objekt $\Delta_\chi$, dessen formale Form (Differenz) von der numerisch berechneten Form (Summe nach Vorzeichenkorrektur) **abweicht**. Die Rettung in App B („Equation (5) is to be read as an operational definition anchored in `np.correlate`") ist ein methodologisches Eingeständnis, das ein analytischer Referee als „die Formel im Theorem stimmt nicht mit dem überein, was berechnet wurde" markieren wird.

**Status.** Verbesserungsfähig (korrekt im numerischen Resultat, aber Definition und Implementation stimmen formal nicht überein).

**Fix-Optionen** (in Reihenfolge der Präferenz):
1. **Definition korrigieren.** Eq. (5) und Eq. (15) auf $\Delta_\chi(t) := (\phi^+ * \phi^+)(t) + (\phi^- * \phi^-)(t)$ ändern — dann ist Theorem 3.4 mit Koeffizient $-2$ formal kohärent. Das ist der saubere Weg.
2. **Beide Konventionen explizit nebeneinander.** Eq. (5) als „mathematical Δ" stehen lassen, Eq. (15a) als „operational $\widetilde\Delta_\chi := R_{\phi^+} - R_{\phi^-}$" definieren, und Theorem 3.4 in beiden Formen angeben.
3. **Status quo + dickerer Disclaimer.** Eine Bemerkung nach Eq. (5): „The convolution form (5) and the autocorrelation form used in §III.E and Appendix B differ by a parity-induced sign on the odd-sector term; both define the same numerical quantity. See Appendix B for the parity-sign analysis."

Option 1 oder 2 sind reviewer-resistent; Option 3 ist akzeptabel für ein Preprint, blocking für ein Journal.

### B.2 [BLOCKING] Proposition 6.1: kein Beweis, selbst-deklarierte Inhaltsleere

**Befund.** Z. 888–899 nennt eine Aussage „Proposition 6.1 (Tautology)". Direkt nach der Aussage steht Z. 901–902: „The content of this statement is *nil*: it is an identity arising directly from the linearity of expectation". Es folgt **kein Beweis**.

**Kritik.** Eine „Proposition" ohne Beweis und mit selbst-deklarierter Inhaltsleere ist terminologisch problematisch. Ein referee mit Sinn für Konventionen wird das aufgreifen. Es gibt hier zwei valide Sichtweisen, die das Paper nicht gleichzeitig behaupten kann:

(a) Es ist eine *triviale* Identität — dann ist „Proposition" zu hoch gegriffen, „Observation" oder „Identity" wäre korrekt.

(b) Es ist eine *nicht-triviale* methodologische Erkenntnis — dann verdient sie einen Beweis, mindestens eine Beweis-Zeile, plus eine Diskussion **warum** sie methodologisch wichtig ist (das ist sie!). Dann ist auch die Aussage „content is nil" missverständlich — der propositionale Inhalt ist trivial, der **methodologische** Inhalt ist substantiell.

**Status.** Verbesserungsfähig (die Mathematik ist korrekt, die Präsentation ist suboptimal).

**Fix.** Zwei-Zeilen-Beweis hinzufügen + Tone-Korrektur:

```latex
\begin{proof}
Linearity of the inner product gives
$\langle\phi^\pm, W^\pm\phi^\pm\rangle = \langle\phi^\pm, W^\pm_{\mathrm{arch}}\phi^\pm\rangle + \langle\phi^\pm, W^\pm_{\mathrm{prime}}\phi^\pm\rangle$
on each sector. Since $\phi^\pm$ are the ground states, the eigenvalues
satisfy $\lambda_1^{\min}(W^\pm) = \langle\phi^\pm, W^\pm\phi^\pm\rangle$.
Subtracting yields the claim. $\square$
\end{proof}
```

Die nachfolgende Bemerkung „content is nil" ändern zu: „Although the propositional content is a trivial identity, its **methodological consequence** is non-trivial: it shows that any apparent agreement between $\gap_{\mathrm{gal}}^{(N)}$ and $S+\archdiff$ at fixed $N$ is structural, not predictive."

Das macht aus einer formal-rhetorischen Schwäche eine inhaltliche Stärke. Idealerweise zusätzlich Step-1-C1 einbauen (Verallgemeinerung auf alle Ritz-Galerkin-Verfahren).

### B.3 [BLOCKING] Post-hoc-Selection bei „stable 9"

**Befund.** Tabelle 2 (Z. 612–636) berichtet wiederholt „Stable" als Kategorie, definiert in der Caption (Z. 616–618) als „excludes $\chi_{21}$ (a known $N$-oscillator)". Die Klassifikation als „$N$-oscillator" entstand jedoch **aus** den Daten (§6, vgl. Z. 207–209 in F3: „...the $9/9$ at $N=200$ was not surprising once $\chi_{21}$ was excluded: the excluded character was the very one whose $N$-dependence was not yet integrated into the test").

**Kritik.** Das ist post-hoc selection. Ein strenger Referee wird sagen: „Wenn man den schwierigsten Fall identifiziert und dann ausschließt, sind 9/9-Statistiken zirkulär."

Das Paper ist hier **bereits ehrlich** — F3 sagt das offen — aber die Tabelle 2 selbst ist nicht entsprechend gekennzeichnet. Ein Reviewer sieht zuerst die Tabelle.

**Status.** Verbesserungsfähig (keine versteckte Manipulation, aber die Selektionslogik muss in der Tabelle selbst stehen).

**Fix.** Caption von Tabelle 2 ergänzen: „The ‚stable' subset excludes $\chi_{21}$ \emph{post hoc}, after its $N$-oscillation was identified as a Galerkin-truncation artefact (see §VI). Statistics on the stable subset are therefore conditional on this identification." Damit ist die Ehrlichkeit, die im Text bereits steht, in der Tabelle abgesichert.

### B.4 [WICHTIG] Galerkin-Konvergenz bei N=600 ohne Stabilitätsargument

**Befund.** Das Paper macht starke Aussagen über $N=600$-Werte und behandelt sie als „empirische Referenz" ($\gap_{\mathrm{emp}}$ in Tabelle 1). Genau die $\chi_{21}$-Episode (N=200 stationär bei +0.282 für 200 *und* 400, dann Kollaps zu −0.004 bei 600) ist ein **Beweis dafür, dass scheinbare Konvergenz fiktiv sein kann**.

**Kritik.** Wie wissen wir, dass $N=600$ nicht selbst ein anderes Plateau ist? Das Paper sagt für $\chi_{21}$ in Remark 6.4 (Z. 1001–1013), dass „the right theoretical object is $\gap_{\mathrm{emp}}(\chi_{21})=-0.00423$ at $N=600$" — aber genau diese Aussage ist durch die $\chi_{21}$-eigene Geschichte erschüttert. Step 1 hat das nur als optionale Idee C3 markiert; aus Expertensicht ist es **kein optionaler Punkt**, sondern eine zentrale strukturelle Frage.

**Status.** Offen (analytisch nicht ausgeräumt, numerisch nicht ausgereizt).

**Fix-Optionen.**
- **Minimal:** Ein Satz in §3 mit Verweis auf Spectral-Convergence-Theorie für glatte symmetrische Integraloperatoren auf $L^2[-L,L]$ (Babuška–Osborn 1991, Trefethen–Bau 1997). Das adressiert Step-1-D5 und macht eine generische Konvergenzaussage.
- **Stark:** Ein einzelner $N=1000$-Lauf für $\chi_{21}$ (Step-1-C3), CCX13 ~30 Min Zusatzkosten. Wenn der Wert bei $-0.0042 \pm \varepsilon$ bleibt, ist die Quasi-Degeneracy-Erklärung **gesichert** statt nur plausibel. Wenn nicht, wäre das Paper neu zu konzeptualisieren — was wahrscheinlicher *vor* dem Preprint passieren sollte.
- **Optimal:** Beides.

Empfehlung: Für das Zenodo-Preprint reicht der minimale Fix; für eine Journal-Einreichung ist der starke Fix angezeigt.

### B.5 [WICHTIG] χ_21-Quasi-Degeneracy: spekulativ, nicht spektral nachgewiesen

**Befund.** §6 erklärt die χ_21-Anomalie als „quasi-degeneracy collapse": die Galerkin-trunkierten Sektor-Eigenwerte für $N \in \{200, 400\}$ entsprechen einem „pseudo-ground state", der nicht der wahre $N\to\infty$-Grundzustand ist. Das ist mechanistisch plausibel, aber im Paper **nicht spektral belegt**.

**Kritik.** Was ein quasi-degeneracy collapse spektral signiert, ist eine kleine Lücke zwischen dem niedrigsten und dem zweitniedrigsten Eigenwert. Das Paper berichtet weder $\lambda_2 - \lambda_1$ für $\chi_{21}$ bei N=200, 400, 600, noch eine Norm-Überlappung der Eigenvektoren $\langle\phi^{(200)} | \phi^{(400)}\rangle$ vs. $\langle\phi^{(400)} | \phi^{(600)}\rangle$. Ohne diese Daten ist die Erklärung **eine Story**, keine Diagnose.

**Status.** Offen (Erklärung plausibel, aber strukturell nicht abgesichert).

**Fix.** Eine kurze Tabelle oder Bemerkung in §6 mit:
- $\lambda_1, \lambda_2$ (zweiter Eigenwert) für $W^\pm_{\chi_{21}, N}$, $N=200, 400, 600$
- Eigenvektor-Überlappung $|\langle\phi^{\pm,(200)} | \phi^{\pm,(600)}\rangle|$

Wenn eine kleine $\lambda_2 - \lambda_1$-Lücke bei N=200, 400 sichtbar ist und die Überlappung zwischen N=200 und N=600 klein ist, ist die Quasi-Degeneracy-Erklärung **belegt**. Das ist eine geringfügige Erweiterung des bestehenden Datensatzes (Skript existiert bereits), die strukturell viel bringt.

### B.6 [MITTEL] Theoreme 3.2/3.3: algebraische Basis ausgelagert

**Befund.** Die zwei zentralen strukturellen Theoreme (Sektor-Kerne, Anti-Diagonal-Differenz) verweisen für ihre Beweise auf [AnalyticKernel], eine interne `.md`-Datei. Step 1 hat das als „Verifizierbarkeitsproblem" angesprochen, aber die strukturelle Konsequenz ist tiefer: **die gesamte algebraische Basis des Papers ruht auf einer nicht-öffentlichen Note.**

**Kritik.** Ein Self-Contained-Argument sollte zumindest eine Beweisskizze enthalten, sonst ist das Paper nicht standalone lesbar.

**Status.** Verbesserungsfähig (Mathematik vermutlich korrekt, aber das Paper macht es einem externen Leser unnötig schwer).

**Fix.** Beweisskizzen in zwei bis drei Zeilen je Theorem hinzufügen:
- **Thm. 3.2:** „Folgt aus $K_\chi(x,y) = K_\chi(-x,-y)$ durch Symmetrisierung mit den Projektoren $P^\pm = (I \pm \mathcal{P})/2$."
- **Thm. 3.3:** „Folgt aus 3.2 durch Subtraktion der ‚+’- und ‚−’-Versionen."

Das ist trivial, macht das Paper aber sofort selbsttragend.

### B.7 [MITTEL] Vier-Achsen-Klassifikation: Konsistenzcheck

**Befund.** §4 verspricht „four orthogonal axes". Tatsächlich:
- **Parity:** konstant +1 (alle Charaktere even) → Filter
- **Root number $W$:** konstant +1 → Filter
- **Archimedean layer $r_\chi$:** informativ (drei Klassen)
- **Gap signature:** informativ (sechs Klassen)

**Kritik.** Step 1 hat nur die Root-Number-Achse als „degenerate" markiert. Aus Expertensicht sind **zwei** der vier Achsen Filter-Konstanten, nicht informative Achsen. Das verfälscht die methodologische Selbstdarstellung.

**Status.** Verbesserungsfähig.

**Fix.** §4 in der Einleitung präzisieren: „We organize the table along four taxonomic dimensions: two filter constants (parity, root number), both equal to +1 across the entire sample, and two informative axes (archimedean weight ratio, gap signature)." Plus Anpassung der Abstract-Formulierung „We catalog the ten characters along four axes" → „... along two filter constants and two informative axes".

### B.8 [GERING] Mathematische Schreibweise: mehrere kleine Inkonsistenzen

- Z. 50: „Numerical Atlas" im Titel; im Inhalt heißt es konsistent „Galerkin Diagnostics" — beides ist OK, aber redundant.
- Z. 750: $r_\chi := |\archdiff_\chi^{(200)}|/|S_\chi^{(200)}|$ — der Index $(200)$ ist eine Definition-Konvention, sollte aber im umgebenden Text gesagt werden („at $N=200$, $\lambda=20000$").
- Z. 762: Tabelle 3 sortiert nicht offensichtlich nach $r_\chi$ (chi_29 hat 1.8% und steht oben, dann 7.4%, 21%, 63%, 64%, 13%, 17%, 22%, 25%, 356% — die Sortierung wechselt zwischen Subklassen). Konsistente Sortierung würde Lesbarkeit erhöhen.
- Z. 84/494: $\sum_{p,m}$ vs. $\sum_{\substack{p\nmid q\\m\geq 1}}$ — konsistent ausschreiben.

Status: Verbesserungsfähig (kosmetisch, aber kumulativ professionell wichtig).

---

## C: Strukturelle Stärken (Reviewer-Resistenz)

Aus Expertensicht macht dieses Paper **drei Dinge** besonders gut, die Reviewer-Angriffe vorwegnehmen:

### C.1 Diagnostische Postur ist konsequent durchgehalten

Das Paper deklariert im Abstract „our contribution is deliberately diagnostic" und hält das durch. F1–F4 widersprechen sich nicht, sondern bauen auf. Das ist **die wichtigste Reviewer-Resistenz-Eigenschaft**: ein Reviewer, der das Paper als „Negativ-Resultat" angreifen will, findet im Text bereits die Antwort, dass dies *Absicht* ist und der Wert genau im präzisen Charakterisieren der Grenze liegt.

### C.2 Tautologie-Argument ist self-deprecating und selbstkonsistent

Das Paper räumt ein, dass das frühere „9/9 sign accuracy at N=200" eine Artefakt-Beobachtung war. Das ist ungewöhnlich für ein wissenschaftliches Paper und entkräftet vorab den schärfsten möglichen Angriff („die Methode funktioniert ja gar nicht"). Ein erfahrener Referee wird diese Selbstkritik *positiv* werten, vorausgesetzt die formalen Defizite (B.1, B.2 oben) werden behoben.

### C.3 Drei-Komponenten-Dekomposition als operativer Rahmen

RH_X = GBC_X + C^(2)_X + Hurwitz_X wird nicht nur als Etikett, sondern als **diagnostisches Werkzeug** verwendet: jedes Ergebnis wird einer Komponente zugeordnet. Das gibt dem Atlas eine theoretische Kohärenz, die ihn von einer reinen Datenpräsentation unterscheidet. Reviewer, die solche Frameworks aus der RH-Literatur kennen (Connes, Berry-Keating, Iwaniec-Kowalski), werden das anerkennen.

### C.4 Sechs Sub-Klassen in zehn Charakteren — genuine taxonomische Substanz

Step 1 hat das richtig gesehen: die Klassifikation (strong-neg, strong-pos, moderate-neg, weak-pos, arch-dominated, oscillator) ist nicht eine ad-hoc Nachsortierung, sondern eine echte Trichotomie-erweitert-um-drei. Ein Referee, der „what is new" fragt, hat hier eine konkrete Antwort.

### C.5 Drei Galerkin-Resolutionen sind ein Stärken-Argument

Die Cross-Validation N=200 → N=400 → N=600 ist vorbildlich. Die χ_21-Anomalie würde bei einer Single-N-Studie unentdeckt bleiben; gerade weil das Paper drei Skalen vergleicht, kann es überhaupt diagnostische Aussagen machen. Das macht das Paper *robuster* als ein hypothetisches einzelnes-N-Paper, nicht schwächer.

---

## D: Readiness-Score

**Aktueller Stand:** **6/10** (übereinstimmend mit Step 1).

**Begründung der Punktabzüge** (im Vergleich zu Step 1 erweitert):

| Defizit | Klasse | Punkte |
|---|---|---|
| §7.3 ↔ App B Widerspruch (B4) | Blocking | −0.5 |
| Abstract-Unvollständigkeit F2 (B3) | Blocking | −0.5 |
| R²=0.99 Tautologie nicht gekennzeichnet (B6) | Blocking | −0.5 |
| **Δ_χ Definition vs. Implementation (B.1 neu)** | **Blocking** | **−0.5** |
| **Prop 6.1 ohne Beweis + „content is nil" (B.2 neu)** | **Blocking** | **−0.5** |
| **Post-hoc selection bei „stable 9" nicht gekennzeichnet (B.3 neu)** | Wichtig | −0.3 |
| Vier-Achsen-Überverkauf, parity konstant (B.7 neu, erweitert Step-1-B5) | Mittel | −0.2 |
| Galerkin-Konvergenz bei N=600 ohne Stabilitätsargument (B.4 neu) | Mittel | −0.2 |
| χ_21-Quasi-Degeneracy nicht spektral belegt (B.5 neu) | Mittel | −0.2 |
| Theoreme 3.2/3.3 ohne Beweisskizzen (B.6 neu) | Mittel | −0.2 |
| Selektionskriterium + λ-Begründung (B7+B8 Step 1) | Gering | −0.2 |
| Internal Supplementary unverifizierbar (B9 Step 1) | Gering | −0.2 |
| Cosmetic (Header, Version, Sortierung Tabellen) | Trivial | −0.2 |

**Erreichbare Werte mit verschiedenen Behebungs-Stufen:**

- **Mit nur Step-1-B-Punkten** (B1, B3, B4, B6 — also Behebung der originalen Blocker): **7.5/10** — preprint-ready.
- **Mit Step-1 + B.1 + B.2 + B.3 (also den drei neuen Blockern aus Expertensicht)**: **8.0/10** — preprint-strong und revision-1-ready.
- **Mit zusätzlich B.4 + B.5 + B.6 (Konsolidierung)**: **8.5/10** — Journal-Submission-fähig.
- **Mit allem inkl. C-Punkten (Tautologie als allgemeines Ergebnis, N=1000-Checkpunkt)**: **9/10** — strong submission.

**Empfehlung Zenodo-Preprint:** 8.0/10 ist das realistische und wünschenswerte Ziel. Die zusätzlichen 0.5 Punkte für ein Journal können in Revision 2 (post-Zenodo) nachgezogen werden.

---

## E: Minimaler Einreichungsplan (Top-5)

Die folgenden fünf Änderungen heben das Paper von ~6/10 auf 8/10. Sie adressieren die Blocker, die Step 1 erkannt hat (4 Punkte) und die zwei zusätzlichen Blocker aus Expertensicht. Geschätzter Aufwand: 2–3 Stunden Bearbeitung, kein Recompute nötig.

### E.1 Δ_χ-Definition kohärent machen [B.1, **MUSS**]

**Schritt 1:** Eq. (5) (Z. 467–469) auf $\Delta_\chi(t) := (\phi^+_\chi * \phi^+_\chi)(t) + (\phi^-_\chi * \phi^-_\chi)(t)$ ändern, mit kurzer Anmerkung „(parity-corrected convolution form; see Appendix B for derivation from `np.correlate`)".

**Alternative Schritt 1:** Eq. (5) belassen + Zeile dahinter: „where the autocorrelation form computed by the script (App B) differs from the convolution form by a parity-induced sign on the odd term and yields the same numerical value."

**Schritt 2:** Theorem 3.4 (Theorem closure) und Eq. (15) entsprechend an die gewählte Form anpassen, sodass Koeffizient $-2$ formal direkt aus der Definition folgt.

### E.2 Proposition 6.1 mit Beweis und Tone-Korrektur [B.2, **MUSS**]

**Schritt 1:** Direkt nach Z. 899 (Ende der Proposition) zwei Beweis-Zeilen einfügen:
```latex
\begin{proof}
Linearity of the inner product on the Galerkin space gives
$\langle\phi^\pm, W^\pm \phi^\pm\rangle = \langle\phi^\pm, W^\pm_{\mathrm{arch}}\phi^\pm\rangle + \langle\phi^\pm, W^\pm_{\mathrm{prime}}\phi^\pm\rangle$.
Since $\phi^\pm$ are eigenvectors with eigenvalue $\lambda_1^{\min}(W^\pm)$,
the gap $\gap_{\mathrm{gal}}^{(N)}$ decomposes as the difference of the
two right-hand sides. $\square$
\end{proof}
```

**Schritt 2:** Z. 901 „The content of this statement is *nil*" ersetzen durch:
„Although the propositional content is a trivial linearity identity, its *methodological consequence* is non-trivial: any apparent agreement between $\gap_{\mathrm{gal}}^{(N)}$ and $S_\chi + \archdiff_\chi$ at fixed $N$ is structural rather than predictive."

**Schritt 3 (optional, hochwertig):** Step-1-C1 als Remark einbauen — drei Zeilen, die die Aussage auf alle Ritz-Galerkin-Diskretisierungen verallgemeinern. Steigert die Zitierbarkeit substantiell.

### E.3 Drei Tabellen/Aussagen mit Tautologie-/Selektion-Caveats [B3+B6+B.3 Step 2]

**Schritt 1:** Tabelle 1 (Z. 700–704), Zeile „$\gap_{\mathrm{gal}}$ sign accuracy ($N=600$): $10/10$ ($R^2=0.99$, tautology)" — den expliziten Hinweis „tautology" so wie er steht **belassen**, aber in der Tabellen-Caption (Z. 678–681) ergänzen: „Six of ten $\gap_{\mathrm{emp}}$ entries (where $N_{\mathrm{src}}=600$) coincide identically with $\gap_{\mathrm{gal}}^{(600)}$ by construction; the apparent $R^2=0.99$ for the gap-vs-empirical comparison reflects this identity, not predictive power."

**Schritt 2:** Tabelle 2 (Z. 612–636), Caption ergänzen: „The ‚stable' subset excludes $\chi_{21}$ \emph{post hoc}, after its $N$-oscillation was identified (§VI). Statistics on the stable subset are therefore conditional on this identification."

**Schritt 3:** Abstract (Z. 96–100) ergänzen: „...drops to $8/10$ owing to two sign failures: $\chi_{21}$ ($N$-oscillator, see §VI) and $\chi_{29}$ (numerical interpolation artifact at near-zero gap)." Das ersetzt die jetzige Formulierung, die nur einen der beiden Fälle benennt.

### E.4 §7.3 ↔ Appendix B Widerspruch glätten [B4 Step 1]

**Schritt 1:** Z. 1098–1115 (§7.3 sign-convention puzzle) umformulieren:

```latex
\subsection{Sign convention: from empirical fix to rigorous derivation}
\label{sec:discussion-sign}

The coefficient $-2$ in \eqref{eq:S-def} was originally fixed
empirically during early versions of this work, against the
Galerkin matrix-assembly convention (function \texttt{build\_W}).
In the present version, the coefficient is rigorously derived
in Appendix~\ref{app:sign} and \cite{AnalyticKernelV2},
Theorem~5.1, by tracking all six relevant sign sources. The
two equivalent conventions—Connes-operator and Weil-$Q$—are
related by an overall sign on the prime part. The atlas uses
the Connes-operator convention throughout, in which
coefficient~$-2$ is the correct one.
```

Dies ersetzt die jetzige Formulierung, die einen historischen Schwebezustand suggeriert.

### E.5 Spectral-Convergence-Statement in §3 + Vier-Achsen-Präzisierung [B.4 Step 2 + B5 Step 1]

**Schritt 1:** In §3 Setup (am Ende von §3.4 Sektor-Zerlegung, vor §3.5) einen Satz ergänzen:

„The Galerkin-truncated operators $W^\pm_{\chi, N}$ are finite-dimensional approximations of the bounded symmetric kernel operator $K^\pm_\chi$ on $L^2[-L,L]$. Standard spectral-convergence theory for compact symmetric operators \cite{BabuskaOsborn} ensures that the discrete eigenvalues converge to the true spectrum as $N\to\infty$, with rates determined by the smoothness of the kernel and the eigenfunction. The empirical $N$-dependence of $\gap_{\mathrm{gal}}^{(N)}(\chi_{21})$ in §VI illustrates that this convergence may exhibit non-monotone transient behavior at moderate $N$ for characters with quasi-degenerate spectra."

(Dafür Babuška–Osborn 1991 oder Trefethen–Bau 1997 in die Bibliographie aufnehmen.)

**Schritt 2:** §4 Einleitung (Z. 707–714) präzisieren: „The raw table of §V admits a structural reading along four taxonomic dimensions: two filter constants (parity and root number, both equal to $+1$ across the entire sample), and two informative axes (archimedean weight ratio and gap signature). We develop each in turn, in order of increasing informativeness."

**Schritt 3:** Abstract (Z. 96–98) entsprechend: „We catalog the ten characters along four axes (two filter constants—parity, root number—and two informative axes—archimedean factor weight, gap signature)..."

---

## Zusammenfassung der Step-2-Bewertung

Das Paper ist **wissenschaftlich solide und methodologisch ehrlich**, hat aber an drei Stellen einen Lücke zwischen formaler Definition und numerischer Implementation oder zwischen Aussage und Beweis (Δ_χ, Prop 6.1, post-hoc selection). Diese drei Lücken sind **nicht-trivial**, aber in 2–3 Stunden Bearbeitung schließbar. Step 1 hat sie unterschätzt — vermutlich, weil sie nur bei sehr aufmerksamer Lektüre der numerischen Anhänge sichtbar werden.

Mit dem oben skizzierten Top-5-Plan erreicht das Paper 8/10 und ist für einen Zenodo-Preprint einreichungsbereit. Zwei zusätzliche Konsolidierungen (B.4 Spectral-Convergence-Beleg im Hauptteil; B.5 spektraler Beleg der Quasi-Degeneracy für χ_21) heben es auf 8.5/10 und machen es Journal-tauglich. Beide sind nachträglich in Revision 2 nachzuholen.

**Empfehlung an die Pipeline:** Top-5 vor Zenodo umsetzen, dann Step 3 (Widerleger) ansetzen. B.4 und B.5 explizit in den Open-Issues-Abschnitt der KONZEPT.md aufnehmen für die spätere Journal-Submission.
