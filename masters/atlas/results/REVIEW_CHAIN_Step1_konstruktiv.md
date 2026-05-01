# REVIEW CHAIN — Step 1: Konstruktiver Reviewer (Creative Referee)
## Paper: Dirichlet Character Atlas — A Weil-Kernel Cartography of Low-Conductor Sector Gaps (v1/v2)
## Datum: 2026-04-29

---

## Abschnitt A: Stärken der Arbeit

**A1 — Klare diagnostische Postur.**
Das Paper nimmt eine seltene und mutige Position ein: es ist ein ehrliches Negativ-Resultat, das explizit als solches deklariert wird. Die Formulierung "our contribution is deliberately diagnostic" ist konsequent durchgehalten — die Befunde F1–F4 widersprechen sich nicht, sie bauen aufeinander auf. Diese Haltung ist wissenschaftlich vorbildlich und für die GRH-Community wertvoller als eine vorschnelle positive Behauptung.

**A2 — Drei-Komponenten-Dekomposition als Diagnoserahmen.**
Die Formel RH_X = GBC_X + C^(2)_X + Hurwitz_X wird nicht nur als Etikett verwendet, sondern als operativer Rahmen: jedes Ergebnis wird darauf rückgebunden, welche Komponente es betrifft. Das gibt dem Atlas strukturelle Kohärenz über die reine Datenpräsentation hinaus.

**A3 — Theoreme 3.2 und 3.3 sind sauber formuliert.**
Die Sektor-Kern-Zerlegung (Thm. 3.2) und die Anti-Diagonal-Differenz (Thm. 3.3) sind klar, mit vollständigen Definitionen und sauberer Zuordnung zu Zitaten. Sie bilden eine solide algebraische Basis, auf der alle numerischen Argumente aufbauen.

**A4 — Proposition 6.1 (Tautologie) ist der methodologische Kernbeitrag.**
Die Tautologie-Proposition ist präzise und elegant. Sie benennt die Grenze numerisch-variativer Methoden (Galerkin, aber verallgemeinerbar auf alle Ritz-Galerkin-Diskretisierungen) in einer Form, die über das Dirichlet-Setting hinausweist. Dies ist der stärkste einzelne Beitrag des Papers.

**A5 — χ_21-Analyse ist vorbildlich offen.**
Die Nicht-Konvergenz von χ_21 zwischen N=200 und N=400 (stationär, dann Vorzeichensprung um Faktor 66 bei N=600) wird nicht unter den Tisch gekehrt, sondern als separater Abschnitt ausgearbeitet. Der Hinweis auf Quasi-Degeneracy-Collapse als Ursache ist physikalisch plausibel und metodologisch wichtig.

**A6 — Vier-Achsen-Klassifikation liefert echten taxonomischen Inhalt.**
Sechs Subklassen in zehn Charakteren (strong-neg, strong-pos, moderate-neg, weak-pos, arch-dominated, composite-oscillator) sind eine genuine Gliederung, nicht eine Nachrationalisierung der Rohdaten. Die Trichotomie prime-dominated / mixed / arch-dominated (Tabelle 3) ist für zukünftige Selektionsentscheidungen direkt nützlich.

**A7 — Siegel-Walfisz-Kompressions-Beobachtung.**
Die empirische Dichotomie zwischen N-stabilen und komprimierbaren Charakteren (Remark 4.1) ist ein eigenständiger Befund, der über die reine Falsifikation der Leading-Order-Formel hinausgeht.

---

## Abschnitt B: Konkrete Verbesserungsvorschläge

**(B1) Header-Kommentar "SKELETON" entfernen (minor, aber dringend)**
Zeile 7–8 lautet: "Status: SKELETON — content sections marked as TODO / placeholders". Das Paper ist vollständig ausformuliert. Dieser Kommentar macht bei einem Einreichungs-Preprint keinen guten Eindruck. Löschen oder durch "Status: Final Draft" ersetzen.

**(B2) Versions-Mismatch auflösen (minor)**
Dateiname ist `v1_en.tex`, die `\date`-Zeile sagt "Draft v2 — FST-Mathematics supplement (revised after N=600 server analysis)". Dies ist ein interner Widerspruch, der in jedem Reviewer-Bericht kommentiert wird. Entscheidung treffen und konsistent durchziehen: entweder Datei in v2 umbenennen oder `\date` auf v1 korrigieren.

**(B3) Abstract zu F2 ist unvollständig — beide Sign-Flips benennen (inhaltlich)**
Der Abstract schreibt: "drops to 8/10 because numerical interpolation errors ... flip the sign of the near-zero-gap character χ_29". Tabelle 1 zeigt jedoch zwei (F)-Marker: χ_21 und χ_29. Der Abstract erklärt einen der zwei Fehlschläge, nicht beide. Ein externer Reviewer wird das sofort bemerken. Vorschlag: explizit sagen "two sign failures: χ_21 (N-oscillator, discussed in §6) and χ_29 (interpolation artifact at near-zero gap)".

**(B4) Widerspruch zwischen §7.3 und Appendix B auflösen (inhaltlich, wichtig)**
§7.3 schreibt: "The coefficient −2 ... was fixed empirically during Session 7 of the development [...] A clean theoretical re-derivation ... is listed as open in Appendix B." Appendix B hingegen erklärt: "The coefficient −2 in eq. (5) is rigorously derived from the Weil explicit formula in [AnalyticKernelV2], Theorem 5.1, by tracking all six relevant sign sources." Das ist ein Widerspruch. Appendix B reflektiert den neueren Stand und löst das Problem auf. §7.3 muss nachgezogen werden: die Formulierung "open" ist falsch, der Widerspruch muss gestrichen und durch "resolved in Appendix B and [AnalyticKernelV2]" ersetzt werden.

**(B5) Root-number-Achse (§4.2) nicht als vollwertige Achse verkaufen (minor)**
Alle zehn Charaktere haben W=+1. §4.2 behandelt deshalb nicht eine Achse des Atlas, sondern eine Konstante. Das ist im Text ehrlich gesagt, aber die Überschrift "Root-number classification" und die Einordnung als eine der "vier Achsen" ist irreführend. Vorschlag: entweder umbenennen zu "Root-number layer (degenerate in this sample)" oder die Vier-Achsen-Aussage auf "drei informative Achsen + eine prospektive Dimension" präzisieren.

**(B6) R²=0.99 für gap_gal bei N=600 als inflationären Wert kennzeichnen (inhaltlich)**
Tabelle 1 zeigt: 6 von 10 Einträgen in der Spalte gap_emp (N_src=600) sind per Konstruktion identisch mit gap_gal^(600) — das ist die Tautologie von §6. Der R²=0.99-Wert für "gap_gal vs. gap_emp" bei N=600 ist deshalb nicht informativ; er misst die Identität, nicht die Vorhersagekraft. Entweder in Tabelle 2 eine Fußnote hinzufügen "includes 6/10 tautological identities; informative comparison restricted to 4 characters" — oder eine separate Zeile für die nicht-tautologischen 4 Charaktere ausweisen. Ohne diese Klarstellung könnte ein Referee fälschlicherweise den R²=0.99 als Erfolg der Formel lesen.

**(B7) Selektionskriterium der zehn Charaktere explizit machen**
Das Paper nennt D ∈ {5,8,12,13,17,21,24,29,33,60}, sagt aber nirgends explizit warum gerade diese zehn. Vermutlich: kleinste Konduktoren mit χ(-1)=+1 und q ≤ 60. Eine Ein-Satz-Begründung in §3.1 oder §1.2 würde Reviewer-Rückfragen verhindern.

**(B8) Wahl λ=20000 begründen**
§3.1 nennt λ=20000 ohne Begründung. Ein kurzer Hinweis — z.B. "chosen as the largest λ at which N=600 Galerkin diagonalization remains computationally feasible on the available hardware within acceptable wall-clock time" — verhindert die Standard-Reviewer-Frage "why not larger λ?".

**(B9) Interne Supplementary-Referenzen als unpubliziert labeln**
Die Bibitems AnalyticKernel, AnalyticKernelV2, AnalyticGroundstate, Woche4Validation, AsymptoticScan verweisen auf .md-Dateien in CORE/zoo-mapping/. Für einen Zenodo-Preprint im öffentlichen Companion-Verbund ist das unbefredigend: ein externer Leser kann diese Quellen nicht verifizieren. Kurze Lösung: Bibitems mit "Internal supplementary note, unpublished, on file with the author" kennzeichnen — oder die Notes auf Zenodo als Supplementary Materials hochladen und mit DOI versehen.

---

## Abschnitt C: Kreative / optionale Ideen

**C1 — Tautologie-Proposition als allgemeines methodologisches Resultat hervorheben**
Proposition 6.1 gilt nicht nur für Galerkin-Diskretisierungen des Weil-Kerns, sondern für jede Ritz-Galerkin-Diskretisierung einer symmetrischen bilinearen Form. Die Aussage "Eigenvektoren eines trunkierten Operators können keine Vorhersagekraft über denselben Operator liefern" ist ein generisches Warnsignal für operator-numerische Methoden in der analytischen Zahlentheorie (und darüber hinaus: FEM, PDE-Spektralmethoden, Variationsprinzipien in der Quantenmechanik). Ein Remark, der explizit auf diese Verallgemeinerung hinweist, würde den Paper-Impact deutlich erhöhen und Zitierbarkeit jenseits der engeren GRH-Community schaffen.

**C2 — Siegel-Walfisz-Kompressionsfaktor 0.35 als Hypothese formulieren**
Die empirische Beobachtung, dass prime-dominated Charaktere einen Kompressionsfaktor ~0.35 zeigen (Remark 4.1), ist bisher nur deskriptiv. Eine explizite Hypothese könnte lauten: "Der Kompressionsfaktor c_χ(N) = gap_gal^(N)/gap_gal^(200) für prime-dominated Charaktere konvergiert gegen 1 als N→∞ mit einer Rate, die von der Siegel-Walfisz-Konstante von χ abhängt." Das ist falsifizierbar und würde den Atlas-Befund quantitativ präzisieren.

**C3 — Einen N=1000-Checkpunkt für χ_21 als Stabilitätsargument einplanen**
Der χ_21-Collapse zwischen N=400 und N=600 wird als Quasi-Degeneracy-Collapse erklärt, bleibt aber ohne Nachweis, dass N=600 der tatsächlich konvergente Wert ist. Ein einzelner N=1000-Lauf für χ_21 (bei ~30 Minuten CCX13-Zeit) würde die Behauptung der Konvergenz bei N=600 enorm stärken — ohne den Paper-Scope wesentlich zu erweitern. Dies könnte als Fußnote oder Appendix-Absatz aufgenommen werden.

**C4 — Route D (χ-twisted CCM) als konkreten Ausblickpunkt schärfen**
§8.1 nennt Route D (χ-twisted CCM-Transfer) als ersten Post-Zookeeper-Task. Hier würde eine konkrete Aussage über den erwarteten Formaldurchmesser helfen: Was ist der kleinste Defekt-Norm-Test, der Route D falsifizieren würde? Ein solches Falsifikationskriterium macht den Ausblick für externe Leser angreifbar und für interne Planung verbindlich.

**C5 — Bezug zur Modular-Form-Theorie andeuten (interdisziplinär)**
Die Halbgruppen-Struktur der Kronecker-Charaktere χ_D mit Konduktoren D = 5,8,12,... ist eng mit dem Theorie-Gebäude der Modularformen über Q(√D) verbunden. Ein kurzer Hinweis, dass die Archimedean-Layer-Trichotomie (prime-dominated / mixed / arch-dominated) möglicherweise in der CM-Arithmetik der entsprechenden quadratischen Körper verankert ist, würde den Atlas für ein breiteres algebraisches Publikum öffnen. Dies ist optional und spekulativ — aber für eine "Cartography"-Arbeit genau das richtige Terrain.

---

## Abschnitt D: Review-Antizipation

**D1 — "Was ist neu?" wird gefragt werden.**
Ein Referee, der den Zookeeper-Paper und den Math-Master nicht kennt, wird fragen: "Wherein does this paper's contribution exceed showing that a certain numerical method does not work?" Antwort: Die Vier-Achsen-Klassifikation, die Tautologie-Proposition (allgemein formulierbar) und das χ_21-Konvergenz-Phänomen sind positive Beiträge, keine Null-Resultate. Diese Unterscheidung muss in §1.2 ("What this paper does") noch klarer herausgestellt werden.

**D2 — "Warum nur zehn Charaktere?" wird Standard-Kommentar sein.**
Antwort liegt in §8.4 (Conduct extension D≤100), ist aber reaktiv formuliert. Besser: proaktiv in §1.2 schreiben: "We restrict to ten real characters as a proof-of-concept sample; extension to D≤100 is feasible (§8.4) but would not change the structural findings."

**D3 — "Sind die Tautologie-Folgerungen trivial?"**
Ein strenger algebraischer Reviewer könnte einwenden, dass Proposition 6.1 eine Identität ist, keine Entdeckung. Die Antwort ist: Sie ist nicht trivial in dem Sinne, dass die früheren (N=200)-Resultate fälschlicherweise als Evidenz für Vorhersagekraft interpretiert wurden. Der Wert liegt im Aufzeigen, wo und warum frühere Evidenz fehlgeleitet war — das ist genuine wissenschaftliche Arbeit, auch wenn die Proposition formal simpel ist.

**D4 — Sign-Convention in §7.3 vs. Appendix B wird bemerkt werden.**
Dieser Widerspruch (§7.3: "open"; Appendix B: "rigorously resolved") wird von einem sorgfältigen Referee kommentiert werden. Priorität: vor Einreichung auflösen (vgl. B4 oben).

**D5 — Konvergenzfrage bei N→∞ allgemein.**
Ein theoretisch orientierter Referee wird fragen: "Ist bekannt, dass die Galerkin-Eigenvektoren überhaupt in einem starken Sinne konvergieren?" Die Antwort ist: Für glatte symmetrische Integraloperatoren auf L^2[-L,L] ja, unter Standard-Spectral-Convergence-Theorie. Ein Satz in §2.3 oder §3 ("The Galerkin eigenvalues converge to the true spectrum by standard spectral theory [ref]") würde diese Frage präventiv beantworten.

**D6 — Standalone-Wert ohne Kenntnis des Companion-Ökosystems.**
Das Paper setzt Vertrautheit mit Math-Master, Zookeeper und RH_II voraus. Für Leser aus der GRH-Community (Weil-Kern, Explicit Formula, Connes-Apparat) ohne Kenntnis des internen FST-Programms könnte §1.1 eine 2-Satz-Zusammenfassung des Three-Component-Frameworks enthalten, die nicht auf externe Preprints verweist. Die jetzige Formulierung lässt sich erst nach Lesen von [MathMaster] §2.4 vollständig nachvollziehen.

---

## Readiness-Score: 6/10

**Begründung:**
Der wissenschaftliche Inhalt ist solide und die diagnostische Postur ist überzeugend durchgehalten. Die Struktur und die wesentlichen Theoreme sind vorhanden. Die Punktabzüge gehen auf folgende behebbare Defizite:

- Widerspruch §7.3 ↔ Appendix B (inhaltlich, muss vor Einreichung aufgelöst werden) → −1 Punkt
- Abstract-Unvollständigkeit bei F2 (zwei Sign-Flips, nicht einer erklärt) → −0.5 Punkte
- Inflationärer R²=0.99-Wert ohne Tautologie-Kennzeichnung → −0.5 Punkte
- Header-SKELETON + Versions-Mismatch (cosmetic, aber professionell störend) → −0.5 Punkte
- Selektionskriterium + λ-Begründung fehlen → −0.5 Punkte
- Interne Supplementary-Quellen nicht öffentlich zugänglich → −0.5 Punkte (je nach Ziel-Journal kritisch)

Mit Behebung der Punkte B1–B6 würde das Paper bei **7.5–8/10** liegen und einreichungsbereit für einen Zenodo-Preprint mit Companion-DOI sein.
