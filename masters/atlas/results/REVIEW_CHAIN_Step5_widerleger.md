# REVIEW CHAIN — Step 5: Widerleger (adversarial)
## Paper: Dirichlet Character Atlas — A Weil-Kernel Cartography of Low-Conductor Sector Gaps (v2, Final Draft)
## Datum: 2026-04-29
## Rolle: Skeptischer Gutachter — keine Lösungsvorschläge, nur Angriffspunkte

---

## Eröffnung

Steps 3 und 4 haben das Paper als `arXiv 9/10` und `Journal 7.5/10` eingestuft. Sie haben die *narrative* Konsistenz geprüft und kosmetische Korrekturen verlangt. Sie haben nicht geprüft, ob die zentralen Identitäten des Papers überhaupt mathematisch tragen. Genau dort setzt dieser Bericht an. Steps 3/4 haben drei Punkte explizit als "ästhetisch", "akzeptabel" oder "Standard" abgelehnt — ich greife genau diese drei Punkte als substantielle Defekte an.

---

## ANGRIFF 1 (kritisch): Die Δ_χ-Definition kollidiert mit der Skript-Implementierung — die "äquivalenten Konventionen" sind nicht äquivalent

**Stelle:** Eq. (5), Z. 467–471; Parität-Sign-Note Z. 472–478; Anhang B Z. 1306–1322.

**Behauptung des Papers:** Der Skript verwendet `np.correlate`, das für φ⁻ "einen paritätsinduzierten Vorzeichenwechsel einführt"; das numerische Resultat ist (φ⁺∗φ⁺) + (φ⁻∗φ⁻) statt der Differenz in (5); "Both conventions yield the same value of S_χ".

**Warum das nicht gilt:** Die Behauptung enthält einen mathematischen Fehler und einen unbewiesenen Sprung.

**(a) Der Paritätsidentitätssatz in Anhang B ist falsch.** Z. 1311–1313 behauptet:
```
R_{φ⁻}(u) := ∫ φ⁻(s)φ⁻(s-u)ds = -(φ⁻∗φ⁻)_math(u)
```
Das ist falsch. `np.correlate(phi, phi, mode='full')` berechnet für reale Arrays exakt R_φ(u) = ∫ φ(s)φ(s-u)ds. Für eine reale (gerade *oder* ungerade) Funktion gilt: R_φ ist *immer eine gerade Funktion* von u und R_φ(0) = ‖φ‖² ≥ 0. Insbesondere ist R_{φ⁻}(0) = ‖φ⁻‖² = +1 (nicht -1) für L²-normiertes φ⁻. Es gibt keinen "paritätsinduzierten Sign-Flip". Die Konvolution (φ⁻∗φ⁻)(u) = ∫ φ⁻(s)φ⁻(u-s)ds ist ein *anderes* Objekt — sie ist gerade für gerade φ und ungerade für ungerade φ, und es gilt (φ⁻∗φ⁻)(0) = -‖φ⁻‖². Die Identität R_{φ⁻} = -(φ⁻∗φ⁻) stimmt also nur an *speziellen* Lags, nicht generell. Anhang B verwechselt Korrelation und Konvolution.

**(b) Die Differenz der beiden "Konventionen" verschwindet nicht generell.** Sei C_χ(t) := (φ⁺∗φ⁺)(t) - (φ⁻∗φ⁻)(t) (Paper-Konvention) und N_χ(t) := R_{φ⁺}(t) - R_{φ⁻}(t) = (φ⁺∗φ⁺)(t) + (φ⁻∗φ⁻)(t) (Skript-Konvention, korrekt verstanden). Dann:
```
N_χ(t) - C_χ(t) = 2(φ⁻∗φ⁻)(t).
```
Bei t=0: N_χ(0) - C_χ(0) = 2(φ⁻∗φ⁻)(0) = -2‖φ⁻‖² = -2. Das ist *kein* Floating-Point-Rauschen. Damit das Paper-Theorem (Eq. 14, gap = S + archdiff) und sein Skript-Pendant denselben Wert von S_χ liefern, muss
```
∑_{p,m} χ(p)^m log(p)/p^{m/2} (φ⁻∗φ⁻)(m log p) = 0
```
gelten. Dieser Satz steht nirgendwo. Anhang B verschiebt ihn an [AnalyticKernelV2] §4 — eine `.md`-Datei, die kein Reviewer einsehen kann.

**(c) Skript-Empirik widerspricht der Behauptung.** Der Skript-Kommentar in `ground_state_difference_analysis.py` Z. 80–87 sagt wörtlich: *"Ein initialer coeff=+2 gab anti-korrelierte Predictions (R=-0.70, 2/10); Flip auf -2 liefert R=+0.70, 8/10."* Das ist keine Vorzeichen-Konvention, sondern **post-hoc Vorzeichen-Anpassung an die Stichprobenkorrelation**. Der Koeffizient -2 wurde gewählt, weil er das gewünschte Vorzeichen lieferte, nicht weil er aus einer unabhängigen Herleitung folgt. §7.3 (Z. 1124–1129) gibt das selbst zu — und die spätere "rigorous derivation" reproduziert genau den vorher empirisch gefitteten Wert, was per Konstruktion keine unabhängige Verifikation ist.

**Pathologie:** Die zentrale numerische Identität des Papers (Eq. 14: gap = S + archdiff) hängt an einer Definitions-Implementierungs-Diskrepanz und einem empirisch gefitteten Vorzeichen. Step 4 hat Option 1 (Δ_χ-Umdefinition zur Summe) als "ästhetischen Aufwand" abgelehnt — der Defekt ist substantiell. Der Reviewer von Experimental Mathematics wird beide Konventionen am gleichen Datenpunkt (z.B. t=0) testen wollen und dabei sofort die Diskrepanz sehen.

---

## ANGRIFF 2 (kritisch): "Empirische Referenz" ist Galerkin-Selbstabbild — die Tab.-2-Statistiken sind Selbstkorrelationen

**Stelle:** Tab. 1 (Z. 700–709), insbesondere die Spalte "gap_emp (N_src)". Tab. 2 (Z. 632–646), Zeile "gap_gal (all 10), N=600: 10/10, R²=0.99".

**Behauptung des Papers:** Das Paper deklariert gap_gal^(600) als "10/10 sign accuracy, R²=0.99" gegenüber gap_emp. §7 (Z. 681–683 in der Tabellen-Caption, Z. 690–691) gibt zu, dass diese Übereinstimmung für sechs von zehn Charakteren *identisch* ist, weil die Referenz selbst aus derselben Galerkin-Diagonalisierung stammt.

**Warum das nicht gilt:** Die Statistik in Tab. 2, Zeilen "gap_gal (all 10): 10/10, R²=0.99" und "gap_gal (stable 9): 9/9, R²=0.99", ist methodologisch leer. Sie misst nicht Übereinstimmung mit einem unabhängigen Referenzwert, sondern Selbstkonsistenz von gap_gal^(600) mit sich selbst. Die Behauptung in §3.3 (Z. 615–618), das Paper liefere "linear correlation R²=0.95 gegen die empirischen Gaps" mit S+archdiff, ist nur dann gehaltvoll, wenn das "Empirische" extern ist — was in 6/10 Fällen *nicht* der Fall ist. Selbst die vier Fälle mit N_src=400 sind nicht extern: Sie sind Galerkin-Werte derselben Methode, lediglich bei kleinerer N.

**Pathologie:** Es existiert in diesem Atlas keine *einzige* Referenz-Gap, die nicht aus derselben Galerkin-Methode generiert wurde. Die Behauptung "predictor agrees with empirical" ist einer Form, in der "empirical" ein zirkulär definierter Begriff ist. Eine echte externe Referenz wäre z.B. (a) hochpräzise Zahlen aus mp-arithmetic mit N ≥ 2000, (b) eine Berechnung über die Spektraldarstellung der Nullstellen via LMFDB, oder (c) ein unabhängiges Galerkin-Verfahren mit anderer Basis. Das Paper liefert nichts davon. Step 3 hat das nicht angegriffen; Step 4 hat das implizit als "Tautologie offen zugegeben" akzeptiert — das genügt nicht: solange die Statistiken mit "gap_emp" als Referenz gefahren werden, sind sie irreführend.

---

## ANGRIFF 3 (kritisch): Tab. 2 "stable 9" ist Selektion auf das Versagensbeispiel — die Statistik ist methodologisch leer

**Stelle:** Tab. 2 Z. 638–639, 643–644 ("stable 9"-Zeilen). Caption Z. 627–629: "'Stable' excludes χ_21 post hoc, after its N-oscillation was identified as a Galerkin-truncation artefact (...). Statistics on the stable subset are therefore conditional on this identification."

**Behauptung des Papers:** Mit einem post-hoc-Caveat sei die "stable 9: 9/9, R²=0.80"-Statistik weiterhin meldenswert.

**Warum das nicht gilt:** Der Caveat macht den Defekt nicht weg. Die Stichprobe ist 10. Das einzige Versagensbeispiel im Vorzeichen ist genau das, das ausgeschlossen wird. Die "9/9"-Statistik ist äquivalent zur Aussage "wenn man den einzigen Sign-Fail entfernt, hat man keinen Sign-Fail mehr". Sie hat null inferentiellen Gehalt: bei beliebiger Wahl des Predictors auf einer Stichprobe von 10 ist nach Ausschluss des Worst-Case-Punkts die Sign-Treffer-Quote trivial maximal.

**Pathologie:** Wenn man χ_29 (zweiter Sign-Fail bei N=600) ebenfalls mit einem post-hoc-Argument (Z. 612–613: "numerical interpolation artifact at a near-zero gap") ausschließt, bleibt eine "8/8"-Statistik. Wenn man konsequent so verfährt — jeden Fail mit einem charakter-spezifischen Caveat versieht — bleibt eine perfekte Statistik aus null Fails durch Selektion. Step 3 hat das in §B.7 (nicht umgesetzt) und Step 4 hat das nicht erneut angegriffen. Die Caption-Klausel reicht nicht. **Die Zeilen "stable 9"** in Tab. 2 sind statistischer Müll und gehören entfernt; sie als legitimes Ergebnis zu führen ist methodisch fragwürdig.

---

## ANGRIFF 4: Theoreme 3.2 und 3.3 sind im Paper unbewiesen — F4 hängt damit in der Luft

**Stelle:** Theorem 3.2 (Z. 339–347), Theorem 3.3 (Z. 353–362), Beweise verweisen auf [AnalyticKernel] (`.md`-Datei). F4 (Z. 211–218): "Weil-kernel architecture is rigorous."

**Behauptung des Papers:** Die strukturellen Theoreme seien "rigorous" und tragen F4.

**Warum das nicht gilt:** Eine Aussage ist nicht "rigorous", weil ein Autor eine `.md`-Datei mit demselben Namen anlegt. Beide Theoreme sind im Paper *Behauptungen ohne Beweis*. Die Verweise auf [AnalyticKernel] §3.2/§4.1 (Z. 339, Z. 353) sind nicht extern verifizierbar — die Datei ist ein lokaler Pfad (`CORE/zoo-mapping/ANALYTIC_PIPELINE.md`), ein Reviewer hat keinen Zugriff. Step 2 hat "Beweisskizzen" gefordert; Step 4 hat sie *vorgeschlagen* (B.1, "5 Minuten"), aber als nicht-blocking eingestuft. Solange sie nicht im Paper stehen, ist F4 nicht aus dem Paper verifizierbar.

**Pathologie 1 (Theorem 3.2):** Das Theorem behauptet eine konkrete Form für K_χ^±(x,y) unter P·κ_χ(·)·P = κ_χ(-·). Aber κ_χ ist eine *Distribution* (Eq. 4: Summe von δ-Funktionen), keine Funktion. Die Identität "Symmetrising with P^±" auf Distributionen erfordert eine Wahl von Test-Funktionen-Raum, in dem alles wohldefiniert ist. Das Paper liefert das nicht.

**Pathologie 2 (Theorem 3.3):** "The character dependence enters only through κ_χ(x+y), since K^arch_χ depends on χ only through the constant log(q/π)" — das ist nicht ganz korrekt. Der archimedische Kernel hängt auch über den Re ψ(1/4 + it/2)-Term von χ ab, *wenn* die Parität wechselt; das wird nur im even-Fall trivial. Die Behauptung im Theorem ist also auf den even-only-Fall des Atlasses spezifisch, was im Theorem-Statement nicht klar gesagt wird.

---

## ANGRIFF 5: χ_21-"N-Oszillator"-Erklärung ist Hypothese, nicht Messung

**Stelle:** §6 (Z. 957–1039), insbesondere die "Quasi-Degeneracy"-Erklärung Z. 996–1004. Tab. 4 (Z. 970–987).

**Behauptung des Papers:** Der Sign-Flip von χ_21 zwischen N=200 und N=600 sei kein Formel-Defekt, sondern ein *quasi-degenerate collapse* der Galerkin-Operatoren bei niedrigem N.

**Warum das nicht gilt:** Tab. 4 zeigt drei Datenpunkte (N=200, 400, 600) und einen Vergleich zu χ_12. Es fehlt jegliche spektrale Diagnose: λ₂-λ₁, Eigenvektor-Überlappungen |⟨φ^(N=200), φ^(N=600)⟩|, Konvergenzdiagramme bei N=500, 700, 800, 1000. Ohne diese Daten ist die "Quasi-Degeneracy"-Erklärung nicht von einer alternativen Hypothese unterscheidbar:

**Alternativhypothese:** Der "Galerkin-Operator" für χ_21 konvergiert *nicht* — das wahre gap_{χ_21}(λ) ist Null oder O(N^{-α}) für irgendein kleines α, und der gemessene Wert -0.0042 bei N=600 ist genauso unstabil wie der Wert +0.2816 bei N=200. Bei N=1000 könnte er erneut springen. Das Paper liefert keinen Beweis, dass die *Konvergenz erreicht wurde*. Step 3 (B.10) hat das angegriffen; Step 4 (B.5) hat es als "Phase 2"-Item akzeptiert. Solange das nicht im Paper steht, ist die Erklärung Spekulation.

**Pathologie:** Die "Bremse" (Rem. 6.4 Z. 1027–1039) — "the right theoretical object is gap_emp(χ_21) = -0.00423 at N=600" — verschiebt die Frage, ohne sie zu beantworten. Wenn χ_21 ein N-Oszillator ist, dann ist auch -0.00423 kein theoretisches Objekt, sondern bloß die nächste Galerkin-Stichprobe. Das Paper liefert keine Argumentation, warum gerade *diese* Stichprobe konvergiert ist.

---

## ANGRIFF 6: Siegel-Walfisz-Kompression als "quantitative image" ist nicht hergeleitet, sondern nachgereichte Erzählung

**Stelle:** Rem. 4.1 (Z. 666–683), insbesondere "This dichotomy is the quantitative image of the Siegel-Walfisz upper bound (...)".

**Behauptung des Papers:** Die zwei N-Verhaltensgruppen (stabil vs. komprimiert) seien das quantitative Bild der Siegel-Walfisz-Schranke.

**Warum das nicht gilt:** Die Siegel-Walfisz-Schranke ist eine *Oberschranke* auf Primsummen für Charaktere mit nicht-Siegel-Nullstellen; sie liefert keine Aussage über N-Konvergenz von Galerkin-Eigenwerten bei λ=const. Das Paper gibt Z. 681–683 selbst zu: "A full theoretical derivation of this compression is open." Das heißt: Die "quantitative image"-Phrase wird ohne Herleitung als Erklärung präsentiert. In wissenschaftlicher Sprache: keine Erklärung, sondern Spekulation mit suggestivem Vokabular.

**Pathologie:** Wenn die Behauptung wahr wäre, müsste der Kompressionsfaktor ~0.35 aus der Siegel-Walfisz-Konstante des jeweiligen Charakters ableitbar sein. Tab. 1 zeigt Faktor ~0.35 für χ_12, χ_17, χ_24, χ_29, χ_60 — ohne theoretische Erklärung, warum gerade dieser Faktor und gerade diese fünf Charaktere. Es ist *Mustererkennung* in einer Stichprobe von 10, präsentiert als strukturelles Resultat.

---

## ANGRIFF 7: Trichotomie prime/mixed/arch-dominated ist Overfitting an die Stichprobe

**Stelle:** §4.3 (Z. 759–815), Tab. 3.

**Behauptung des Papers:** Die zehn Charaktere zerfallen in drei strukturelle Klassen — prime-dominated (r ≤ 25%), mixed (60%–70%), arch-dominated (>100%).

**Warum das nicht gilt:** Die Schwellenwerte sind an die Daten gefittet. Aus Tab. 3:
- r ∈ {1.8%, 7.4%, 13%, 17%, 21%, 22%, 25%}: sieben Werte unter 25%
- r ∈ {63%, 64%}: zwei Werte im "mixed"-Bereich
- r = 356%: ein Wert oben

Es gibt keinen einzigen Datenpunkt im Bereich r ∈ [25%, 60%] oder r ∈ [70%, 100%]. Das ist genau, was Overfitting zeigt: man wählt Schwellen so, dass keine Stichprobe nahe der Grenze liegt. Die "drei Klassen" sind schlicht die drei Cluster, die in der Stichprobe aufgetreten sind — keine strukturelle Vorhersage.

**Pathologie:** Eine echte Klassifikation würde testen, ob ein neuer Charakter (z.B. χ_37, χ_40, χ_41) in eine der drei Klassen *vorhergesagt* fällt. Das Paper liefert keine Vorhersage, sondern Beschreibung der zehn vorhandenen Werte. Mit sechs Sub-Klassen (Z. 893–894) auf zehn Charakteren beträgt das Klassen-zu-Datenpunkte-Verhältnis 0.6 — gemäß üblicher Heuristiken (z.B. mindestens 5 Datenpunkte pro Klasse für stabile Inferenz) ist das stark overfittet.

---

## ANGRIFF 8: "Fortuitous smallness" der N=200-Übereinstimmung ist nicht hergeleitet — alternative Lesart bleibt offen

**Stelle:** §5 (Z. 925–935), insbesondere "the '9/9 sign accuracy at N=200' (...) was due to *fortuitous smallness* of this discrepancy at low N, not to the formula's predictive power."

**Behauptung des Papers:** Die N=200-Übereinstimmung war Glück.

**Warum das nicht gilt:** Das ist die *bevorzugte* Lesart, nicht eine bewiesene. Eine alternative Lesart ist mindestens ebenso plausibel:

**Alternative Lesart:** Bei N=200 war die Formel prädiktiv, weil die Galerkin-Approximation die *richtige* niedrigfrequente Struktur erfasst hat. Bei N=600 tritt eine *neue* Pathologie auf — die im Skript erwähnten 10⁻³-Floating-Point-Differenzen werden vorzeichenrelevant, weil |gap| bei einigen Charakteren in dieselbe Größenordnung fällt. Das wäre kein "fortuitous", sondern ein *fehlerhaft konditioniertes* Berechnungsschema bei hohem N und kleinem Gap.

Das Paper trennt die zwei Lesarten nicht. Es nimmt einfach die "fortuitous"-Interpretation an, ohne sie gegen die Konditionszahl-Hypothese zu testen (z.B. durch Berechnung mit erweiterten Floats oder durch Vergleich mit einem Referenzwert in mp-arithmetic).

**Pathologie:** Wenn die alternative Lesart zutrifft, ist das *gesamte* Atlas-Ergebnis bei N=600 mit Unsicherheit ±10⁻³ behaftet — was viele der "kleinen" Gaps (χ_17, χ_21, χ_24, χ_29, χ_60, alle mit |gap| ≤ 0.015) auf Vorzeichen-Niveau in Frage stellt. Das Paper ignoriert diese Möglichkeit.

---

## ANGRIFF 9: GRH-Annahme für composite D=21,33,60 ist unbelegt

**Stelle:** §2.1 Z. 268–270: "(GRH for the tested conductors, which is classical for all D ≤ 10⁴ in this range via [LMFDB])."

**Behauptung des Papers:** GRH gilt für die getesteten Conductors per LMFDB.

**Warum das nicht gilt:** Die Aussage "GRH ist klassisch für D ≤ 10⁴" ist nicht präzise. LMFDB enthält numerisch *verifizierte Nullstellen auf der kritischen Geraden bis zu einer endlichen Höhe* T — typischerweise T ∈ [10⁴, 10⁷] je nach Konduktor. Es gibt keinen klassischen Satz, der GRH für L(s, χ_D) mit D ≤ 10⁴ unbedingt beweist. Die Behauptung im Paper ist eine ungenaue Verkürzung. Bei composite Conductors D=21 (q=84 als Modul mit Reduktion auf prime Conductor 21), D=33 (q=33), D=60 (q=60) müsste das Paper genau angeben, *bis zu welcher Höhe* die Nullstellen verifiziert wurden, und ob die im Atlas relevanten Höhen (γ ≤ T_kritisch) tatsächlich abgedeckt sind. Bei λ=20000 und L ≈ 9.9 ist die Test-Funktionen-Bandbreite ~2L ≈ 20, was Nullstellen bis ungefähr T ~ 10⁴ erfasst — die LMFDB-Coverage ist zu prüfen.

**Pathologie:** Wenn die Annahme GRH lokal verletzt ist (eine Nullstelle abseits der kritischen Geraden), bricht die Definition der Weil-Quadratform zusammen — die Summenformel Eq. (2) hat dann komplexe Beiträge, und die Eigenwertanalyse verliert ihre Grundlage. Das Paper hat dies nicht in einer eigenen Sektion abgesichert, sondern mit einem Halbsatz erledigt.

---

## ANGRIFF 10: Slope 0.72 als "N=200-spezifische Konstante" widerspricht der Verwendung der N=200-Statistiken im Paper

**Stelle:** §7.2 (Z. 1097–1122), insbesondere "The 0.72 is therefore an N=200-specific Galerkin-truncation constant, not a theoretical invariant of the formula."

**Behauptung des Papers:** Slope 0.72 sei ein N=200-Artefakt.

**Warum das nicht gilt — der Selbstwiderspruch:** Wenn Slope 0.72 ein N=200-Artefakt ist, weil bei N=200 die Galerkin-Approximation noch nicht konvergiert war, dann ist *jede* andere N=200-Statistik im Paper genauso ein Artefakt. Insbesondere:
- Tab. 2, Zeile "S+archdiff (all 10), N=200: 9/10, R²=0.42" — N=200-Artefakt.
- Tab. 2, Zeile "S+archdiff (stable 9), N=200: 9/9, R²=0.80" — N=200-Artefakt.
- Tab. 3 (Z. 768–794) "Archimedean weight ratio at N=200" — sollte demnach bei N=600 neu berechnet werden, denn die Verhältnisse r_χ ändern sich mit N.

Das Paper kann nicht gleichzeitig sagen "N=200-Slope ist ein Artefakt" *und* "N=200-Trichotomie aus Tab. 3 ist strukturell". Step 4 (C.6) hat das als "Slope-0.72 integrieren oder kürzen" bei den nicht-blocking Items abgelegt — das ist zu schwach: der Selbstwiderspruch betrifft das gesamte Atlas-Klassifikationsschema.

**Pathologie:** Eine konsistente Position wäre, *alle* N=200-Statistiken aus dem Paper zu entfernen und nur N=600 zu zeigen. Das würde aber die Behauptung "drei Galerkin-Resolutionen" entkräften, und es würde Tab. 3 (auf N=200 basierend) und damit die Trichotomie auflösen. Das Paper steht zwischen den Stühlen.

---

## ANGRIFF 11: Sign-Konvention-Doppeldeutigkeit — "Sign accuracy" testet Selbstkonsistenz, nicht externes Vorzeichen

**Stelle:** Anhang B (Z. 1297–1304): "Two equivalent conventions (...) give oppositely signed gaps."

**Behauptung des Papers:** Connes- und Weil-Q-Konventionen sind äquivalent, geben aber entgegengesetzte Vorzeichen.

**Warum das nicht gilt:** Wenn das Vorzeichen einer berechneten Größe konventionsabhängig ist, dann ist eine "Sign accuracy"-Statistik nur dann gehaltvoll, wenn (a) eine externe, konventionsfreie Referenz existiert, und (b) der Predictor und die Referenz in derselben Konvention berechnet werden. Im Paper:
- Predictor S+archdiff: in Connes-Konvention (-2 Koeffizient).
- Referenz gap_emp: aus derselben Galerkin-Diagonalisierung, *in derselben Konvention*.

"Sign accuracy" testet also lediglich, ob der Predictor mit sich selbst (in einer gewählten Konvention) konsistent ist — was nach Prop. 6.1 (Tautologie) für die Matrixrealisierung trivial ist. Die einzige nicht-triviale Sign-Aussage wäre Übereinstimmung mit einer *unabhängig* berechneten Größe, z.B. dem Vorzeichen, das aus der spektralen Auswertung der LMFDB-Nullstellen folgt.

**Pathologie:** Wenn man die andere ("Weil-Q") Konvention wählt, dreht sich jedes Vorzeichen um. Die Statistik "8/10 sign accuracy at N=600" wird dann zu "2/10 sign accuracy" — beides gleichermaßen aussagekräftig, weil die Statistik konventionsabhängig ist. Das Paper hat keine Lösung für diesen Status; Anhang B verweist nur auf Konsistenz innerhalb des gewählten Frameworks.

---

## ANGRIFF 12: λ=20000 ist ein einziger Datenpunkt — Asymptotik-Behauptung in (I2) nicht gestützt

**Stelle:** §7.1 (Z. 1056–1075), insbesondere (I2): "The gap gap_χ(λ) is O(1) in λ rather than Ω(√λ) as in the trivial case (...)".

**Behauptung des Papers:** gap_χ(λ) = O(1) — eine Asymptotik in λ.

**Warum das nicht gilt:** Eine Asymptotik in λ kann nicht aus Messungen bei *einem* Wert λ=20000 abgeleitet werden. Das Paper räumt ein (Z. 845–854 §4.4), dass χ_5 und χ_12 in [AsymptoticScan] über drei Größenordnungen bei N=200 untersucht wurden — aber [AsymptoticScan] ist eine `.md`-Datei, und die Resultate werden im Atlas selbst nicht gezeigt. Ohne λ-Skala-Daten für alle zehn Charaktere ist die "O(1) in λ"-Aussage Spekulation.

**Pathologie:** Wenn gap_χ(λ) tatsächlich O(λ^{1/4}) oder O(log λ) wäre, würden die Werte bei λ=20000 und λ=200000 unterschiedliche Vorzeichen oder Größenordnungen haben. Das Paper liefert keinen λ-Scan im Atlas selbst. Step 1 (B.8) und Step 4 (B.6) haben das als "Begründung des Cutoffs" diskutiert, aber niemand hat angegriffen, dass die Asymptotik-*Behauptung* in (I2) auf einem einzelnen λ-Wert ruht.

---

## Zusammenfassende Position des Widerlegers

Die Liste der zwölf Angriffspunkte zerfällt in drei Substanzklassen:

**Klasse A — strukturelle Defekte, die das Paper nicht halten kann:**
1. Δ_χ-Definition vs. np.correlate-Implementierung (mathematischer Fehler in App. B + post-hoc Vorzeichen-Fitting).
2. "Empirische Referenz" ist Galerkin-Selbstabbild — alle Statistiken mit "gap_emp" sind zirkulär.
3. "Stable 9"-Selektion ist methodologisch leer.

Diese drei Punkte sind nicht mit kosmetischen Korrekturen zu beheben. Das Paper braucht entweder (a) einen externen, konventionsfreien Referenzwert, (b) eine eigenständige Herleitung der Δ_χ-Identitäten, und (c) Verzicht auf die "stable 9"-Statistiken.

**Klasse B — substantielle Lücken, die ein Reviewer angreifen wird:**
4. Theoreme 3.2 und 3.3 unbewiesen im Paper.
5. χ_21 ohne spektrale Diagnose.
6. Siegel-Walfisz-Behauptung ohne Herleitung.
7. Trichotomie ist Overfitting.
8. "Fortuitous smallness" ist eine Lesart, keine Beweiskette.

**Klasse C — Selbstwidersprüche und Annahmen ohne Beleg:**
9. GRH-Annahme für composite Conductors unpräzise.
10. Slope-0.72-Argument widerspricht der eigenen Tab.-3-Verwendung.
11. Sign-Statistik konventionsabhängig.
12. O(1)-Asymptotik aus einem λ-Wert.

**Status:** Die Klasse-A-Punkte sind blockierend für *jede* seriöse Submission, einschließlich arXiv-Preprint. Step 4 hat die Fehleinschätzung, dass D.3 (Δ_χ-Umdefinition) "ästhetisch" sei — der Defekt ist *substantiell*. Step 3 hat den Selbstwiderspruch in Angriff 10 nicht gesehen. Beide Steps haben den Skript-Empirik-Befund (post-hoc Sign-Fitting in `ground_state_difference_analysis.py` Z. 80–87) übersehen, der §7.3 entkräftet.

Mit diesen Punkten ist die Readiness-Einstufung "arXiv 9/10" nicht haltbar. Eine ehrliche Bewertung wäre: arXiv-Submission akzeptabel als *Diskussionspaper mit offenen Defekten*, nicht als Beitrag mit prädiktiver Substanz; Journal-Submission ohne Bearbeitung der Klasse-A-Punkte unangemessen.

---

**Relevante Pfade:**
- Paper: `DIRICHLET_CHARACTER_ATLAS_v1_en.tex`
- Skript: `_scripts/ground_state_difference_analysis.py` (Z. 80–87 in `compute_gap_exact`, Z. 54–62 in `autocorr`)
- Step 3: `_results/REVIEW_CHAIN_Step3_konstruktiv2.md`
- Step 4: `_results/REVIEW_CHAIN_Step4_experte2.md`
