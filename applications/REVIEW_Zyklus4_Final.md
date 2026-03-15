# Review Zyklus 4 -- Finaler Pre-Publication Check

Datum: 2026-03-15
Reviewer: Claude (automatisiert)

---

## 1. N4-Fix (FST-I DE: "derived" -> "reinterpreted" in Conclusion)

**STATUS: NICHT UMGESETZT -- P0**

Die EN-Version verwendet korrekt "reinterpreted" (Zeile 547):
> "The Born rule is **reinterpreted** via envariance as the unique equilibrium distribution."

Die DE-Version sagt stattdessen "abgeleitet" (Zeile 547):
> "Die Born-Regel wird aus Envariance als einzige Gleichgewichtsverteilung **abgeleitet**."

Gleiches Problem im Abstract -- EN Zeile 44: "reinterpreted", DE Zeile 44: "abgeleitet wird".

**Fix:** DE Zeile 547: "abgeleitet" -> "uminterpretiert" (oder "neuinterpretiert"). DE Zeile 44: "abgeleitet wird" -> "uminterpretiert wird".

---

## 2. N5-Fix (FST-I DE: Intro teleologische Formulierung)

**STATUS: NICHT UMGESETZT -- P1**

EN Zeile 60 (korrekt, abgemildert):
> "the universe is modeled as a complex, non-equilibrium thermodynamic system **whose large-scale dynamics are consistent with high entropy production rates** over cosmological epochs"

DE Zeile 60 (noch teleologisch):
> "das Universum als ein komplexes, nicht-gleichgewichtiges thermodynamisches System modelliert, das **kontinuierlich danach strebt, seine Entropieproduktionsrate ... zu maximieren**"

**Fix:** DE Zeile 60 an EN angleichen: "dessen grossskalige Dynamik mit hohen Entropieproduktionsraten ueber kosmologische Epochen konsistent ist".

---

## 3. Verbleibende EN/DE-Divergenzen

### FST-I Conclusion, letzter Absatz (Zeile 555)

- **EN** (abgemildert): "physical laws at different scales **may share** a common optimization logic ... **a hypothesis that requires formal validation** through explicit renormalization-group mappings"
- **DE** (assertiv): "Das Universum **erfindet keine** neuen Regeln ... es **renormiert** die gleiche Optimierungslogik"

**Schweregrad: P1** -- DE macht eine Aussage als Tatsache, EN formuliert es als Hypothese.

### FST-III Abstract (Zeile 45-47)

- **EN:** Ein Satz: "Cosmological extensions are deferred to a companion paper and are not part of the biological theory presented here."
- **DE:** Laengere Formulierung mit Verweis auf `\ref{sec:cosmology}` PLUS ein zusaetzlicher "Hinweis zum Gueltigkeitsbereich"-Paragraph (Zeile 47), der in EN komplett fehlt.

**Schweregrad: P2** -- Die DE-Version ist vorsichtiger (besser), aber die Versionen sollten inhaltlich gleich sein. Empfehlung: Den "Hinweis"-Paragraph auch in EN einfuegen.

### FST-III Conclusion (Zeile 635 EN vs. 602 DE)

- **EN:** "The cosmological extension of this framework is developed separately in \citep{FST-Unified2026}." (assertiv)
- **DE:** "Die spekulative Erweiterung auf kosmologische Skalen (Abschnitt~\ref{sec:cosmology}) verbindet sich mit dem umfassenderen FST-Forschungsprogramm, erfordert aber dedizierte formale Entwicklung" (vorsichtiger)
- **EN letzter Satz:** "Whether the same logic extends to cosmological scales is explored in the companion FST framework paper" (verweist auf companion paper)
- **DE letzter Satz:** "Ob die gleiche Logik bis zu galaktischen Skalen reicht, bleibt eine offene Frage." (kein Verweis)

**Schweregrad: P2** -- Inkonsistenter Tonfall. DE ist vorsichtiger (besser), EN sollte angeglichen werden.

---

## 4. Broken refs / fehlende bibitems

- **FST-I EN/DE:** Keine Probleme. Alle `\citep`/`\cite` haben passende `\bibitem`.
- **FST-II EN/DE:** Keine Probleme.
- **FST-III EN/DE:** 3 **unbenutzte bibitems** in beiden Versionen:
  - `DarkMatterMOND`
  - `RARGalaxies`
  - `SchoolBiologicalPhysics`

  Kein Kompilierfehler, aber unnoetig. Empfehlung: Entfernen oder einbinden.

- **Broken `\ref{}`:** Keine gefunden (alle Labels existieren).

---

## 5. Publishing-Readiness Score (EN-Versionen)

| Paper | Score | Kommentar |
|-------|-------|-----------|
| FST-I EN | **8/10** | Solide. Die abgemilderten Formulierungen (N4, N5) sind in EN korrekt. Die Conclusion ist vorsichtig genug. Abzug: 98.8%-Claim ist ein single data point; die "first semi-analytic calculation" sollte nicht ueberbewertet werden. |
| FST-II EN | **8/10** | Kompakt, gut strukturiert. Abstracts und Conclusions EN/DE konsistent. Keine offenen Fixes. Abzug: "deterministic optimization process" in Intro (Zeile 59) koennte als ueberstarke Behauptung gelesen werden. |
| FST-III EN | **7/10** | Abzug: (a) 3 unbenutzte bibitems, (b) Cosmology-Section als Teil des Papers aber gleichzeitig "deferred to companion paper" -- das ist inkonsistent, (c) der Abstract sagt "not part of the biological theory" aber das Paper enthaelt offenbar einen `sec:cosmology`. |

---

## Zusammenfassung der offenen P0/P1 Items

| # | Paper | Prio | Problem | Fix |
|---|-------|------|---------|-----|
| 1 | FST-I DE | P0 | Born-Regel "abgeleitet" statt "uminterpretiert" (Abstract Z.44 + Conclusion Z.547) | Wort ersetzen |
| 2 | FST-I DE | P1 | Intro Z.60 teleologisch ("strebt danach zu maximieren") | An EN angleichen |
| 3 | FST-I DE | P1 | Conclusion Z.555 assertiv statt hypothetisch | An EN angleichen |
| 4 | FST-III EN | P2 | Fehlender "Scope Note" Paragraph im Abstract | Aus DE uebernehmen |
| 5 | FST-III EN/DE | P2 | Conclusion-Tonfall inkonsistent | Angleichen |
| 6 | FST-III EN/DE | P2 | 3 unbenutzte bibitems | Entfernen |
