# REVIEW Zyklus 2: Formatierung, Fachsprache, Literatur, DE-EN-Gleichwertigkeit

**Datum:** 2026-03-15
**Reviewer:** Claude (akademisches Lektorat)
**Dateien:** FST-I/II/III jeweils `_v2_en.tex` und `_v2_de.tex`

---

## 1. Fachsprache (pro Paper)

### FST-I (Thermodynamic Game Theory)

| Zeile (EN) | Problem | Vorschlag |
|---|---|---|
| EN Z.60 | "mathematically rigorous explanation" -- uebertrieben fuer ein Framework-Paper, das Analogien formuliert | "mathematically grounded interpretation" |
| EN Z.555 | **Duplikation im Conclusion:** Zwei aufeinanderfolgende Saetze beginnen mit "If successful" und widersprechen sich teilweise ("not one of disconnected..." vs. "may share a common..."). Der erste Satz ist apodiktisch, der zweite hedged. | Einen der beiden Saetze streichen; der zweite (hedged) ist konsistenter mit dem Ton des Papers |
| DE Z.60 | "mathematisch rigorose Erklaerung" -- gleiches Problem wie EN | "mathematisch fundierte Interpretation" |
| DE Z.549 | **Kein Duplikat:** DE Conclusion hat nur eine Version ("eines einzigen Spiels, das unter Grobkoernung wiederholt gespielt wird"). Aber dieser Text ist apodiktischer als die EN-Fassung. | An EN-Formulierung (hedged) angleichen |
| EN Z.76 | "structural analogy" -- korrekte, abgeschwaechte Formulierung | OK |
| DE Z.76 | "strukturell notwendige Instanz" -- deutlich staerker als EN ("structural analogy") | An EN angleichen: "strukturelle Analogie" |
| DE Z.424 | "rigorosen mathematischen Praezedenzfall" -- EN Z.428 sagt "structural analogy motivated by..." | An EN angleichen: "strukturelle Analogie, motiviert durch..." |
| DE Z.532 | "rigoros in der eichtheoretischen Reformulierung...etabliert" -- EN Z.536 sagt "developed as structural analogy" | An EN angleichen |

### FST-II (Chemical Game Theory)

| Zeile | Problem | Vorschlag |
|---|---|---|
| DE Z.490 | "Maximierung dissipierter Energie" -- EN Z.504 sagt "stability-selected dissipation". **Dies ist eine inhaltlich kritische Abweichung:** die EN wurde bewusst von "maximization" auf "stability-selected" geaendert. | Korrigieren zu: "stabilitaetsselektierte Dissipation" |
| DE Z.464 | "Systeme sind dynamisch stabil, weil sie die Entropieproduktion maximieren" -- EN Z.466 sagt "because they satisfy the resolvent stability criterion under dissipative constraints, which correlates with replication efficiency" | Die DE-Version reflektiert die alte Maximierungs-Logik, nicht die neue Stabilitaets-Perspektive. Angleichen an EN. |
| DE Z.479 | "rigoros in der eichtheoretischen Reformulierung...etabliert" -- EN Z.491 sagt "developed as structural analogy" | An EN angleichen |

### FST-III (Biological Game Theory)

| Zeile | Problem | Vorschlag |
|---|---|---|
| DE Z.577 | Punkt 4 in "Was ist wirklich neu?": "Kosmologische Verbindungen" -- EN Z.587 sagt "MEPP as stability filter" | Inhaltlich andere Aussage. EN hat den Punkt aktualisiert. DE muss nachgezogen werden. |
| DE Z.589 | "rigoros in der eichtheoretischen Reformulierung...etabliert" -- EN Z.599 sagt "motivated by the structural parallels in..." | An EN angleichen |
| DE Z.502 | "vier Phasen" -- EN Z.516 sagt "three phases" (Phase 4 wurde in EN entfernt bzw. in die Conclusion verlagert) | Angleichen an EN: nur 3 Phasen, Phase 4 streichen oder als kurzen Verweis belassen |

---

## 2. Formatierungsprobleme (pro Paper)

### FST-I

| Zeile | Problem | Schwere |
|---|---|---|
| DE Z.71, 272, 379, 408, 416, 417, 424, 480, 511, 524, 532 | **Nicht-escaped Umlaute** (UTF-8 statt LaTeX-Escaping): "lebensfaehig" statt `lebensf\"ahig`, "Parametersaetze" statt `Parameters\"atze` etc. Insgesamt 13 Stellen. Bei `inputenc utf8` kompiliert es, aber es ist inkonsistent mit dem Rest der Datei, die `\"a` etc. verwendet. | P1 |
| DE Z.33-35 | **Author-Block:** Verwendet `\texttt{lukisch@users.noreply.github.com}` statt ORCID. EN nutzt korrekt ORCID. | P0 |
| EN Z.78 | `\paragraph{Terminology.}` -- Vorhanden in EN, **fehlt in DE**. | P1 |
| EN Z.224 | Martyushev (2010) Absatz mit MEPP-Einschraenkungen -- **fehlt in DE**. Auch `\bibitem{Martyushev2010}` fehlt in DE-Bibliographie. | P0 |
| EN Z.536-538 | Resolvent-Perspektive mit numerischen Ergebnissen (alpha-scan, 98.8%, Vergleich mit Landscape/Anthropic) -- **in DE (Z.532) stark gekuerzt**, kein numerischer Vergleich. | P1 |
| EN Z.555 | **Duplikat-Satz** in der Conclusion: "If successful, the picture that emerges..." gefolgt von "If successful, this picture suggests..." -- ein Satz muss gestrichen werden. | P1 |

### FST-II

| Zeile | Problem | Schwere |
|---|---|---|
| DE Z.47, 61, 74, 163, 234, 236, 338, 342, 356, 376, 378, 406, 479, 490 | **Nicht-escaped Umlaute** (UTF-8): "Homochiralitaet", "verdraengt", "Frequenzabhaengige" etc. Insgesamt 16 Stellen. | P1 |
| DE Z.36-38 | **Author-Block:** `\texttt{lukisch@users.noreply.github.com}` statt ORCID. | P0 |
| EN Z.468-476 | **Gesamte Subsection "The MEPP Controversy and Scope Limitations" fehlt in DE.** DE springt direkt von "Pross' DKS" (Z.462) zu "Was ist wirklich neu?" (Z.466). | P0 |
| DE Z.234, 236 | Tabellen-Inhalt mit `&` -- die `&`-Zeichen in der Tabelle koennten bei falschen column-specifiers Probleme machen. Pruefung ob Kompilierung korrekt laeuft. | P2 |

### FST-III

| Zeile | Problem | Schwere |
|---|---|---|
| DE Z.118, 293, 589 | **Nicht-escaped Umlaute** (UTF-8): 4 Stellen. | P1 |
| DE Z.34-36 | **Author-Block:** `\texttt{lukisch@users.noreply.github.com}` statt ORCID. | P0 |
| EN Z.271 | Proposition hat Label "(Conjecture)" -- **DE Z.268 hat kein "(Conjecture)"-Label**. | P0 |
| EN Z.305 | **Status-Paragraph fehlt in DE.** EN hat: "This proposition is stated as a conjecture supported by qualitative arguments. A formal proof would require..." | P0 |
| EN Z.603-612 | **Gesamte Subsection "MEPP: Scope and Limitations" fehlt in DE.** | P0 |
| EN Z.612 | **Circularity-Caveat-Paragraph fehlt in DE.** | P0 |
| DE Z.472-495 | **Cosmological Extensions:** DE hat ausfuehrlichen Hunt-These-Abschnitt (3 Subsections, ~24 Zeilen). EN Z.505-509 hat nur einen kurzen Verweis auf das Framework-Paper. **DE ist hier ausfuehrlicher als EN -- intentional?** | P1 |
| DE Z.516-518 | **Phase 4 "Kosmische Integration (Spekulativ)"** existiert nur in DE, nicht in EN. | P1 |

---

## 3. Literatur/Zitierfehler

### Pro Paper

#### FST-I

| Problem | Details | Schwere |
|---|---|---|
| **Fehlender bibitem in DE:** `Martyushev2010` | EN Z.631 definiert `\bibitem[Martyushev(2010)]{Martyushev2010}`, referenziert in EN Z.224. In DE existiert weder bibitem noch Zitation. Da der zugehoerige Absatz fehlt, kein broken ref -- aber inhaltliche Luecke. | P0 |
| **Bibitem-Key Alias:** `DynamicalOrigin2005` | Zitiert als "Crooks(2007)" in bibitem -- Key suggeriert 2005, bibitem sagt 2007. Inkonsistenz. | P2 |
| **Bibitem-Key Alias:** `Lindblad2015` | Zitiert als "Lindblad(1976)" -- Key suggeriert 2015. | P2 |
| **Bibitem-Key Alias:** `MEPP2023` | Zitiert als "Martyushev & Seleznev(2006)" -- Key suggeriert 2023. | P2 |
| **Bibitem-Key Alias:** `MFG2025` | Zitiert als "Lasry & Lions(2007)" -- Key suggeriert 2025. | P2 |
| **Bibitem-Key Alias:** `QuantumFoundations2023` | Zitiert als "Kotecha(2024)" -- Key suggeriert 2023. | P2 |
| **Bibitem-Key Alias:** `QuantumCognition2015` | Zitiert als "Busemeyer & Bruza(2012)" -- Key suggeriert 2015. | P2 |
| **Zukunftsdatum:** `Sandhu et al.(2026)` (FST-III) | Zitiert als 2026 -- akzeptabel wenn Preprint, aber sollte als "forthcoming" oder "preprint" markiert sein. | P2 |
| **Zukunftsdatum:** `Safron et al.(2025)` (FST-III) | 2025 -- pruefbar ob publiziert. | P2 |

#### FST-II

| Problem | Details | Schwere |
|---|---|---|
| **Bibitem-Key Alias:** `SchmitzReview2011` | Zitiert als "Plasson et al.(2007)" -- Key suggeriert Schmitz 2011. | P2 |
| **Bibitem-Key Alias:** `StrategicPersistence` | Zitiert als "Pross(2012b)" -- generischer Key. OK aber inkonsistent mit Stil. | P2 |

#### FST-III

| Problem | Details | Schwere |
|---|---|---|
| Keine fehlenden bibitems | DE und EN haben identische Bibliographien. | -- |

### Cross-Paper-Probleme

| Problem | Details | Schwere |
|---|---|---|
| **Alle 3 DE-Versionen:** Email statt ORCID | `lukisch@users.noreply.github.com` -- dies ist eine GitHub-noreply-Adresse, keine Kontaktadresse. EN nutzt korrekt ORCID. | P0 |
| **Self-Citations:** `Geiger(2026b/c/d/e)` | Die Zuordnung der Buchstaben (b,c,d,e) muss konsistent ueber alle Papers sein. Aktuell: FSTI=2026a (in II,III), FSTII=2026b (in I,III), FST-RH3=2026c, FST-Unified=2026d, FST-RH2=2026e. Das scheint konsistent zu sein. | -- |
| **Inkonsistente Resolvent-Sprache** | EN sagt durchgehend "structural analogy" und "developed as structural analogy"; alle 3 DE-Versionen sagen "rigoros...etabliert". **Dies ist ein systematisches Problem ueber alle 3 DE-Papers.** | P0 |

---

## 4. DE-EN Gleichwertigkeit (pro Paper)

### FST-I: DE liegt hinter EN zurueck

| Was fehlt in DE | EN-Zeile | Schwere |
|---|---|---|
| **Terminology-Paragraph** mit Definition von "resolvent-stable" und Einordnung als "structural analogy" | EN Z.78 | P0 |
| **Martyushev (2010) MEPP-Einschraenkungen** im MEPP-Abschnitt (bilineare Entropieproduktion, LTE-Bedingung) | EN Z.224 (letzter Satz) | P0 |
| **Detaillierte Resolvent-Perspektive** mit numerischen alpha-scan-Ergebnissen, Vergleich mit Landscape/Anthropic | EN Z.536-538 (vs. DE Z.532 stark gekuerzt) | P1 |
| **Bibitem `Martyushev2010`** in Bibliographie | EN Z.631 | P0 |

| Was in DE anders formuliert ist | DE-Zeile | EN-Zeile | Bewertung |
|---|---|---|---|
| "strukturell notwendige Instanz" vs. "structural analogy" | DE 76 | EN 76 | DE zu stark |
| "rigorosen mathematischen Praezedenzfall" vs. "structural analogy motivated by" | DE 424 | EN 428 | DE zu stark |
| "rigoros...etabliert" vs. "developed as structural analogy" | DE 532 | EN 536 | DE zu stark |
| Conclusion: ein Satz vs. Duplikat-Saetze in EN | DE 549 | EN 555 | EN hat Bug (Duplikat) |

### FST-II: DE liegt hinter EN zurueck

| Was fehlt in DE | EN-Zeile | Schwere |
|---|---|---|
| **Gesamte Subsection "The MEPP Controversy and Scope Limitations"** mit Grinstein & Linsker-Kritik und Martyushev-Bedingungen | EN Z.468-476 | P0 |
| **Racemate-Paragraph** in Resolvent-Perspektive (illustrierendes Beispiel fuer L-/D-Enantiomere) | EN Z.493 | P1 -- ist in DE vorhanden in Z.479, aber dort integriert statt separiert |

| Was in DE anders formuliert ist | DE-Zeile | EN-Zeile | Bewertung |
|---|---|---|---|
| Pross-DKS: "maximieren" vs. "satisfy the resolvent stability criterion" | DE 464 | EN 466 | **Inhaltlich kritisch:** DE reflektiert alte Maximierungs-Logik |
| Conclusion: "Maximierung dissipierter Energie" vs. "stability-selected dissipation" | DE 490 | EN 504 | **Inhaltlich kritisch** |
| "rigoros...etabliert" vs. "developed as structural analogy" | DE 479 | EN 491 | DE zu stark |

### FST-III: DE liegt hinter EN zurueck (und hat Zusatzmaterial)

| Was fehlt in DE | EN-Zeile | Schwere |
|---|---|---|
| **(Conjecture)** im Proposition-Titel | EN Z.271 | P0 |
| **Status-Paragraph** zur FEP-MEPP-Proposition | EN Z.305 | P0 |
| **Gesamte Subsection "MEPP: Scope and Limitations"** mit Martyushev-Referenz und Einschraenkungen | EN Z.603-612 | P0 |
| **Circularity-Caveat** Paragraph | EN Z.612 | P0 |
| "MEPP as stability filter" als Punkt 4 in "Was ist wirklich neu?" | EN Z.587 | P1 |

| Was DE zusaetzlich hat (nicht in EN) | DE-Zeile | Bewertung |
|---|---|---|
| **Ausfuehrlicher Hunt-These-Abschnitt** mit 3 Subsections (Hunt-These, MOND, Kosmologische Evolution) | DE Z.472-495 | EN hat dies auf 4 Zeilen Verweis gekuerzt. DE muss entscheiden: kuerzen oder behalten. |
| **Phase 4: Kosmische Integration (Spekulativ)** | DE Z.516-518 | Nur in DE; EN hat 3 Phasen. |

---

## 5. Priorisierte Fix-Liste

### P0 -- Muss vor Upload gefixt werden

1. **[ALLE DE] Author-Block: Email durch ORCID ersetzen**
   - FST-I DE Z.35, FST-II DE Z.38, FST-III DE Z.36
   - `\texttt{lukisch@users.noreply.github.com}` ersetzen durch `ORCID: \url{https://orcid.org/0009-0005-7296-1534}`

2. **[ALLE DE] "rigoros...etabliert" durch "structural analogy"-Sprache ersetzen**
   - FST-I DE Z.424, Z.532
   - FST-II DE Z.479
   - FST-III DE Z.589
   - Systematisch "rigoros" durch "als strukturelle Analogie entwickelt" ersetzen, wo EN "structural analogy" sagt.

3. **[FST-I DE] Terminology-Paragraph einfuegen** (nach Z.76)
   - Definition von "resolvent-stable" als strukturelle Analogie
   - Uebersetzung von EN Z.78

4. **[FST-I DE] Martyushev (2010) MEPP-Einschraenkungen einfuegen** (nach Z.222)
   - Uebersetzung des letzten Satzes von EN Z.224
   - `\bibitem{Martyushev2010}` in Bibliographie ergaenzen

5. **[FST-II DE] Subsection "MEPP-Kontroverse und Geltungsbereich" einfuegen** (nach Z.464)
   - Uebersetzung von EN Z.468-476

6. **[FST-III DE] "(Konjektur)" in Proposition-Titel einfuegen** (Z.268)
   - `[FEP-MEPP-Dualit\"at via Markov-Decke]` aendern zu `[FEP-MEPP-Dualit\"at via Markov-Decke (Konjektur)]`

7. **[FST-III DE] Status-Paragraph einfuegen** (nach Z.287)
   - Uebersetzung von EN Z.305

8. **[FST-III DE] Subsection "MEPP: Geltungsbereich und Einschraenkungen" einfuegen** (nach Z.589)
   - Uebersetzung von EN Z.603-612, inklusive Circularity-Caveat

9. **[FST-II DE] Pross-DKS Absatz aktualisieren** (Z.464)
   - "weil sie die Entropieproduktion maximieren" ersetzen durch Resolvent-Stabilitaets-Formulierung

10. **[FST-II DE] Conclusion aktualisieren** (Z.490)
    - "Maximierung dissipierter Energie" ersetzen durch "stabilitaetsselektierte Dissipation"

### P1 -- Sollte vor Upload gefixt werden

11. **[ALLE DE] Nicht-escaped Umlaute konsistent machen**
    - FST-I DE: 13 Stellen (Z.71, 272, 379, 408, 416, 417, 424, 480, 511, 524, 532)
    - FST-II DE: 16 Stellen (Z.47, 61, 74, 163, 234, 236, 338, 342, 356, 376, 378, 406, 479, 490)
    - FST-III DE: 4 Stellen (Z.118, 293, 589)
    - Entweder alle auf `\"a`-Stil umstellen (konsistent) oder `\usepackage[utf8]{inputenc}` beibehalten und ALLE Umlaute als UTF-8 schreiben (dann auch die bereits escaped).

12. **[FST-I DE Z.76] "strukturell notwendige Instanz" abschwaechen** zu "strukturelle Analogie"

13. **[FST-I EN Z.555] Duplikat-Satz in Conclusion entfernen**
    - Zweiten "If successful"-Satz behalten (hedged), ersten streichen.

14. **[FST-III DE] "Was ist wirklich neu?" Punkt 4 aktualisieren** (Z.577)
    - Von "Kosmologische Verbindungen" zu "MEPP als Stabilitaetsfilter"

15. **[FST-III DE] Cosmological Extensions vs. EN abstimmen** (Z.472-495)
    - Entscheidung: ausfuehrlichen Hunt-These-Text beibehalten (dann auch in EN ergaenzen) oder auf EN-Niveau kuerzen (4 Zeilen Verweis auf Framework-Paper)

16. **[FST-III DE] Phase 4 entfernen oder an EN angleichen** (Z.516-518)
    - EN hat 3 Phasen, DE hat 4. Angleichen.

17. **[FST-I DE Z.532] Numerische alpha-scan-Ergebnisse ergaenzen**
    - EN Z.536 hat detaillierte Zahlen (98.8%, viable window etc.)

### P2 -- Wuenschenswert (Cleanup)

18. **[ALLE] Bibitem-Keys bereinigen:** Irreführende Jahreszahlen in Keys
    - `DynamicalOrigin2005` (Crooks 2007), `Lindblad2015` (Lindblad 1976), `MEPP2023` (M&S 2006), `MFG2025` (L&L 2007), `QuantumFoundations2023` (Kotecha 2024), `QuantumCognition2015` (B&B 2012), `SchmitzReview2011` (Plasson 2007)
    - Nicht funktionsbrechend, aber verwirrend fuer Wartung

19. **[ALLE] Zukunftsdaten in bibitems pruefen:** `Sandhu et al.(2026)`, `Safron et al.(2025)` -- als "forthcoming" oder "preprint" markieren falls noch nicht publiziert

20. **[FST-I EN Z.60] "mathematically rigorous" abschwaechen** zu "mathematically grounded"

---

## Zusammenfassung

Die EN-Versionen wurden systematisch ueberarbeitet mit:
- Abschwaechen von "rigoros" zu "structural analogy"
- Einfuegen von MEPP-Einschraenkungen (Martyushev 2010, Scope & Limitations)
- Hinzufuegen von Circularity-Caveats und Status-Paragraphen
- Kuerzen spekulativer Abschnitte (Hunt-These)

Die DE-Versionen haben diese Aenderungen **nicht nachvollzogen**. Es gibt **10 P0-Probleme**, die vor einem Upload gefixt werden muessen. Die kritischsten sind:
1. Fehlende MEPP-Einschraenkungen/Caveats in allen 3 DE-Papers
2. Inkonsistente "rigoros"-vs.-"structural analogy"-Sprache in allen 3 DE-Papers
3. Email statt ORCID in allen 3 DE-Papers
4. Fehlender "(Conjecture)"-Status in FST-III DE

**Geschaetzter Aufwand:** 2-3 Stunden fuer P0+P1 Fixes.
