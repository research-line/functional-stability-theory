# TODO -- FST Unified Framework

> Stand: 2026-03-15 (ALLE 6 substanziellen TODOs geschlossen + 5 triviale Fixes von 2026-03-14)

## Offene substanzielle TODOs (aus Review 2026-03-14)

**ALLE 6 substanziellen TODOs am 2026-03-15 geschlossen (Claude Opus 4.6 Session).**

### Mathematische Luecken

- [x] **REV-1: RFEP Definitionsbereich spezifizieren (Section 2)**
  FIX: Remark 2.2 (Domain of definition) eingefuegt. Fordert $e^{H_\sigma} \in L^1(\mu_\sigma)$
  explizit und begruendet warum dies in allen Instanziierungen erfuellt ist.

- [x] **REV-2: Lemma 3.1 (Euler-Maclaurin Rank-Finiteness) begruenden (Section 3.1)**
  FIX: Vollstaendige Beweisskizze eingefuegt. Argument: Boundary-Korrektur B_p haengt von
  2p Datenpunkten ab (Randwerte von a^{(2k-1)}), jeder Bernoulli-Term ist Rang-1-Korrektur.
  Remainder R_N -> 0 asymptotisch. Spektrale Abweichung = finite-rank boundary correction.

- [x] **REV-3: Theorem 4.2 (Pattern A Universality) praezisieren (Section 4)**
  FIX: Von "Theorem" zu "Principle" umdeklariert. "if and only if" zu "if" abgeschwaecht.
  Konverse als "conjectural" markiert mit Verweis auf Proposition 2.3.

### Strukturelle Luecken

- [x] **REV-4: Section 5 (Primary Implementation FST-RH) ausbauen**
  FIX: Von 6 Zeilen auf 3 Subsections erweitert: (5.1) Arithmetic Kernel K_sigma mit
  Hub-and-Spoke-Geometrie, (5.2) Shift Parity Lemma, (5.3) Resolvent Gap Theorem,
  plus Pattern A Instantiation Summary.

- [x] **REV-5: Conclusion abschwaechen / praezisieren (Section 6)**
  FIX: "definitive resolution" durch "unified structural approach" ersetzt. Klare
  Abgrenzung zwischen bewiesenem (RH) und programmatischem Teil (YM, NS, TU, DE).
  Benennung der spezifischen offenen Hypothesen pro Problem.

- [x] **REV-6: Formaler Zusammenhang RFEP -> Pattern A herstellen**
  FIX: Proposition 2.3 (RFEP implies Pattern A structure) mit Beweisskizze eingefuegt.
  Kumulanten-Expansion zeigt: 0. Ordnung neutral (deterministisch), 1. Ordnung degeneriert
  (lineare Korrektur cancelt), 2. Ordnung dominant (Var(H^{(2)}) > 0 bricht Symmetrie).

## Triviale Fixes (bereits angewendet, 2026-03-14)

- [x] **FIX-1:** Tippfehler K5: "emergent stays" -> "emergent states"
- [x] **FIX-2:** Abkuerzungen TU (Turbulence) und DE (Dark Energy) eingefuehrt
- [x] **FIX-3:** ORCID im Author-Block ergaenzt
- [x] **FIX-4:** Bibliographie: Titel der FST-RH Papers an tatsaechliche Titel angepasst,
  DOI-Placeholders mit TODO-Kommentaren versehen (Papers noch nicht auf Zenodo publiziert)
- [x] **FIX-5:** Unbenutztes Makro `\Ksym` auskommentiert
- [x] **FIX-6:** Referenzen \cite{} fuer Null-Tail und Rank-Finiteness in Section 5 ergaenzt

## Erledigte Anpassungen (aus RH-Analyse 2026-03-14)

### Pattern A Generalisierung

- [x] **RH-FW-1: Pattern A von "Rank-2 Gauge" zu "Second-Order Resolvent Dominance" generalisieren**
  FRAMEWORK_ABSTRACT.md: Pattern A Triple auf "Second-Order Resolvent Dominance" umgeschrieben.
  Drei-Stufen-Schema: Leading neutral -> First-order flat -> Second-order dominant + gauge-neutralized.
  -> Erledigt: FRAMEWORK_ABSTRACT.md Abschnitt "Pattern A" aktualisiert

- [x] **RH-FW-2: FST_Unified_Framework.tex aktualisieren**
  "Rank-2 Gauge" durch "Second-Order Resolvent Dominance" ersetzt in Section 5 (Primary Implementation).
  -> Erledigt: FST_Unified_Renormalized_Energy_Framework.tex

- [x] **RH-FW-3: Pattern A Universality Theorem formulieren**
  Theorem 4.2 (Pattern A Universality) mit 3 Bedingungen (neutral leading branch,
  first-order flatness, second-order dominance) und Instantiierungstabelle (RH/YM/NS/TU/DE).
  -> Erledigt: FST_Unified_Renormalized_Energy_Framework.tex Section 4

### Euler-Maclaurin als universales Werkzeug

- [x] **RH-FW-4: Euler-Maclaurin-Endpunktanalyse dokumentieren**
  Subsection 3.1 mit Lemma (Euler-Maclaurin Rank-Finiteness) und 5 isomorphen Endpunkt-Strukturen.
  -> Erledigt: FST_Unified_Renormalized_Energy_Framework.tex Subsection 3.1

### Stabilitaet durch Resolvent-Asymmetrien

- [x] **RH-FW-5: Neue Subsection "Universal Architecture"**
  Subsection 4.1 "Universal Architecture: Stability through Controlled Resolvent Asymmetries".
  Erklaert warum naive Ansaetze scheitern, numerische Evidenz traege ist, alle Probleme
  "fast stabil" aussehen, und das Free-Energy-Framework universell greift.
  -> Erledigt: FST_Unified_Renormalized_Energy_Framework.tex Subsection 4.1

### Querbezuege aktualisieren

- [x] **RH-FW-6: Bibliographie erweitern**
  Neue Referenzen: Geiger 2026 (FST-RH I-III), Kato 1995, Slepian-Pollak 1961.
  -> Erledigt: FST_Unified_Renormalized_Energy_Framework.tex (thebibliography)

- [x] **RH-FW-7: Ordnungsschema-Tabelle in Paper einfuegen**
  Remark 4.3 (Resolvent Order Schema) mit O(1)/O(1/L)/O(1/L^2)-Tabelle.
  -> Erledigt: FST_Unified_Renormalized_Energy_Framework.tex Section 4
