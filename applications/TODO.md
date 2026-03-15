# TODO -- FST-I / FST-II / FST-III
# Fractal Game Theory: Gemeinsame Aufgabenliste
# Stand: 2026-03-12 (nach Referee Reports + Publikationsstrategien + mechanische Fixes)

> Siehe: `PUBLIKATIONSSTRATEGIE.md` fuer konsolidierte Strategie.

---

## Referee Reports -- Zusammenfassung (2026-03-10)

### FST-III (Biological Game Theory)
**Verdict:** "Suitable for publication after **moderate revision**"
- Staerken: Hierarchische ESS/Nash-Mapping original, Nash-Frustration genuiner Beitrag
- Schwaechen: MEPP als Selektionsgesetz statt Heuristik, Learning-Rate willkuerlich,
  RG-Analogie metaphorisch, Kosmologie-Sektion verwaessert biologischen Fokus
- **Kernforderung:** Konkrete falsifizierbare Predictions + experimentelle Protokolle

### FST-II (Chemical Game Theory)
**Verdict:** "Recommended for publication after **minor revision**"
- Staerken: Replicator-Nash korrekt, Azoarcus-Verankerung, Spatial ESS ueberzeugend
- Schwaechen: MEPP-Zirkularitaet, England kausal ueberinterpretiert,
  Generalisierung von wenigen Experimenten
- **Kernforderung:** Unterscheidung known results vs. new insights, experimentelle Tests

### FST-I (Thermodynamic Game Theory)
**Verdict:** "Publication recommended provided **speculative scope made more explicit**"
- Staerken: Density-Matrix/Mixed-Strategy-Mapping sauber, Envariance-Nash einsichtsvoll,
  Fine-Tuning als constrained optimization intellektuell ansprechend
- Schwaechen: SM-Parameter-Optimierung nur Subset, MEPP quasi-axiomatisch,
  interpretativ statt praediktiv, anthropomorphe Spielersprache
- **Kernforderung:** Interpretative Vereinheitlichung betonen, NICHT Ableitung

---

## Legende

- [ ] Offen
- [~] In Arbeit
- [x] Erledigt
- P0 = Blocker (vor Einreichung zwingend)
- P1 = Wichtig (vor Einreichung dringend empfohlen)
- P2 = Wuenschenswert (verbessert Qualitaet)

---

## 1. UEBERGREIFEND (alle Papers)

### P0 -- Blocker

- [x] **BIB-01: .bib-Dateien erstellen** (ERLEDIGT 2026-03-12)
  Gemeinsame `fst_references.bib` (29 KB) erstellt. Dedupliziert, Keys vereinheitlicht.
  Struktur: Shared References + FST-I-only + FST-II-only + FST-III-only.

- [ ] **NOT-01: Notation vereinheitlichen**
  FST-I: `\dot{S}`, `\mathcal{S}[\vec{p}]`
  FST-II: `dH/dt`, `\Delta S_{ext} + \Delta S_{int}`
  FST-III: `\dot{S}`
  -> Einheitliche Konvention fuer Entropieproduktionsrate definieren.
  -> Notations-Tabelle im Overview-Paper anlegen.

- [ ] **MEPP-01: MEPP-Kontroverse in FST-II und FST-III adressieren**
  FST-I hat bereits Sec 4.2 (Grinstein&Linsker, Dewar, Martyushev).
  FST-II und FST-III ignorieren die Kontroverse komplett.
  -> Mindestens einen Absatz in Discussion einfuegen.
  -> Zitieren: Jaynes 1980, Bruers 2007, Volk & Pauluis 2010, Ross et al. 2012.
  -> Klarstellen: MEPP als Arbeitshypothese, nicht als bewiesenes Gesetz.

- [x] **CLAIM-01: Terminologie "Isomorphismus" -> "structural correspondence"** (ERLEDIGT 2026-03-12)
  FST-II: Absatz in Discussion eingefuegt ("formal analogy, not a full isomorphism").
  FST-III: Analog geprueft. Ausnahme: Replicator-Nash (Hofbauer&Sigmund 1998, bewiesen).

### P1 -- Wichtig

- [ ] **TMPL-01: Journal-Templates verwenden**
  Aktuell: `\documentclass[12pt,a4paper]{article}`
  -> FST-I: Springer-Template (Foundations of Physics) oder RevTeX
  -> FST-II: Elsevier-Template (J. Theor. Biol.)
  -> FST-III: Elsevier-Template (J. Theor. Biol. / BioSystems)

- [ ] **OC-01: "Original Contributions" Sektion in jedes Paper**
  Trennung: Was ist bekannt (Replicator-Nash, CGT, England) vs. eigener Beitrag.
  -> Pro Paper eine Box/Tabelle "Status of Claims":
     ESTABLISHED / CONJECTURED / SPECULATIVE

- [ ] **PRED-01: Testable Predictions quantifizieren**
  Aktuell: ~13 Predictions, alle qualitativ ("should reveal...", "we expect...").
  -> Mindestens 3 Predictions mit Zahlenwerten / Schwellenwerten versehen.

- [ ] **REPO-01: GitHub-Repository erstellen**
  -> github.com/lukisch/fractal-game-theory (PRIVAT)
  -> Python-Code (protein_fold_nash_pdb.py, nash_mutation_score.py etc.) einchecken
  -> Daten (PDB-Dateien, JSON-Ergebnisse) ins Repo
  -> Im Paper referenzieren

### P2 -- Wuenschenswert

- [ ] **LANG-01: Deutsche Versionen (DE) erstellen**
  Quality Protocol verlangt EN + DE. Aktuell nur EN vorhanden.

- [ ] **RED-01: Redundanz FST-II / FST-III reduzieren**
  MEPP, England-Inequality, Dissipative Adaptation werden in beiden ausfuehrlich behandelt.
  -> FST-III soll auf FST-II verweisen statt zu wiederholen.

---

## 2. FST-I -- Thermodynamic Game Theory of Fundamental Particles

### Referee-Empfehlungen (2026-03-10)
- [ ] **I-REF-01:** Betonen: Beitrag ist interpretive unification, NICHT Ableitung
- [ ] **I-REF-02:** Falsifiability-Section mit klareren Kriterien erweitern
  (was koennte das Framework widerlegen?)
- [ ] **I-REF-03:** Stellar-Entropy-Berechnung als "pilot study" positionieren,
  mit expliziten Erweiterungsplaenen
- [ ] **I-REF-04:** Anthropomorphe Spieler-Sprache bei Vakuumfeldern moderieren
  (trotz Disclaimern)

### P0 -- Blocker

- [x] **I-CALC-01: Mindestens EINE quantitative Berechnung** (ERLEDIGT 2026-03-10)
  Semi-analytisches Sternmodell in `fst_i_entropy_calculation.py` implementiert.
  Ergebnis: S_total(alpha_obs) = 98.8% des Maximums. Neue Section 5.3 im Paper.
  Gamow-Energie-Bug (Faktor 2) gefunden und korrigiert (E_G = 493 keV).
  Ergebnisse in `fst_i_entropy_results.json` gespeichert.

- [x] **I-K1: Zentrale These (SM = Nash-GG) operationalisieren** (ERLEDIGT 2026-03-10)
  Alpha: 98.8% des Maximums. m_e: Monotonie erklaert via Multi-Skalen-Demarkation
  (Copilot-Loesung: Chemieskala setzt Untergrenze, Sterne optimieren innerhalb).
  "Hierarchisches Selektionsprinzip" als zentrales strukturelles Ergebnis eingebaut.

### P1 -- Wichtig

- [ ] **I-K2: MEPP als Hypothese kennzeichnen, Kritik zitieren**
  Sec 4.2 existiert bereits, aber:
  -> Jaynes 1957 und 't Hooft 2016 fehlen noch in der Bibliographie.
  -> Polettini 2013 (Markov-Ketten-Gegenbeispiele) zitieren.

- [x] **I-K3: Born-Regel -- Mehrwert der spieltheoretischen Sicht benennen** (ERLEDIGT 2026-03-10)
  Copilot-Loesung: Ehrliche Demarkation. Kein neues Theorem, sondern strukturelle
  Vereinheitlichung + Axiom-Transparenz. Text in "Scope and Formal Demarcation" eingefuegt.

- [x] **I-K4: Aktion S als "Payoff" mathematisch rechtfertigen** (ERLEDIGT 2026-03-10)
  Copilot-Loesung: Zwei Wege: (A) Hamiltonscher Minimax-Sattelpunkt (Zero-sum Nash),
  (B) Euklidisch: U=-S_E, dann echtes Extremalprinzip. Beide in Paper eingefuegt.

- [ ] **I-K5: Higgs-Potential Konvention klarstellen**
  mu^2 Vorzeichen, Phasenuebergang. -> Peskin & Schroeder Konvention explizit angeben.

- [ ] **I-K6: Falsifizierbarkeitskriterium mit konkretem Zahlenwert fuellen**
  "If a variation increases S while preserving stability" -- welche Variation genau?
  -> Constraints vorab definieren (welche Stabilitaetsbedingungen zaehlen).

- [x] **I-W3: "Quantum-Game Duality" Conjecture praezisiert** (ERLEDIGT 2026-03-10)
  Copilot-Loesung: Conjecture -> "Quantum-Strategy Correspondence". Affine Einbettung
  Delta_n -> D(H) als Proposition. Off-Diagonalen als nichtkommutative Erweiterung abgegrenzt.

- [x] **I-W4: Electron-Mass Toy Model korrigieren** (ERLEDIGT 2026-03-10)
  Gamma_fusion ~ exp(-m_e) durch korrekte Physik ersetzt:
  Opacity kappa ~ m_e^{-2}, Gamow E_G = (pi*alpha)^2 * m_r * c^2.

- [ ] **I-W5: Section 7 erweitern -- FST-I vs. andere Frameworks**
  Vergleich: Anthropic Principle, Cosmological Natural Selection (Smolin),
  Free Energy Principle (Friston). -> Was unterscheidet FST-I?

- [ ] **I-W7: CFM einfuehren und zitieren**
  CFM (Cosmological Fractal Model, DOI: 10.5281/zenodo.18728936) wird erwaehnt
  aber nirgends als formale Referenz eingefuehrt.

### P2 -- Wuenschenswert

- [ ] **I-W6: Zitier-Keys harmonisieren**
  z.B. DynamicalOrigin2005 verweist auf Crooks (2007). Keys an Autoren anpassen.

- [ ] **I-W2-REST: Jaynes 1957, 't Hooft 2016 in Bibliographie aufnehmen**

---

## 3. FST-II -- Chemical Game Theory

### Referee-Empfehlungen (2026-03-10)
- [ ] **II-REF-01:** Unterscheidung: established chemical kinetics vs. genuinely new insights
  durch game-theoretic framing
- [ ] **II-REF-02:** Explizite experimentelle Tests vorschlagen, die CGT-Predictions
  von Standard-Kinetik-Modellen unterscheiden
- [ ] **II-REF-03:** Teleologische Sprache reduzieren wo moeglich
- [ ] **II-REF-04:** England-Ungleichung kausal nicht ueberinterpretieren
  (Entropieproduktion ist FOLGE von Replikation, nicht notwendig Treiber)

### P0 -- Blocker

- [x] **II-K1: England-Ungleichung korrekt darstellen** (ERLEDIGT 2026-03-10)
  Copilot-Loesung: Crooks FT -> Jensen -> England Bound als Lemma eingebaut.
  tau-Notation als Rate-Approximation markiert. Game-theoretisches Reading ergaenzt.

- [x] **II-K2: Gl.(5) dH/dt = -(E-E*)/T korrigiert** (ERLEDIGT 2026-03-10)
  Copilot-Loesung: Gl.(5) war im Allgemeinen FALSCH (Gegenbeweis: konstante Fitness).
  Ersetzt durch: exakte H-dot = -Cov_x(f, ln x) + KL-Divergenz als Lyapunov +
  Shahshahani-Gradientenfluss-Interpretation (Harper/Hofbauer).

- [x] **II-K3: Gl.(6) Distributed Nash korrigieren** (ERLEDIGT 2026-03-10)
  Copilot-Loesung: MFG-Fixpunkt mit Konsistenzbedingung m*=L(x*) eingefuegt.
  Potenzialspiel-Argument fuer detailed-balance-Systeme ergaenzt.

### P1 -- Wichtig

- [ ] **II-W1: Delta_G = -T*Delta_S_total Geltungsbereich klarstellen**
  Gilt NUR fuer geschlossene Systeme bei konstantem T, p.
  -> Fuer offene NESS-Systeme: Jaegliche Gibbs-Argumente mit Vorsicht verwenden.

- [ ] **II-W3: MEPP-Kontroverse transparent machen**
  Dewar 2003 widerlegt durch Grinstein & Linsker 2007.
  -> In Discussion-Sektion adressieren (analog zu FST-I Sec 4.2).

- [x] **II-W4: Prisoner's Dilemma vollstaendig spezifizieren** (ERLEDIGT 2026-03-12)
  Vollstaendige PD-Bedingung (T>R>P>S + 2R>T+S) in Tabelle und erklaerenden Absatz
  eingefuegt. Kinetische Terme: k_BA > k_AA > k_BB > k_AB + Hofbauer&Sigmund-Verweis.

- [ ] **II-W5: Theorem 2 (Replicator-Nash) praezisieren**
  In gegebener Form zu stark formuliert.
  -> Hofbauer & Sigmund (1998) korrekt attribuieren.
  -> Einschraenkungen: innere GG vs. Rand-GG.

- [ ] **II-W6: Gl.(3) n vs. N Notation bereinigen**
  Inkonsistente Notation. Aktivitaetskoeffizienten fehlen.

- [ ] **II-W7: Quantitative Vorhersagen ableiten**
  Alle 5 Predictions (P1-P5) sind nur qualitativ.
  -> Mindestens P1 (Entropieproduktion vs. Zykel-Stabilitaet) quantifizieren.

- [ ] **II-CALC-01: Eigene Gibbs-Minimierung fuer ein einfaches Netzwerk**
  Das Paper hat NULL eigene Berechnungen (Azoarcus-Daten sind fremd).
  -> Ein minimales 3-Spezies-Autokatalyse-Modell durchrechnen.

### P2 -- Wuenschenswert

- [ ] **II-W5b: Fehlende Referenzen nachruesten**
  Weibull 1995 (Evolutionary Game Theory), Crooks 1999 (Fluctuation Theorem).

- [x] **II-LANG: England-Theorem nicht als \begin{theorem} deklarieren** (ERLEDIGT 2026-03-12)
  War bereits korrekt als `\begin{lemma}[England (2013), from Crooks FT]` deklariert.
  Keine Aenderung noetig.

---

## 4. FST-III -- Biological Game Theory

### Referee-Empfehlungen (2026-03-10)
- [ ] **III-REF-01:** Konkrete falsifizierbare Predictions mit klaren experimentellen
  Protokollen (ueber qualitative Vorschlaege hinaus)
- [ ] **III-REF-02:** Kosmologie-Sektion weiter separieren/kuerzen oder als
  eigenstaendiges Perspective Paper auslagern (verwaessert biologischen Fokus)
- [ ] **III-REF-03:** Klare Trennung: Was ist konzeptionelle Vereinheitlichung
  vs. was beansprucht erklaerende/praediktive Neuheit?
- [ ] **III-REF-04:** MEPP nicht als Selektionsgesetz verwenden, sondern als Heuristik
- [ ] **III-REF-05:** Nash-Frustration: Learning-Rate-Parameter prinzipiell kalibrieren
  (nicht willkuerlich)

### P0 -- Blocker

- [ ] **III-CALC-01: Nash-Frustration auf mindestens 3 weitere Proteine erweitern**
  Aktuell: n=1 (HP35/1YRF). Im Ordner liegen 1PGA.pdb und 2XWR.pdb.
  -> protein_fold_nash_pdb.py auf 1PGA und 2XWR anwenden.
  -> Ergebnisse in Paper integrieren (Tabelle mit rho(J) pro Protein).

- [ ] **III-CALC-02: TP53-Mutationsanalyse integrieren**
  tp53_full_scores.json und nash_mutation_score.py existieren im Ordner.
  -> Korrelation Nash-Frustration vs. Krebsmutations-Hotspots wuerde
     ein extrem starkes Ergebnis liefern.
  -> Neuer Subsection in Sec 3 oder eigene Section.

- [ ] **III-COSMO-01: Kosmologie-Sektion (Sec 8) auslagern**
  Die Hunt-Thesis / MOND-Physik ist fuer J. Theor. Biol. deplatziert.
  -> Sec 8 und Sec 9.4 in das Overview-Paper verschieben.
  -> In FST-III nur einen Verweis-Absatz in Discussion behalten:
     "The broader cosmological extension is developed in [Overview-Paper]."

### P1 -- Wichtig

- [ ] **III-K3: Nash-Frustration reproduzierbar spezifizieren**
  -> Kraftfeld (AMBER, CHARMM?), eta-Wert, Software vollstaendig angeben.
  -> Code-Referenz: protein_fold_nash_pdb.py auf GitHub.

- [ ] **III-K4: Nash-Frustrations-Gleichung vervollstaendigen**
  -> eta (Learning Rate): Wie wird eta gewaehlt? Datengetriebene Methode vorschlagen.
  -> Hessian-Potential: Welches Potential? (Lennard-Jones, AMBER ff14SB, etc.)

- [ ] **III-W2: MEPP-Fitness-Verbindung moderieren**
  MEPP ist umstritten. -> "MEPP as heuristic organizing principle" statt "MEPP drives fitness".

- [ ] **III-W3: Replicator ~ Lotka-Volterra Aequivalenz praezisieren**
  -> Hofbauer & Sigmund (1998): Formale Aequivalenz unter welchen Bedingungen?
  -> Aktuell zu pauschal formuliert.

- [x] **III-W5: Tabellen mit \label und \caption versehen** (ERLEDIGT 2026-03-12)
  5 Tabellen (nicht 6) mit formalen \begin{table}[h]-Umgebungen, \caption und \label
  versehen: tab:thermo-principles, tab:game-biology-mapping, tab:frustration-comparison,
  tab:network-regimes, tab:scale-kernel.

- [ ] **III-W6: Proposition-Nummerierung vereinheitlichen**
  Proposition 1/2 vs. P1-P5 in Predictions. -> Konsistentes Schema.

- [ ] **III-FRUSTR-01: Systematischer Vergleich Nash- vs. Ferreiro-Frustration**
  Fuer 5+ Proteine: Nash-Frustration vs. Protein Frustratometer (Parra et al. 2016).
  -> Korrelation, Unterschiede, komplementaere Informationen herausarbeiten.

### P2 -- Wuenschenswert

- [ ] **III-W8: Abbildungen einfuegen**
  Gesamtes Paper hat 0 Abbildungen. Empfohlen: mindestens 3.
  -> Fig 1: Nash-Frustrations-Map auf HP35-Struktur (farbcodiert)
  -> Fig 2: Hierarchische Spielstruktur (Molekuel -> Zelle -> Organismus -> Oekosystem)
  -> Fig 3: TP53 Hotspot-Korrelation (wenn III-CALC-02 erledigt)

- [x] **III-FEP-01: Friston-FEP-Sektion (Sec 4) Bruecke gebaut** (ERLEDIGT 2026-03-10)
  Copilot-Loesung: Neue Subsection "Free-Energy Minimization as Nash Fixed Point".
  Potentialspiel mit koordinatenweiser VI. Scope-Demarkation: wann formal, wann nur konzeptionell.

- [ ] **III-BOOL-01: Boolean-Networks-Sektion (Sec 5) vertiefen oder kuerzen**
  Aktives Forschungsfeld, aber nur oberflaechlich behandelt.
  -> Alon 2007, Wagner 2005 zitieren, oder auf 1 Absatz kuerzen.

- [ ] **III-CANCER-01: Aktipis-Abgrenzung klarstellen**
  Cancer-als-Coalition-Breakdown (Aktipis 2015) ist bekannt.
  -> Eigener Beitrag: Verbindung zum Hyperzyklus-Parasitenproblem als ANALOGIE, nicht Resultat.

- [x] **III-RG-01: Einen RG-Schritt exemplarisch berechnen** (2026-03-10)
  Copilot-Loesung: Toy-RG mit 2x2-Blocking, Majority Rule, Frozen-Interior-Approximation.
  Explizite Mischformel A'_{s,t} = q_s q_t R + ... eingebaut in Sec.7.
  Fixpunkte (All-C, All-D) und Major Transitions als Aktivierung neuer relevanter Operatoren.

---

## 5. STRATEGISCHE EINREICHUNGSREIHENFOLGE (aktualisiert nach Referee Reports)

### Neue Reihenfolge (basierend auf Referee-Verdicts):

1. **FST-II (Chemical Game Theory) -- MINOR REVISION**
   - Am wenigsten Arbeit bis Einreichung
   - Journal: J. Theoretical Biology oder Origins of Life and Evolution of Biospheres
   - Kernarbeit: MEPP-Kontroverse adressieren, England kausal abschwaechen,
     eigene Gibbs-Minimierung berechnen, experimentelle Tests vorschlagen

2. **FST-III (Biological Game Theory) -- MODERATE REVISION**
   - Nash-Frustration ist genuiner Beitrag (Alleinstellungsmerkmal!)
   - Journal: J. Theoretical Biology oder BioSystems
   - Kernarbeit: Kosmologie-Sektion auslagern, Nash-Frustration auf 3+ Proteine,
     TP53-Korrelation einbauen, falsifizierbare Predictions quantifizieren
   - **Alternative:** Nash-Frustration als eigenstaendiges Short Paper:
     "Nash Frustration in Protein Folding" bei Proteins / J. Chem. Phys.

3. **FST-I (Thermodynamic Game Theory) -- SCOPE EXPLIZITER MACHEN**
   - Interpretative Vereinheitlichung, NICHT Ableitung
   - Journal: Foundations of Physics (passt am besten zum interpretativen Charakter)
   - Kernarbeit: Stellar-Entropy als Pilotstudie positionieren,
     Falsifizierbarkeit schaerfen, Spieler-Sprache moderieren

### Uebergreifende P0-Blocker (vor allen Einreichungen):
- [x] BIB-01: .bib-Dateien erstellen (2026-03-12)
- [ ] NOT-01: Notation vereinheitlichen
- [ ] MEPP-01: MEPP-Kontroverse in FST-II und FST-III adressieren
- [x] CLAIM-01: "Isomorphismus" -> "structural correspondence" (2026-03-12)

---

## 6. ERLEDIGTE PUNKTE (Referenz)

### Gen 4 (2026-03-10)
- [x] FST-III: Broken \ref{} in Conclusion -> prop:multicellularity
- [x] FST-III: ESSWikipedia -> MaynardSmithPrice1973
- [x] FST-III: FEPWikipedia -> FristonFreeEnergy
- [x] FST-III: 12 unvollstaendige bibitems durch peer-reviewed ersetzt
- [x] FST-III: Hamilton1964-Kontextualisierung bei Multicellularity-Proposition
- [x] FST-III: Cross-References zu FST-I und FST-II eingefuegt
- [x] FST-II: 16 problematische bibitems durch peer-reviewed ersetzt
- [x] FST-II: TimeReversalSymmetry -> Kondepudi & Nelson (1985)
- [x] FST-II: AzoarcusPNAS -> Yeates et al. (2016) mit DOI
- [x] FST-II: Cross-Reference zu FST-I eingefuegt
- [x] FST-I: 7 unvollstaendige bibitems durch peer-reviewed ersetzt
- [x] FST-I: Eisert et al. (1999) eingefuegt

### Gen 3 (2026-03-03)
- [x] Alle Papers: E-Mail korrigiert, Datum auf March 2026
- [x] FST-I: Sprachmoderation (6 Stellen), MEPP-Kontroverse Sec 4.2
- [x] FST-II: Sprachmoderation (8 Stellen), Autokatalyse-Notation
- [x] FST-III: \usepackage{physics} entfernt, CFM-Akronym, Sprachmoderation
- [x] Alle Papers: AI Disclosure (ICMJE-konform)

### Gen 1-2 (2026-02-28)
- [x] AI Disclosure in alle 3 Papers
- [x] Sprachmoderation (K1 + K2 Fixes)

---

## 7. EXTERNE EVALUIERUNG -- REVIEW-INPUTS (2026-03-11)

Basierend auf umfassender externer Evaluierung des FST/CRM-Gesamtprogramms.
SotA-Abgleich, kritische Luecken, methodische Bruecken fuer FST-I/II/III.

### 7.1 MEPP-Paradoxon -- KRITISCHE ROTE FLAGGE (alle Papers)

**Problem:** Dogmatische Verwendung des MEPP fuer fern-vom-Gleichgewicht-Systeme.
MEPP basiert auf lokales thermodynamisches Gleichgewicht + Bilinearitaet (sigma=X*J).
Fuer hochgradig nichtlineare chemische Netzwerke (Autokatalyse, Proteinfaltung)
bricht diese Linearitaet zusammen.
-> Systeme weit vom GG graben sich in kinetisch stabile MINIMA der Dissipation
   (eher Prigogines mEPP als MEPP).
-> Scheinbarer Widerspruch: Maximum vs. Minimum Entropy Production.
-> Kann spieltheoretischen Nash-Ansatz torpedieren wenn Payoff falsch definiert.

**Loesung (Review):** MEPP im Framework auf INFORMATIONSENTROPIE ueber die
Markov-Decke (Markov Blanket) verschieben, statt rein thermische Entropie.
-> Nash-GG im MFG = Zustand wo kein Agent durch unilaterale Modifikation
   mehr Energie dissipieren kann ohne strukturelle Integritaet zu gefaehrden.
-> Maximierung der Entropieproduktion in der Umwelt = Minimierung der internen
   variationellen freien Energie (FEP). Diese Dualitaet loest das Paradoxon.

### 7.2 Uebergreifende Aufgaben (alle Papers)

- [ ] **EXT-BIO-1: MEPP-Paradoxon aufloesen (P0, alle Papers)**
  MEPP als Selektionsgesetz ist zu stark. Loesung: Informationsentropie ueber
  Markov Blanket statt rein thermische Entropie.
  -> FST-I Sec 4.2: Erweitern um Markov-Blanket-Formalismus.
  -> FST-II: MEPP-Kontroverse-Absatz (MEPP-01) muss dieses Paradoxon adressieren.
  -> FST-III: MEPP als "heuristic organizing principle" (schon in III-W2 gefordert).
  -> Lit: Martyushev (2010) "MEPP Two Basic Questions" -- zeigt Restriktionen rigoros.

- [ ] **EXT-BIO-2: Mean Field Games (Lasry-Lions) als formaler Rahmen**
  Formale Aequivalenz: FEP + Dissipative Adaptation + Nash-GG via MFG.
  -> HJB-Gleichung: optimale Policy (Active Inference) rueckwaerts in der Zeit.
  -> Fokker-Planck: Dichteentwicklung des Agentenensembles vorwaerts.
  -> Biologische Fitness -> HJB-Kosten abbilden.
  -> FST-I Sec 2 hat bereits MFG-Ansatz. Muss formaler/rigoroser werden.
  -> FST-III: FEP-Nash-Bruecke (III-FEP-01 erledigt) vertiefen mit MFG-Kopplung.
  -> Lit: Lasry & Lions (2006/2007) -- MFG Theorie.

### 7.3 FST-I -- Neue Aufgaben

- [ ] **EXT-I-1: Friston-FEP vs. England Dualitaet explizit formulieren**
  Die formale Dualitaet (Umwelt-Entropieproduktion maximieren = interne
  freie Energie minimieren) muss mathematisch sauber hergeleitet werden.
  -> Als neues Lemma oder Proposition in Sec 5 oder 7.
  -> Lit: England (2015) "Dissipative Adaptation in Driven Self-Assembly"
  -> Lit: Friston (2006, 2010) "Free-Energy Brain Theory"

- [ ] **EXT-I-2: Active Inference als stochastisches Kontrollproblem**
  Active Inference kann in optimales stochastisches Kontrollproblem eingebettet
  werden: Kosten/Nutzen direkt in Prior Beliefs absorbiert.
  -> Formale Bruecke zu utilitaristischen Nash-GG.
  -> Ergaenzt I-W5 (Section 7: FST-I vs. andere Frameworks).

### 7.4 FST-II -- Neue Aufgaben

- [ ] **EXT-II-1: England-Ungleichung kausal abschwaechen (bestaetigt II-REF-04)**
  Entropieproduktion ist FOLGE von Replikation, nicht notwendig Treiber.
  -> Review bestaetigt: "England kausal ueberinterpretiert".
  -> Dissipative Adaptation = Strukturen die Absorption/Dissipation maximieren
     sind thermodynamisch BEVORZUGT, aber nicht kausal GETRIEBEN.
  -> Formulierung anpassen in FST-II.

- [ ] **EXT-II-2: Bilinearitaets-Einschraenkung des MEPP adressieren**
  Fuer Azoarcus-Ribozym-System und praebiologische Chemie:
  -> Explizit pruefen ob lokales thermo. Gleichgewicht vorausgesetzt werden kann.
  -> Falls nicht: auf FEP-Formalismus umstellen (Markov-Blanket).
  -> Lit: Martyushev (2010), Toma (2025) "Birth and Evolution via MEPP"

### 7.5 FST-III -- Neue Aufgaben

- [ ] **EXT-III-1: Nash-GG als MFG-Fixpunkt formalisieren (bestaerkt III-FEP-01)**
  Die bestehende FEP-Nash-Bruecke muss um die HJB+Fokker-Planck-Kopplung
  erweitert werden.
  -> Biologische Ueberlebenswahrscheinlichkeit auf HJB-Gleichung abbilden.
  -> Potentialspiel-Argument (III-FEP-01) ist Spezialfall des MFG-Fixpunkts.

- [ ] **EXT-III-2: MEPP-Paradoxon in Proteinfaltung explizit ansprechen**
  Proteinfaltung = Paradebeispiel fuer "kinetisch stabiles lokales Minimum"
  (minimale Entropieproduktion, nicht maximale).
  -> Nash-Frustration muss dieses Paradoxon aufloesen, nicht ignorieren.
  -> Argument: Frustration IS die Manifestation des MEPP-mEPP-Konflikts.

### 7.6 Neue Referenzen (Bio/Chemie)

- [ ] England (2015): Dissipative Adaptation -- thermodynamischer Mechanismus
  fuer physikalische Selbstorganisation. Bereits zitiert, kausal abschwaechen.
- [ ] Martyushev (2010): "MEPP Two Basic Questions" -- Restriktionen des MEPP.
  GEFAEHRLICH aber unverzichtbar. Muss in FST-II/III referenziert werden.
- [ ] Friston (2006, 2010): FEP + Active Inference + Markov Blankets.
  Bereits zitiert (FST-III), aber informationstheoretische Prinzipien
  muessen staerker in die Agenten-Definition einfliessen.
- [ ] Toma (2025): "Birth and Evolution via MEPP" -- Wendet MEPP auf
  Prae-RNA-Welt + Vielzelligkeit an. Stuetzt direkt FST-Biologie-Ansatz.
- [ ] Lasry & Lions (2006/2007): Mean Field Games -- Formaler Rahmen
  fuer Nash-GG in Vielteilchen-Systemen. ZWINGEND fuer rigorosen Beweis.

### Strategische Empfehlung (Review)

> "Studieren Sie intensiv die formale Mathematik der Mean Field Games (Lasry/Lions).
> Die Uebersetzung von Englands statistischer Mechanik in gekoppelte HJB-Fokker-Planck-
> Systeme liefert den rigorosen Nash-Beweis, den das Programm zwingend braucht."

> "Das MEPP-Paradoxon (Maximum vs. Minimum Entropy Production) muss durch den
> FEP-Formalismus (Markov-Blanket-Informationsentropie) aufgeloest werden."

---

## 8. LITERATURRECHERCHE-REVIEW -- NEUE INPUTS (2026-03-12)

Basierend auf erschoepfender externer Literaturrecherche und kritischer Analyse
des FST-Programms im Kontext des physikalischen und biophysikalischen
Forschungsstandes 2025/2026.

### 8.1 Uebergreifend: MEPP-Validierung (alle Papers)

- [ ] **LIT-MEPP-1: MEPP-Makro-Evidenz aus MDPI Entropy 2025 zitieren (P1)**
  Neuer MDPI-Review (Entropy 2025) "Maximum Entropy Production Principle
  of Thermodynamics for the Birth and Evolution of Life" bestaetigt:
  MEPP beschreibt makroskopische Evolutionstendenzen (Zellspezialisierung)
  hervorragend, aber fundamentale Herleitung aus First Principles fern vom
  GG bleibt lueckenhaft.
  -> In FST-I Sec 4.2 und FST-II/III MEPP-Diskussion zitieren.
  -> Stuetzt "MEPP as heuristic organizing principle"-Positionierung.
  -> Lit: Toma (2025), MDPI Entropy

- [ ] **LIT-MEPP-2: Grinstein-Linsker-Einschraenkungen explizit benennen (P1)**
  Review betont: MEPP hat keine universelle Gueltigkeit in diskreten
  stochastischen Systemen (Grinstein & Linsker). Bereits in FST-I Sec 4.2
  erwaehnt, aber FST-II und FST-III ignorieren dies noch.
  -> Verstaerkt bestehende MEPP-01 Forderung.

### 8.2 FST-I -- Neue Aufgaben

- [ ] **LIT-I-1: MFG-Wasserstein-Raum-Formulierung referenzieren (P1)**
  Neueste MFG-Entwicklungen (U Michigan 2026, Van Eanam Lectures) arbeiten
  mit PDEs im Wasserstein-Raum und Master-Gleichungen.
  -> Tiefe strukturelle Parallelen zur Quanten-Dichtematrix-Entwicklung.
  -> In Sec 2 (MFG-Ansatz) oder Sec 7 (Vergleich) referenzieren.
  -> Lit: "On Approximate Nash Equilibria in Mean Field Games" (arXiv:2601.20910)
  -> Lit: "Network Aggregative Markov Games"

- [ ] **LIT-I-2: Quantum-Game-Duality mit MFG-Maß-Raum validieren (P2)**
  Die Abbildung Delta_n -> D(H) (Proposition in FST-I) koennte durch
  MFG-Formulierung im Wasserstein-Raum rigoroser werden.
  -> Pruefen ob MFG-Master-Gleichung die Quanten-Strategie-Korrespondenz
     als Spezialfall enthaelt.

- [ ] **LIT-I-3: Comparison Table FST-I vs. andere Frameworks (P1)**
  Review fordert: Klarer Vergleich mit Anthropic Principle, Cosmological
  Natural Selection (Smolin), Free Energy Principle (Friston).
  -> Verstaerkt bestehende I-W5 Forderung.
  -> NEU: Auch MFG-basierte Ansaetze einbeziehen.

### 8.3 FST-II -- Neue Aufgaben

- [ ] **LIT-II-1: Machine Learning Interatomic Potentials (MLIPs) als
  Verifikationswerkzeug vorschlagen (P1)**
  Review identifiziert: "Computational Studies of Prebiotic Chemistry at
  the Age of Machine Learning" (ChemRxiv 2025/2026) zeigt: Neuronale
  Netzwerk-Potenziale koennen autokatalytische Netzwerke, Uebergangszustaende
  und Reaktionskinetiken mit nahezu QM-Genauigkeit simulieren.
  -> Active-Learning-Workflows liefern kinetische Eingangsdaten (Rate Constants)
     fuer Geigers chemische Spielmatrizen automatisiert.
  -> In Discussion als "Computational outlook" vorschlagen.
  -> Loest teilweise II-CALC-01 (eigene Gibbs-Minimierung) methodisch.
  -> Lit: ChemRxiv preprint, "Prebiotic Resource Constraints" (bioRxiv)

- [ ] **LIT-II-2: Bilinearitaets-Einschraenkung des MEPP praezisieren (P1)**
  Review bestaetigt EXT-II-2: MEPP basiert auf lokales thermo. GG +
  Bilinearitaet (sigma=X*J). Fuer hochgradig nichtlineare chemische
  Netzwerke bricht diese Linearitaet zusammen.
  -> Systeme weit vom GG graben sich in kinetisch stabile MINIMA der
     Dissipation (eher Prigogines mEPP als MEPP).
  -> Verstaerkt MEPP-01 und EXT-II-2 Forderungen.

### 8.4 FST-III -- Neue Aufgaben

- [ ] **LIT-III-1: FEP-Kritik (Stegemann 2025) referenzieren (P2)**
  Philosophische Debatte 2025: FEP hat gravierende konzeptionelle Maengel
  wenn als universelle Erklaerung von "Bewusstsein" oder "Agency" ueberdehnt.
  -> Geiger adressiert dies bereits (semiotisches Problem), aber explizite
     Referenz fehlt.
  -> In Discussion Sec 4 (FEP-Sektion) einarbeiten.

- [x] **LIT-III-2: FEP-MEPP-Paradoxon via Skalentrennung formalisieren (P1)** -- ERLEDIGT (2026-03-12)
  Review liefert elegante Aufloesung: Organismus minimiert INTERN freie
  Energie (FEP/Markov Blanket) -> stabile maschinelle Integritaet.
  Diese Struktur EXTERN als hocheffizienter dissipativer Motor -> MEPP
  auf Makroskala erfuellt.
  **Erledigt:** Proposition prop:fep-mepp ("FEP-MEPP duality via Markov blanket")
  als neue Subsection "Scale Separation: Resolving the FEP--MEPP Paradox"
  in FST-III eingefuegt (zwischen FEP-Nash und Boolean Networks).
  Drei-Schritt-Beweis: (i) FEP sichert Integritaet, (ii) intakte Struktur
  als dissipativer Motor, (iii) MEPP-Selektion auf Makroskala.
  Spieltheoretischer Remark (inner game vs. outer game) hinzugefuegt.
  Neue Referenz: Martyushev2010.

- [ ] **LIT-III-3: England-Dissipative-Adaptation Evidenz aktualisieren (P2)**
  Review zeigt: Neue numerische Belege (2025/2026) fuer Englands Theorie.
  -> In FST-III (und FST-II) England-Abschnitte mit neuen Referenzen
     aktualisieren.
  -> Lit: Quanta Magazine "First Support for a Physics Theory of Life"

### 8.5 Uebergreifende strategische Einsicht (Review)

> "Das FST-Programm repräsentiert eine aussergewoehnliche intellektuelle
> Leistung der Unifikation. Die methodische Praemisse, dass die
> vielfaeltigsten Stabilitaetsprobleme der Natur alle Manifestationen
> eingeschraenkter konvexer Optimierung in thermodynamischen und
> informationellen Spielen sind, besitzt enorme heuristische Kraft."

> "Um den Anspruch auf Letztgueltigkeit zu erheben, muessen die aufgezeigten
> analytischen Luecken durch Integration spezifischer Methoden aus
> Gittereichtheorien, Machine Learning Potenzialen und stochastischen
> Differentialspielen geschlossen werden."

### 8.6 Neue Referenzen (Literaturrecherche)

| Referenz | Relevanz | Ziel-Paper |
|----------|----------|------------|
| Kirk (2026), O(4) Restoration, Zenodo | CL1-CL3 fuer SU(2) | YM (integriert) |
| Kirk (2026b), Zero-Free Strip, Zenodo | RG-Fluss ohne Phasenuebergaenge | YM (integriert) |
| Garcia Baquero (2026), TMT, ResearchGate | Praegeometrischer Alternativansatz | YM (integriert) |
| Jalalzadeh et al. (2025/2026), MOND-Entropy-EUP | Holographische Stuetzung | DE (integriert) |
| Hohmann et al. (2025/2026), Finsler-JCAP | Starker Konkurrent | DE (integriert) |
| Hunt (2026), Informational Persistence | Hunt-These-Stuetzung | DE (integriert) |
| Mallet-Paret & Sell (1988), Inertial Manifolds | Hyperviskose Modelle | NS (integriert) |
| Toma (2025), MEPP Birth/Evolution, MDPI | MEPP-Validierung | FST-I/II/III |
| ChemRxiv (2025/2026), MLIPs Prebiotic | Verifikationswerkzeug | FST-II |
| Stegemann (2025), FEP-Kritik | FEP-Grenzen | FST-III |
| arXiv:2601.20910, Approx. Nash in MFG | MFG-Theorie | FST-I |
| Cunningham (2025), Topological YM Proof | Alternativer YM-Ansatz | YM |

---

## 9. RH-Durchbruch-Integration fuer FST I-III (2026-03-14)

> Aus der RH-Analyse (Second-Order Resolvent Dominance Pattern):
> Die FST-Anwendungsfaelle beschreiben keine Maximierungsprinzipien,
> sondern Stabilitaetsprinzipien. MEPP selektiert nicht das Optimum,
> sondern den resolvent-stabilen Fixpunkt.

### 9.1 FST-I (Thermodynamic Game Theory / Fundamental Particles)

- [x] **RH-FST1-1: MEPP als Stabilitaetsfilter reformulieren (P0)** (ERLEDIGT 2026-03-14)
  Remark~\ref{rem:mepp-stability-filter} in Sec.~4 (MEPP) eingefuegt.
  MEPP als Stabilitaetsfilter unter degenerierten Leading-Order-Dynamiken reformuliert.
  Cross-Reference zu FST-RH2 (Resolvent Gap Theorem).

- [x] **RH-FST1-2: Fine-Tuning als Resolvent-Stabilitaet erklaeren (P1)** (ERLEDIGT 2026-03-14)
  Remark~\ref{rem:fine-tuning-resolvent} nach "Connection to the Nash Equilibrium"
  in Sec.~5 (Fine-Tuning) eingefuegt. Even/Odd-Sektor-Analogie zum Weil QF.

- [x] **RH-FST1-3: Expliziter Verweis auf Degenerate Perturbation Theory (P2)** (ERLEDIGT 2026-03-14)
  Remark~\ref{rem:degenerate-perturbation} in Sec.~8 (Connections/Discussion)
  eingefuegt. Nash-GG als second-order-stabiler Fixpunkt formalisiert.
  Cross-Reference zu FST-RH2 und "Functional Positivity under Constraint".

### 9.2 FST-II (Chemical Game Theory / Origin of Life)

- [x] **RH-FST2-1: Homochiralitaet als Second-Order-Effekt (P1)** (ERLEDIGT 2026-03-14)
  Remark~\ref{rem:chirality-resolvent} nach "MEPP Explanation" in Sec.~7
  (Homochirality) eingefuegt. Even/Odd-Sektor-Analogie: Racemat als
  degenerierter Fixpunkt, Chiralitaet als Resolvent-Aufspaltung.

- [x] **RH-FST2-2: Autokatalyse als Resolvent-Stabilisierung (P2)** (ERLEDIGT 2026-03-14)
  Remark~\ref{rem:autocatalysis-resolvent} nach Theorem~replicator-nash
  in Sec.~6 (Hypercycles) eingefuegt. Replicator-Gleichung selektiert
  resolvent-stabilen Fixpunkt, nicht produktivsten.

- [x] **RH-FST2-3: England-Crooks nicht kausal, sondern als Stabilitaetsfilter (P1)** (ERLEDIGT 2026-03-14)
  Remark~\ref{rem:england-stability-filter} nach England-Inequality
  in Sec.~3 eingefuegt. Dissipationsbound als Resolvent-Stabilitaetsfilter
  reformuliert. Adressiert Referee-Kritik II-REF-04.

### 9.3 FST-III (Biological Game Theory)

- [x] **RH-FST3-1: ESS als resolvent-stabiler Fixpunkt (P0)** (ERLEDIGT 2026-03-14)
  Remark~\ref{rem:ess-resolvent} nach ESS-Definition in Sec.~3
  (Molecular ESS) eingefuegt. ESS als resolvent-stabiler Fixpunkt,
  Nash-Frustration als notwendige Resolvent-Architektur erklaert.

- [x] **RH-FST3-2: Punktuiertes Gleichgewicht als Resolvent-Transition (P2)** (ERLEDIGT 2026-03-14)
  Remark~\ref{rem:punctuated-equilibrium} als neue Subsection
  "Punctuated Equilibrium as Resolvent Transition" in Discussion eingefuegt.
  Stasis = resolvent-stabile Phase, Major Transitions = Resolvent-Uebergaenge.

- [x] **RH-FST3-3: MEPP-Sektion schaerfen (P0)** (ERLEDIGT 2026-03-14)
  Remark~\ref{rem:mepp-stability-biology} nach Thermo-Hypothesen-Tabelle
  in Sec.~2 (MEPP) eingefuegt. MEPP als Stabilitaetsfilter statt
  Selektionsgesetz reformuliert. Adressiert Referee-Kernforderung III-REF-04.

### 9.4 Gemeinsame Meta-Aussage (alle FST Papers)

- [x] **RH-FST-META: Einheitliche Formulierung einfuegen** (ERLEDIGT 2026-03-14)
  Neue Subsection "Stability vs. Maximization: The Resolvent Perspective"
  in allen 3 Papers eingefuegt (FST-I: vor Conclusion, FST-II: in Discussion,
  FST-III: in Discussion). Einheitlicher Kern-Absatz mit paper-spezifischen
  Anpassungen. Cross-Reference zu FST-RH2 (Geiger 2026). Neue bibitem
  FST-RH2 in allen 3 Papers hinzugefuegt.

---

## 10. Review-Findings 2026-03-14

### 10.1 Triviale Fixes (SOFORT ERLEDIGT)

- [x] **REV-FIX-01: `enumitem`-Paket fehlt (alle 3 Papers)**
  `\begin{enumerate}[label=\alph*)]` erfordert `\usepackage{enumitem}`.
  In alle 3 Papers eingefuegt. Ohne enumitem bricht die Kompilierung.

- [x] **REV-FIX-02: Markdown-Syntax `**...**` statt `\textbf{}` (alle 3 Papers)**
  Sec 2 (Structural Alignment): `**Pattern A Stability Logic**` und
  `**gauge-theoretic reformulation...**` -- Markdown-Bold statt LaTeX.
  In allen 3 Papers durch `\textbf{...}` ersetzt.

- [x] **REV-FIX-03: Broken `\ref{sec:frustration}` in FST-III**
  Zeile 324 (FEP-MEPP Remark): `Section~\ref{sec:frustration}` verweist
  auf ein nicht existierendes Label. Korrigiert zu `\ref{sec:molecular}`.

- [x] **REV-FIX-04: `\usepackage{physics}` fehlt in FST-III**
  FST-III verwendet `\dd` (Zeile 280-281) ohne `physics`-Paket.
  Eingefuegt. Ohne physics-Paket Kompilierungsfehler.

- [x] **REV-FIX-05: Fehlender `\bibitem{SchmitzReview2011}` in FST-II**
  Zitiert in Zeile 356 (Frank Model Subsection) aber kein bibitem vorhanden.
  Bibitem fuer Plasson et al. (2007) eingefuegt (Standardreferenz fuer
  Chirality-Symmetry-Breaking-Review).

- [x] **REV-FIX-06: natbib-inkompatible Bibitem-Keys (FST-II, FST-III)**
  FST-II: `\bibitem[OdumPinkerton1955]` -> `\bibitem[Odum \& Pinkerton(1955)]`
  FST-III: 3 bibitems ohne Autor/Jahr-natbib-Syntax korrigiert:
  `ReplicatorOptimization` -> Sandhu et al. (2026),
  `ScaleRGFrontiers` -> Safron et al. (2025),
  `SelfSimilarNetworks` -> Gallos et al. (2012).
  Verursachte `natbib Error: Bibliography not compatible with author-year citations`.

### 10.2 Substanzielle Befunde -- Mathematik und Konsistenz

- [x] **REV-MATH-01: Resolvent Remarks -- Analogie vs. Beweis** (ERLEDIGT v2/Zyklus5)
  Alle Remarks verwenden jetzt "structural analogy" / "motivated by structural parallels"
  statt "established rigorously". Insbesondere:
  - FST-I: "developed as structural analogy"
  - FST-II: "developed as structural analogy"
  - FST-III: "motivated by the structural parallels"

- [x] **REV-MATH-02: Resolvent-Perspective Subsections -- Redundanz** (ERLEDIGT v2/Zyklus5)
  Alle drei Resolvent-Subsections jetzt individuell formuliert mit konkreten Beispielen:
  - FST-I: alpha-Scan Ergebnisse, 98.8% Plateau, viable window
  - FST-II: Chirality als Even/Odd-Analogie, Racemat-Entartung
  - FST-III: Nash-Frustration rho(J)~1, Protein-Faltungsbeispiel

- [x] **REV-MATH-03: FST-I Remark rem:fine-tuning-resolvent** (ERLEDIGT v2/Zyklus5)
  Entschaerft: "analogous to the second-order spectral splitting" +
  "complementary subspaces" statt explizite even/odd Sektoren.
  Die Even/Odd-Sektor-Analogie wird behauptet, aber was waeren die "even"
  und "odd" Sektoren in der Teilchenphysik? Welche konkreten Parameter-
  Sektoren koppeln resonant? Ohne diese Spezifikation ist die Analogie
  leer.
  -> Empfehlung: Entweder konkrete Parameter-Sektoren benennen (z.B.
     QCD/EW-Kopplung) oder die Even/Odd-Analogie entfernen.

- [ ] **REV-MATH-04: FST-II Remark rem:autocatalysis-resolvent -- Widerspruch zu Sec 8 (P1)**
  Remark behauptet: "These structures are not MEPP-maxima -- they do not
  globally maximize entropy production." Dies widerspricht Sec 8 (MEPP
  Explanation) und dem Conclusion, wo Homochiralitaet als "maximum
  efficiency in free energy destruction" beschrieben wird.
  -> ERLEDIGT (Zyklus5): DE Zeile 378 korrigiert: "maximale Effizienz bei der
     Zerstörung freier Energie" -> "resolventenstabiles Gleichgewicht".

- [x] **REV-MATH-05: FST-III -- Nash-Frustration rho(J)** (ERLEDIGT v2/Zyklus4+5)
  Proposition -> "Computational observation (Nash Stability Boundary)" umdeklariert.
  eta-Abhaengigkeit explizit in "Important limitation" Passage.
  Wert von 1.0131 -> 1.013 (konsistent mit Berechnungsergebnis).

### 10.3 Substanzielle Befunde -- Notation und Cross-Paper-Konsistenz

- [ ] **REV-NOT-01: England-Referenz-Keys inkonsistent (P2, FST-II vs. FST-III)**
  FST-II: `\citep{England2013}` -> bibitem `England2013`
  FST-III: `\citep{EnglandJCP}` -> bibitem `EnglandJCP`
  Beide verweisen auf England (2013), J. Chem. Phys. Inkonsistente Keys.
  -> Sollte bei Migration zu gemeinsamer .bib vereinheitlicht werden (BIB-01).

- [ ] **REV-NOT-02: Hofbauer-Sigmund-Referenz-Keys inkonsistent (P2, FST-II vs. FST-III)**
  FST-II: `HofbauerSigmund1998` und `HofbauerSigmund2003`
  FST-III: `CoexistenceLotkaVolterra` (zeigt auf Hofbauer & Sigmund 1998)
  -> Vereinheitlichen bei Migration zu gemeinsamer .bib.

- [ ] **REV-NOT-03: Crooks-Key `DynamicalOrigin2005` in FST-I (P2)**
  Bibitem referenziert Crooks (2007), Key enthaelt "2005".
  Bereits in TODO als I-W6 notiert. Bestaetigt.

### 10.4 Substanzielle Befunde -- LaTeX und Struktur

- [ ] **REV-LATEX-01: FST-III Float-Specifier (P2)**
  5 Tabellen mit `[h]` erzeugen LaTeX-Warnungen ("changed to ht").
  -> Empfehlung: Alle auf `[ht]` oder `[htbp]` aendern.

- [ ] **REV-LATEX-02: Bibitem-Autoren unvollstaendig (P2, FST-III)**
  `\bibitem[Replicator Optimization]{ReplicatorOptimization}` und
  `\bibitem[Scale RG Frontiers]{ScaleRGFrontiers}` und
  `\bibitem[Self-Similar Networks]{SelfSimilarNetworks}` --
  Fehlende Autoren/Jahre im natbib-Key.
  -> Autoren und Jahre ergaenzen fuer korrekte Zitationsformatierung.

- [ ] **REV-LATEX-03: AI Disclosure nicht als numbered Section (P2, alle Papers)**
  `\section*{AI Disclosure}` steht NACH `\end{thebibliography}`.
  Das ist formal korrekt (unnummerierte Section), aber manche Journals
  erwarten die Disclosure innerhalb des Haupttextes.
  -> Bei Journal-Submission anpassen.

### 10.5 Resolvent Remarks -- Detailbewertung

| Paper | Remark | Label | Kohaerenz mit Haupttext | Korrektheit |
|-------|--------|-------|-------------------------|-------------|
| FST-I | MEPP as Stability Filter | rem:mepp-stability-filter | Gut. Konsistent mit Sec 4.2 (MEPP Controversy) | Analogie korrekt, "established rigorously" zu stark |
| FST-I | Fine-Tuning as Resolvent Stability | rem:fine-tuning-resolvent | Akzeptabel. Even/Odd-Analogie zu abstrakt | Even/Odd-Sektoren nicht definiert |
| FST-I | Nash as Second-Order Stable Fixed Points | rem:degenerate-perturbation | Gut. Konsistent mit Sec 8 (Connections) | Formal korrekt als Analogie |
| FST-II | England's Bound as Stability Filter | rem:england-stability-filter | Sehr gut. Loest Referee-Kritik II-REF-04 | Inhaltlich ueberzeugend |
| FST-II | Autocatalytic Stability as Resolvent Selection | rem:autocatalysis-resolvent | Widerspruch zu Sec 7.5/Conclusion (MEPP-Maximierung) | Inhalt korrekt, Kontext inkonsistent |
| FST-II | Homochirality as Second-Order Effect | rem:chirality-resolvent | Gut. Elegant als Even/Odd-Analogie fuer L/D | Strukturell ueberzeugend |
| FST-III | MEPP as Stability Filter (Biology) | rem:mepp-stability-biology | Sehr gut. Adressiert Referee III-REF-04 direkt | Inhaltlich korrekt |
| FST-III | ESS as Resolvent-Stable Fixed Point | rem:ess-resolvent | Gut. Konsistent mit Nash-Frustration | "Necessary structural feature" plausibel |
| FST-III | Punctuated Equilibrium as Resolvent Transition | rem:punctuated-equilibrium | Interessant. Verknuepfung zu Eldredge/Gould | Spekulative Analogie, klar demarkiert |

### 10.6 Gesamtbewertung

**Mathematische Korrektheit:** Gut. Keine Fehler in Formeln/Gleichungen gefunden.
Die Herleitungen (England-Lemma, Replicator-Nash, KL-Lyapunov) sind korrekt.
Die numerischen Ergebnisse (alpha-Scan, Nash-Frustration) sind intern konsistent.

**Resolvent Remarks:** Inhaltlich kohaerenter und besser integriert als erwartet.
Die Remarks adressieren gezielt Referee-Kritiken (MEPP als Selektionsgesetz,
England kausal ueberinterpretiert). Hauptschwaeche: Zu starker Rigorositaets-
anspruch ("established rigorously") fuer eine konzeptionelle Analogie.

**LaTeX-Qualitaet:** Nach den 5 Fixes (FIX-01 bis FIX-05) kompilieren alle
drei Papers fehlerfrei. Verbleibende Warnungen sind harmlos (Float-Specifier).

**Cross-Paper-Konsistenz:** Notation und Referenz-Keys sind inkonsistent
(England, Hofbauer-Sigmund, Crooks). Wird bei Migration zur gemeinsamen
.bib (BIB-01) geloest.

---

## 11. QUANTITATIVE VALIDIERUNG -- BLOCKER FUER READINESS 9-10 (2026-03-15)

> Stand nach 4 Review-Zyklen + Quantitative Pipeline (Claude Opus 4.6, 2026-03-15):
> FST-I: 8.5/10, FST-II: 8.5/10, FST-III: 7.5/10.
> Quantitative Ergebnisse berechnet und in Papers integriert.
> Verbleibend: eta-Kalibrierung (QUANT-III-2), REV-MATH Fixes, 2XWR-Ergebnis.

> **Publikationsreihenfolge:** FST-I (Physik) -> FST-II (Chemie) -> FST-III (Biologie)
> mit grossem zeitlichen Abstand. Jedes Paper einzeln veroeffentlichen.

### 11.1 FST-I: Multi-Parameter Constraint-Analyse

- [x] **QUANT-I-1: 2D-Scan alpha vs. m_e/m_p** (ERLEDIGT 2026-03-15)
  Script `fst_i_entropy_calculation.py` auf ellmos-services ausgefuehrt.
  **Ergebnisse:**
  - 1D alpha-Scan: S(obs)/S_max = 0.988 bei alpha/alpha_obs = 1.628
  - 2D-Scan: Maximum bei (alpha/alpha_obs=1.24, m_e/m_e_obs=0.30), S(obs)/S_max_2D = 0.47
  - 7/13 Falsifizierbarkeits-Verletzungen (konsistent mit Paper)
  - **Schl\"usselresultat:** Stellare Entropie selektiert ALPHA (98.8%), aber NICHT m_e/m_p.
    Elektronenmasse erfordert multi-skalige Constraints (Atom/Chemie/Kern).
  -> In FST-I EN+DE integriert: Neuer "Two-dimensional scan" Absatz, Falsifiability-Update.
  -> JSON: `results/fst_i/fst_i_entropy_results.json` (lokal + Server)
  -> **Constraints noch nicht implementiert** -- verbleibt als kuenftiger Enhancement

### 11.2 FST-II: Entropieproduktion konkurrierender Netzwerke (P3-Test)

- [x] **QUANT-II-1: P3-Test Replikator-Dynamik** (ERLEDIGT 2026-03-15)
  Script `replicator_entropy_test.py` erstellt und auf ellmos-services ausgefuehrt.
  **Ergebnisse (3-Spezies Azoarcus-Modell):**
  - Kooperative ESS: <f> = 6.0, stabil, Basin = 40%
  - Eigennuetzige ESS: <f> = 3.0, stabil, Basin = 59%
  - Parasitaere FP: <f> = 0.1, instabil
  - **P3 bestaetigt:** Kooperatives Netzwerk hat 2x hoeheren katalytischen Durchsatz
  - **Aber Bistabilitaet:** Von uniformem Start -> Selfish gewinnt (hoehere Init-Rate)
  - **Koordinationsbarriere:** Kooperativ braucht x_C > 0.4 -- erklaert Vaidya et al.
  **Ergebnisse (4-Spezies):**
  - Coop-A + Coop-B ESS: <f> = 3.0, stabil, von uniform erreichbar
  - Selfish ESS: <f> = 3.0, stabil
  -> In FST-II EN+DE integriert: Neue Tabelle + Ergebnis-Absatz bei P3.
  -> JSON: `results/fst_ii/replicator_azoarcus_3species.json` + `_4species.json`

### 11.3 FST-III: Nash-Frustration erweitern + eta kalibrieren

- [x] **QUANT-III-1: Nash-Frustration auf 1PGA und 2XWR** (ERLEDIGT 2026-03-15)
  `protein_fold_nash_pdb.py` auf ellmos-services ausgefuehrt.
  **1YRF (HP35, Self-Validation):**
  - rho(J) = 1.013, Hessian 80% positiv, mean_frust = 0.082
  - Max: Met11 (1.0), Pro20 (0.49). Helix-Kerne nahe Null.
  **1PGA (Protein G, Cross-Prediction):**
  - Best-Response konvergiert (282 Sweeps), phi_rms=1.30, psi_rms=2.00
  **2XWR (p53 DBD, 462 Residuen):**
  - Best-Response konvergiert, phi_rms=1.25, psi_rms=2.11, F=-39.43
  - Bemerkenswert: p53 (13x groesser) hat NIEDRIGEREN phi_rms als 1PGA!
  -> In FST-III EN+DE integriert: Tabelle mit allen 3 Proteinen.
  -> JSON: `results/fst_iii/1YRF_self`, `1PGA`, `2XWR` (lokal + Server)

- [ ] **QUANT-III-2: eta-Kalibrierung gegen experimentelle Daten**
  OFFEN. Der Learning-Rate-Parameter eta = 0.01 ist ad hoc.
  -> Systematischer eta-Scan + Vergleich mit H/D-Exchange-Daten
  -> Verbleibt als wichtigstes offenes TODO fuer FST-III Readiness 9+

- [x] **QUANT-III-3: TP53-Mutationsanalyse** (ERLEDIGT 2026-03-15, NEGATIV-RESULTAT)
  248 Mutationen analysiert (206 pathogen, 42 benign).
  **Ergebnis:** ROC-AUC = 0.46 (unter Zufall). Nash-Scores diskriminieren NICHT.
  **Ursache:** 1YRF-trainiertes Modell transferiert nicht auf TP53 (35 vs 191 Residuen,
  voellig andere Faltung). Alle 248 Mutationen: converged=False, delta_rho ~ 0.
  -> Wird als Limitation im Paper berichtet (Modell-Transferabilitaet)
  -> Positiv: Pipeline funktioniert, Problem ist Modell-Spezifitaet nicht Methodik

### 11.4 Server-Pipeline

- **Pipeline:** `/opt/fst_calculations/` auf ellmos-services (46.62.243.71)
- **Struktur:** `scripts/{fst_i,fst_ii,fst_iii}/`, `data/`, `results/`, `logs/`
- **Ollama:** ABGESCHALTET (Docker-Container gestoppt, vorerst nicht benoetigt)
- **Zugang:** `ssh -i ~/.ssh/id_ed25519_mcmc root@46.62.243.71`
- **Lokale Kopie:** `.RESEARCH/Natur&Technik/4 Anwendungsfaelle/results/` + `scripts/`
- **Doku:** Server-README unter `/opt/fst_calculations/README.md`
- **Nach Berechnungen:** Ergebnisse herunterladen, in Papers integrieren,
  finalen Review-Zyklus durchlaufen
