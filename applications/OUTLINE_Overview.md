# OVERVIEW: Fractal Game Theory -- A Scale-Invariant Principle from Particles to Cosmology

## Ziel
Uebersichtsartikel, der das gesamte Programm in einem Paper zusammenfasst.
Referenziert CRM I-IV (DOI: 10.5281/zenodo.18728936) als Beweis der kosmologischen Skala.
Zentrales formales Werkzeug: Renormierungsgruppe (RG) als Bruecke zwischen
Spieltheorie und Skalenuebergaengen.

## Journal: Foundations of Physics / Entropy

---

## Gliederung

### 1. Introduction
- Das Fine-Tuning-Problem als Symptom fehlender Einheit
- Variationsprinzipien als roter Faden der Physik (Fermat, Hamilton, Lagrange)
- These: Ein einzelnes Optimierungsprinzip (MEPP + Nash-GG) erklaert
  Strukturbildung auf allen Skalen -- und die Renormierungsgruppe liefert
  den mathematischen Beweis der Skaleninvarianz
- Ueberblick ueber das Programm (CRM I-IV + FST I-III)

### 2. The Universal Game: Axioms and Framework
- Axiom 1: Anfangsbedingung (maximaler Gradient, minimale Struktur)
- Axiom 2: Entwicklungsprinzip (MEPP -- Maximum Entropy Production Principle)
- Formale Definition: Spieler, Strategien, Payoff auf jeder Skala
- Tabelle: Skalenuebergreifende Zuordnung (Level 0-2c)

### 2b. The Universal Stability Criterion
**[NEU -- formaler Kern des gesamten Programms]**

- 2b.1 Scale-Tripel: Fuer jede Skala L definiere:
  - Zustand: x_L (Konfiguration / Populationsanteile / Feldparameter)
  - Spiel: G_L = (S_L, Pi_L) (Strategieraum + Payoff-Funktion)
  - Dynamik: x_dot_L = D_L(x_L; Pi_L)

- 2b.2 Nash/ESS-Stabilitaetskriterium (einheitlich):
  x_L* ist stabil, wenn es ein lokal attraktiver Fixpunkt der Dynamik
  ist UND gegen Invasion/Abweichung robust bleibt.
  -> Stabiliaet ist messbar: Fixpunkt + Robustheit (Jacobian-Spektrum)
  -> Keine "gefuehlte" Stabilitaet, sondern Lyapunov-Kriterium

- 2b.3 Gemeinsamer Payoff-Proxy:
  - Thermo-Proxy: Entropieproduktion S_dot oder freie Energie-Destruktion
  - Game-Proxy: mittlerer Payoff Pi_bar oder Fitness f_bar
  - **FST-Kernsatz:** Unter offenen Nichtgleichgewichtsbedingungen korreliert
    langfristig stabile Nash/ESS-Struktur mit hoher S_dot
    (oder hoher integrierter Dissipation ueber Zeit).

### 3. The Renormalization Group as Mathematical Engine
**[NEU -- zentrales Theorie-Kapitel]**

- 3.1 RG in der Physik: Warum Skalen zusammenhaengen
  - Wilson (1971): Kadanoff-Blocktransformation
  - Fixpunkte, relevante/irrelevante Operatoren, Universalitaetsklassen
  - Kritische Phaenomene: Am kritischen Punkt ist das System skaleninvariant
    -> fraktale Struktur emergiert

- 3.2 RG-Spieltheorie: Block-Spiele
  - Gitter von N Spielern mit lokaler Interaktion (Auszahlungsmatrix A)
  - Kadanoff-Schritt: Fasse Bloecke von b^d Spielern zu "Super-Spielern" zusammen
  - Effektive Auszahlungsmatrix A' = R_b[A] (Renormierungstransformation)
  - Fixpunkt: A* = R_b[A*] -> das Spiel sieht auf jeder Skala gleich aus
  - Nash-GG des Mikro-Spiels -> Nash-GG des Makro-Spiels (emergent)
  - Referenzen: Szabo & Fath (2007) "Evolutionary games on graphs",
    Hauert & Szabo (2005), Perc & Szolnoki (2010)

- 3.2b RG als formale Abbildung mit Invarianten
  - R: (x_L, Pi_L, D_L) -> (x_{L+1}, Pi_{L+1}, D_{L+1})
  - Invariante A: Fixpunktstruktur (Anzahl/Typ stabiler Attraktoren)
  - Invariante B: Best-response-Geometrie (Vorzeichenstruktur des Jacobians um x*)
  - Invariante C: Potential-/Lyapunov-Existenz (ob monotone Groesse existiert)
  - **RG-Forderung:** R erhaelt mindestens eine der Invarianten A-C
  - Das ist die Stelle, wo "Metapher" zur "Theorieform" wird

- 3.3 RG-Fluss und laufende Kopplungen
  - Beta-Funktion: beta(g) = dg/d(ln mu) beschreibt wie Kopplungen
    mit der Skala laufen
  - In der Spieltheorie: beta(A) = dA/d(ln L) beschreibt wie die
    Auszahlungsmatrix mit der Beobachtungsskala L laeuft
  - UV-Fixpunkt (kleine Skalen) vs. IR-Fixpunkt (grosse Skalen)
  - Phasenuebergaenge = Skalenwechsel zwischen verschiedenen Nash-GG

- 3.4 Universalitaetsklassen als Nash-Klassen
  - Systeme mit verschiedener Mikrostruktur aber gleichem Makroverhalten
    -> gleiche Universalitaetsklasse
  - Hypothese: Teilchen, Chemie, Biologie, Kosmologie gehoeren zur
    SELBEN Universalitaetsklasse im Sinne der RG-Spieltheorie
  - Die "kritischen Exponenten" sind die universellen Strukturkonstanten
    (z.B. das Verhaeltnis 4/3 in CRM = kritischer Exponent?)

### 4. Level 0: Spacetime Substrate
- Raumzeit-Bits als Spieler (LQG-Spin-Netzwerke, Kausalmengen)
- Kooperatives Alignment-Spiel -> tanh-Saettigung als Mean-Field-Loesung
- Verbindung zu Jacobson (1995): GR als thermodynamische Zustandsgleichung
- **RG-Interpretation:** tanh ist die Mean-Field-Loesung des Ising-Modells.
  Das Ising-Modell hat einen bekannten RG-Fixpunkt. Die kosmologische
  Saettigung Phi_0 ist der Ordnungsparameter am Fixpunkt.

### 5. Level 1: Particles as Nash-Optimal Toolkit
- Standardmodell als Loesung eines Optimierungsproblems
- Fine-Tuning -> Entropieproduktions-Maximierung
- Quantenmechanik als Mixed-Strategy-Gleichgewicht
  - Superposition = gemischte Strategie
  - Pfadintegral = Summe ueber alle Strategien
  - Stationaere Phase = Nash-GG
  - Destruktive Interferenz = Eliminierung nicht-Nash-Strategien
- **RG-Interpretation:** Die laufenden Kopplungskonstanten des SM
  (alpha_s, alpha_EM, alpha_W) sind die Beta-Funktionen des Level-1-Spiels.
  Grand Unification bei ~10^16 GeV = UV-Fixpunkt = primordialer Nash-GG
  vor Symmetriebrechung.
- Higgs-VEV als IR-Fixpunkt des elektroschwachen Spiels
- Conjecture: 19 freie Parameter als Nash-GG (Ref: FST-I)

### 6. Level 1b: Chemistry as Autocatalytic Games
- Eigen-Schuster Hyperzyklen als kooperative Spiele
- Nash-Gleichgewichte in Reaktionsnetzwerken
- Praebiotic Chemistry: vom Thermodynamik-Gradienten zur Autokatalyse
- Chiralitaet als Symmetriebrechung durch ESS
- **RG-Interpretation:** Autokatalytische Zyklen = Block-Spieler auf
  chemischer Skala. Ein stabiler Hyperzyklus ist ein RG-Fixpunkt:
  Er sieht "von aussen" wie ein einzelner Spieler aus (Protozelle),
  obwohl er intern aus vielen kooperierenden Molekuelen besteht.
  -> Kadanoff-Block = Kompartiment = Protozelle

### 7. Level 1c: Biology as Evolutionary Stable Strategies
- Maynard Smith ESS als biologisches Nash-GG
- England (2013): Thermodynamische Selektion
- Protozellen: Von Chemie zu Leben ueber MEPP
- Molekulare + Verhaltensebene: dieselbe Logik (Ref: FST-III)
- **RG-Interpretation:** Die "Major Transitions in Evolution"
  (Szathmary & Maynard Smith 1995) sind RG-Skalenschritte:
  - Gene -> Chromosomen (Block 1)
  - Prokaryoten -> Eukaryoten (Block 2)
  - Einzeller -> Mehrzeller (Block 3)
  - Individuen -> Gesellschaften (Block 4)
  Jeder Uebergang: Spieler der unteren Ebene werden zum "Super-Spieler"
  der naechsten Ebene. Die Spiellogik (Kooperation vs. Defektion)
  bleibt invariant.

### 8. Level 2: Cosmology (CRM Results)
- Zusammenfassung CRM I-IV (mit DOI-Referenzen)
- Paper I: Delta_chi2 = -12.2 (Pantheon+)
- Paper II: Baryonisches Universum, MOND, Joint-Fit
- Paper III: R + gamma*R^2, hi_class, Delta_chi2 = -3.7
- Paper IV: Vektorsektor, a_0 = cH_0/(2pi), SPARC
- Die kosmologische Skala als empirisch validierter Ankerpunkt
- **RG-Interpretation:** alpha_M(a) ist die laufende Kopplungskonstante
  der Gravitation im CRM. Sie beschreibt wie stark die geometrische
  Rueckkopplung von der Skala abhaengt:
  - Kosmologisch (Hubble-Skala): alpha_M ~ 0.001 (schwache Modifikation)
  - Galaktisch (kpc-Skala): Vektorsektor dominiert -> MOND-Regime
  - Sonnensystem (AU-Skala): Chameleon-Screening -> GR wiederhergestellt
  Dies ist exakt ein RG-Fluss mit drei Regimen:
  UV (Sonnensystem) -> intermediär (Galaxie) -> IR (Kosmologie)

### 9. The Fractal Structure: Synthesis
- 9.1 Formaler Isomorphismus: Tabelle der Strukturgleichheit

  | Element | Teilchen | Chemie | Biologie | Kosmologie |
  |---------|----------|--------|----------|------------|
  | Spieler | Feldanregungen | Molekuele | Organismen | Geometrie+Materie |
  | Strategie | Quantenzustand | Konfiguration | Phaenotyp | Metrik |
  | Payoff | Action S | Freie Energie G | Fitness W | Lagrangian L |
  | GG | Grundzustand | Autokatalyse | ESS | Attraktor |
  | RG-Block | Hadron | Protozelle | Mehrzeller | Galaxienhaufen |
  | Laufende Kopplung | alpha_s(mu) | K_cat(T) | W(N) | alpha_M(a) |
  | Fixpunkt | Confinement | Metabolismus | Oekosystem | de Sitter |
  | Symmetriebrechung | EW -> Higgs | Racemisch -> Chiral | Asexuell -> Sexuell | Homogen -> Struktur |

- 9.2 Self-Similarity: Was ist wirklich "fraktal"?
  - Strikte Selbstaehnlichkeit: A* = R_b[A*] (exakt gleiche Spielstruktur)
  - Self-affinity: Gleiche Logik, aber skalenabhaengige Parameter
  - Hypothese: Die Natur zeigt Self-affinity -- gleiche Spiellogik (MEPP),
    aber die "Kopplungskonstanten" laufen mit der Skala (RG-Fluss)
  - Das erklaert warum auf jeder Ebene NEUE Phaenomene auftreten
    (Emergenz), obwohl die zugrundeliegende Optimierungslogik invariant ist

- 9.3 Die Deduktionskette (Bottom-Up)
  1. Raumzeit-Substrat (Alignment-Spiel) -> Geometrie
  2. Geometrie + Quantenfelder (Stabilitaetsspiel) -> Teilchen
  3. Teilchen + Thermodynamik (Kooperationsspiel) -> Molekuele
  4. Molekuele + MEPP (Autokatalysespiel) -> Leben
  5. Leben + Selektion (ESS-Spiel) -> Biodiversitaet
  6. Materie + Geometrie (Feedback-Spiel) -> Kosmologie (CRM)
  Jeder Pfeil = ein Kadanoff-Blockschritt der RG

### 9b. The FST Test Protocol
**[NEU -- einheitliches Testprotokoll fuer alle Skalen]**

- 4-Schritte-Protokoll (fuer jede Skala identisch):
  1. Spiel spezifizieren: Spieler/Strategien/Payoff Pi
  2. Dynamik waehlen: Replikator, Best-Response, Gradient Flow, MFG
  3. Stabilitaet messen: Fixpunkte + Jacobian-Spektrum / Invasionskriterium / Lyapunov
  4. Thermo koppeln: S_dot entlang Trajektorien und an Attraktoren vergleichen

- Minimal-Deliverable pro Skala:
  | Skala | Spieler | Dynamik | Stabilitaetstest | Thermo-Proxy |
  |-------|---------|---------|-----------------|--------------|
  | Protein | Aminosaeurereste | Gradient Flow | Jacobian um Fold | DeltaG_fold |
  | Hypercycle | Molekuelspezies | Replikator | Parasiten-Invasion | S_dot(Reaktionsfluss) |
  | GRN | TF + Bindestellen | Boolean/ODE | Attraktorbecken | Energie/Info-Proxy |
  | Oekosystem | Arten | Lotka-Volterra | Nash + Resilienz | Dissipationsrate |
  | Kosmologie | Geometrie+Materie | Friedmann+CRM | alpha_M Lauf | S_dot(Sterne) |

### 10. Testable Predictions
- Pro Skala mindestens eine falsifizierbare Vorhersage:

  | Skala | Vorhersage | Test |
  |-------|-----------|------|
  | Teilchen | S_dot(SM) > S_dot(SM') fuer jede 10%-Variation | Numerische Berechnung |
  | Chemie | Autokatalytische Zyklen mit hoeherem S_dot sind stabiler | In-vitro-Experiment |
  | Biologie | Metabolische Rate korreliert mit Reproduktionserfolg | Empirische Daten |
  | Kosmologie | alpha_M laeuft mit Skala (bereits getestet: CRM) | hi_class + Planck |
  | Uebergreifend | Das Verhaeltnis 4/3 (mu -> 4/3 in CRM) tritt als kritischer Exponent auch auf anderen Skalen auf | Literaturvergleich |

- **Zusaetzlich (spekulative Erweiterung):**
  | Zahlentheorie | Riemann-Nullstellen als Eigenfrequenzen eines MEPP-selektierten Operators | Spektraldeterminante = Xi-Funktion (offen) |
  -> Siehe separates Conjecture-Paper: FST-RH

- Die staerkste Vorhersage: Wenn die fraktale Struktur real ist, dann
  muessen die Uebergaenge zwischen Skalen (Major Transitions) als
  Phasenuebergaenge beschreibbar sein -- mit berechenbaren kritischen
  Exponenten.

### 11. Discussion and Limitations
- 11.1 Ist das Analogie oder Isomorphismus?
  - Ehrliche Antwort: Derzeit Self-affinity (gleiche Logik, verschiedene
    Manifestationen). Strikte RG-Selbstaehnlichkeit ist eine HYPOTHESE.
  - Was noetig ist: Explizite Berechnung der RG-Transformation fuer
    mindestens einen Skalenuebergang (z.B. Chemie -> Biologie)

- 11.2 Verhaeltnis zum Anthropischen Prinzip
  - AP sagt: "Die Parameter sind so, weil wir existieren" (tautologisch)
  - FST sagt: "Die Parameter sind so, weil sie S_dot maximieren" (testbar)
  - Entscheidender Unterschied: FST macht quantitative Vorhersagen

- 11.3 Verhaeltnis zu anderen Ansaetzen
  - Verlinde (2011): Entropische Gravitation (Teilmenge von FST Level 0+2)
  - Friston (2010): Free-Energy-Prinzip (Teilmenge von FST Level 1c)
  - Wolfram (2002): Computational Universe (aehnlicher Geist, andere Methode)
  - Tegmark (2014): Mathematical Universe (metaphysisch, nicht testbar)
  - FST vereinigt diese Ansaetze unter einem Dach

- 11.4 Was fehlt
  - Quantitative Bruecke Teilchen -> Chemie (haerteste Luecke)
  - Explizite RG-Berechnung fuer biologische Skalen
  - Stringtheorie-Kompatibilitaet (offen)

### 12. Conclusion
- Ein Prinzip (MEPP/Nash), ein Formalismus (RG), alle Skalen
- Die Natur spielt kein "neues Spiel" auf jeder Ebene -- sie spielt
  dasselbe Spiel, renormiert auf die jeweilige Beobachtungsskala
- CRM I-IV liefern den empirischen Beweis auf kosmologischer Skala
- FST I-III werden die anderen Skalen adressieren
- Forschungsprogramm: 5+ Jahre, interdisziplinaer

---

## Schluesselreferenzen (erweitert)

### Renormierungsgruppe + Spieltheorie
- Wilson, K.G. & Kogut, J. (1974). The renormalization group and the
  epsilon expansion. Physics Reports 12(2), 75-199.
- Szabo, G. & Fath, G. (2007). Evolutionary games on graphs.
  Physics Reports 446, 97-216. DOI: 10.1016/j.physrep.2007.04.004
- Perc, M. & Szolnoki, A. (2010). Coevolutionary games -- a mini review.
  BioSystems 99, 109-125.
- Hauert, C. & Szabo, G. (2005). Game theory and physics.
  American Journal of Physics 73, 405.

### Skaleninvarianz und Fraktale
- Mandelbrot, B. (1982). The Fractal Geometry of Nature. Freeman.
- Barabasi, A.L. & Albert, R. (1999). Emergence of scaling in random
  networks. Science 286, 509-512.
- West, G.B. (2017). Scale: The Universal Laws of Life, Growth, and
  Death in Organisms, Cities, and Companies. Penguin.

### Thermodynamik + Information
- Jacobson, T. (1995). Thermodynamics of Spacetime. PRL 75, 1260.
- Verlinde, E. (2011). On the Origin of Gravity. JHEP 04, 029.
- England, J.L. (2013). Statistical physics of self-replication. JCP 139.
- Friston, K. (2010). The free-energy principle. Nature Rev. Neurosci. 11.
- Dewar, R.C. (2003). Information theory and the fluctuation theorem.

### Evolution + Major Transitions
- Maynard Smith, J. (1982). Evolution and the Theory of Games. Cambridge.
- Szathmary, E. & Maynard Smith, J. (1995). The Major Transitions in
  Evolution. Freeman.
- Eigen, M. & Schuster, P. (1979). The Hypercycle. Springer.
- Kauffman, S.A. (1993). The Origins of Order. Oxford UP.
- Nowak, M.A. (2006). Evolutionary Dynamics. Harvard UP.

### CRM (eigene Arbeiten)
- Geiger, L. (2026). CRM Papers I-IV. DOI: 10.5281/zenodo.18728936.
  GitHub: github.com/lukisch/cfm-cosmology
