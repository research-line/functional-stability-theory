# FST Folgebeweise -- Forschungs-Roadmap

## Erstellt: 2026-03-16
## Quellen: 3 LLM-Reviews (Copilot, Gemini, Deep Research Pro) + 6 Copilot-Einzelreviews + Externe Evaluierung + Literatur-Review
## Gesamt: 131 offene TODOs aus 12 Review-Dokumenten

---

## Uebersichtstabelle: Loesungsvorschlaege nach Paper und Quelle

| Paper | Copilot Gesamt | Gemini | Deep Research | Copilot Einzel | Ext. Eval. | Lit-Review | Summe |
|-------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| BSD | 5 | 3 | 6 | 5 | -- | -- | 19 |
| Dark Energy | 4 | 4 | 4 | 5 | -- | 3 | 20 |
| Hodge | 3 | 3 | 3 | 5 | -- | -- | 14 |
| NS Skeleton | 3 | 5 | 4 | 6 | 3 | 2 | 23 |
| NS LogDistance | 1 | 2 | 2 | 10 | -- | -- | 15 |
| P vs NP | 3 | 3 | 3 | 4 | -- | -- | 13 |
| Turbulenz | 3 | -- | -- | 4 | 1 | -- | 8 |
| Yang-Mills | 3 | 2 | 3 | 8 | 3 | 3 | 22 |
| Framework | 3 | 3 | 4 | 2 | 1 | -- | 13 |
| Cross-Paper | -- | -- | 1 | 2 | 1 | -- | 4 |
| **Gesamt** | **28** | **25** | **30** | **51** | **9** | **8** | **~151 (131 nach Deduplizierung)** |

> Hinweis: Einige Ansaetze wurden von mehreren Quellen unabhaengig vorgeschlagen. Die deduplizierte Gesamtzahl ist 131.

---

## 1. BSD Positivity (7.5/10)

### Kritische Luecken
- **L1:** Higher Gross-Zagier fuer Rang >= 2 (keine expliziten Heegner-Punkte)
- **L2:** Euler-Systeme fuer Rang >= 2 (Kolyvagin-Hierarchie endet bei Rang 1)
- **L3:** Sha-Endlichkeit fuer Rang >= 2 (Teil der BSD-Vermutung selbst)
- **Review-Finding BSD-1:** Theorem 4.1 ist zirkulaer (setzt offene Vermutungen als Axiome voraus)
- **Review-Finding BSD-2:** Step 2 der Proof-Strategie nicht bewiesen
- **Review-Finding BSD-3:** Exclusion E2 hat logische Luecke

### Loesungsansaetze (konsolidiert aus allen Reviews)

| # | Ansatz | Quelle | Prioritaet | Machbarkeit |
|---|--------|--------|:-----------:|:-----------:|
| 1 | Heegner-Simplex (r-dim. Objekt statt r Punkte auf Kuga-Sato/Shimura) | Copilot | Hoch | Mittel |
| 2 | Diagonal-Zyklen als matrixwertige Euler-Systeme (Exterior-Power-Systeme) | Copilot | Hoch | Mittel |
| 3 | Selmer-Splitting ueber CM/Endomorphismen (zerlegbare Selmer-Gruppen) | Copilot | Mittel | Mittel |
| 4 | Hoehere Chow-Gruppen CH^r auf Picard-Modulflaechen (derivierte Schnittpaarungen) | Gemini | Mittel | Niedrig |
| 5 | Derivierte Hecke-Algebren: Taylor-Wiles-Defekt = Rang | Gemini | Niedrig | Niedrig |
| 6 | Hodge-Arakelov-Volumina: Sha als Mass eines topologischen Volumens | Gemini | Mittel | Niedrig |
| 7 | Arakelov-motivische Volumenformel: L^(r)(E,1) als arithmetisches Volumen | Deep Research | Hoch | Niedrig |
| 8 | p-adische Deformation nicht-diagonaler Heegner-Zykel auf Shimura-Varietaeten | Deep Research | Hoch | Mittel |
| 9 | Galois-Sektoren-Zerlegung (Transfer-Operator, analog YM) | Deep Research | Mittel | Niedrig |
| 10 | Kohomologische Sweeping Processes fuer Sha (konvexe Analysis auf Galois-Kohomologie) | Deep Research | Mittel | Niedrig |
| 11 | O-minimale Strukturen fuer Sha-Endlichkeit (geometrische Topologie erzwingt Endlichkeit) | Deep Research | Hoch | Niedrig |
| 12 | Shimura-Varietaeten mit CM-Punkten fuer Rang-2 Heegner-Analoga identifizieren | Copilot Einzel | Hoch | Mittel |
| 13 | Abel-Jacobi-Abbildungen in Familien, non-degeneracy der Regulator-Matrix | Copilot Einzel | Hoch | Mittel |
| 14 | p-adische Eisenstein-/Beilinson-Flach-Elemente in Hida-Familien | Copilot Einzel | Hoch | Mittel |
| 15 | "Rank-r Kolyvagin-relations" in numerischen Beispielen (Rang-2 Kurven) testen | Copilot Einzel | Hoch | Hoch |
| 16 | Numerische Experimente: L^(r)(E,1), Regulatoren, Sha-Schaetzungen fuer Rang 2-4 | Copilot Einzel | **Sofort** | **Hoch** |
| 17 | Diagonal-Cycle + Euler-System Hybrid fuer Rang-2 CM-Kurven | Copilot/Gemini | Mittel | Mittel |
| 18 | Plektische Gross-Zagier-Formel | Gemini | Mittel | Niedrig |
| 19 | Universeller Arakelov-Regulator als Dach ueber p-adisch/archimedisch | Copilot | Mittel | Niedrig |

### Empfohlener naechster Schritt
**Numerische Experimente fuer Rang-2-Kurven** (L^(r), Regulatoren, Sha-Schaetzungen) -- hoher Ertrag, moderater Aufwand, liefert datengetriebene Hypothesen. Parallel: Theorem 4.1 als "Conditional/Framework Theorem" relabeln (BSD-1).

### Schluessel-Literatur
- Kim (2022-2024), Darmon-Rotger, Burungale-Tian (2026), Burns-Sano-Sakamoto (2021-2024)

---

## 2. Dark Energy (9.0/10)

### Kritische Luecken
- **L4:** Rigoroses RG-Matching (explizite beta_Lambda-Berechnung offen)
- **L6:** Solar-System Screening (Skalaron-Masse lokal ausreichend gross?)
- **Review-Status:** Minimalitaetslemma GESCHLOSSEN (2026-03-15), RG-Kern offen

### Loesungsansaetze (konsolidiert aus allen Reviews)

| # | Ansatz | Quelle | Prioritaet | Machbarkeit |
|---|--------|--------|:-----------:|:-----------:|
| 1 | FRG als Constraint: Kovarianz + Subtraktion + Phi-Stabilitaet => H0^2 | Copilot | Hoch | Mittel |
| 2 | "Renormierung = Eichfixierung" im FST-Sinn | Copilot | Hoch | Mittel |
| 3 | Spektrale Regularisierung (Hubble-adaptiver Regulator) | Copilot | Mittel | Mittel |
| 4 | Zwei-Skalen gamma(R): kosmologisch weich, lokal steif => Screening | Copilot | **Sofort** | **Mittel** |
| 5 | Konkrete f(R)-Klasse: thermodynamisches Minimum + PPN-Grenze + keine Ghosts | Copilot | **Sofort** | **Mittel** |
| 6 | Topologisches RG-Matching: quartische Divergenz als topologische Invariante (Atiyah-Singer) | Gemini | Mittel | Niedrig |
| 7 | Holographisches Bulk-Boundary-Matching (AdS/CFT-inspiriert) | Gemini | Mittel | Niedrig |
| 8 | Dichtekorrelierter FST-Constraint: m_s koppelt nichtlinear an T^mu_mu | Gemini | Hoch | Mittel |
| 9 | Geometrische Phasenuebergaenge: Sonnensystem als spontaner Symmetriebruch | Gemini | Niedrig | Niedrig |
| 10 | Thermodynamische Holographie: UV-Divergenzen in de-Sitter-Horizont absorbiert | Deep Research | Hoch | Niedrig |
| 11 | Asymptotische Sicherheit: Vakuumenergie am UV-Fixpunkt, IR-Fluss = FST-Penalty | Deep Research | Hoch | Mittel |
| 12 | Vainshtein im Poeschl-Teller-Modell (Galileon-Terme) | Deep Research | Hoch | Mittel |
| 13 | Topologisches Screening (Aharonov-Bohm bei Dichtegradienten) | Deep Research | Niedrig | Niedrig |
| 14 | Semiklassische Effektivwirkung: zeta-Regulierung, heat-kernel, R ln(Box) R | Copilot Einzel | **Sofort** | **Mittel** |
| 15 | Beta-Funktion fuer Lambda in EFT mit f(R): Wilsonian RG + Matching | Copilot Einzel | Hoch | Mittel |
| 16 | Hu-Sawicki-artige f(R): IR = rho^2/k^2, dicht = grosses m_s (chameleon) | Copilot Einzel | **Sofort** | **Mittel** |
| 17 | Euclid/DESI/LISA Forecasts fuer eta=5/4 | Copilot Einzel | Hoch | Hoch |
| 18 | Fehlerbudget: Standard-Model-Spektren, Renorm-Scheme-Abhaengigkeit | Copilot Einzel | Mittel | Mittel |
| 19 | Finsler-Diskriminante als "smoking gun test" schaerfen (Euclid DR1 2027) | Lit-Review | Hoch | Mittel |
| 20 | MOND-Entropy-EUP Verbindung vertiefen (Jalalzadeh et al.) | Lit-Review | Mittel | Mittel |
| 21 | DVFT als verwandten Ansatz einordnen (Abgrenzung in Discussion) | Lit-Review | Niedrig | Hoch |

### Empfohlener naechster Schritt
**Semiklassische Effektivwirkung ableiten** (zeta-Regulierung, nichtlokale Terme) + **Zwei-Skalen gamma(R) konstruieren** fuer Screening. Beides zusammen loest L4+L6 auf einmal.

### Schluessel-Literatur
- Wetterich (2024), Vasak/Garcia Baquero (2026, CCGG), Hu-Sawicki (2007)

---

## 3. Hodge Positivity (7.0/10)

### Kritische Luecken
- **L1/L2:** AP = AbsHodge? (Abgrenzung zur Standard-Vermutung unklar)
- **Review-Finding HOD-4:** Proof Sketch Theorem 5.1 Step 5 ist Programm, nicht Beweis
- **Review-Finding HOD-5:** AP1 Vorzeichen-Fehler (KRITISCH -- muss (-1)^p Q_L >= 0 lauten)
- **Review-Finding HOD-1:** Proposition 2.3 logisch unvollstaendig

### Loesungsansaetze (konsolidiert aus allen Reviews)

| # | Ansatz | Quelle | Prioritaet | Machbarkeit |
|---|--------|--------|:-----------:|:-----------:|
| 1 | p-adische Slope-Constraints als staerkere AP' (Newton vs. Hodge Polygon) | Copilot | Hoch | Niedrig |
| 2 | Arakelov-Positivitaet (adelische Metriken, Heights als globaler Verstaerker) | Copilot | Hoch | Niedrig |
| 3 | Motivische Positivitaet (Variationsprinzip auf Motiven, Testfall identifizieren) | Copilot | Mittel | Niedrig |
| 4 | Prismatische Positivitaet AP4 (Scholze, perfektoide Raeume) | Gemini | Hoch | Niedrig |
| 5 | Bridgeland-Stabilitaet auf derivierten Kategorien (Positivitaet kategorisch) | Gemini | Mittel | Niedrig |
| 6 | Spektrale Luecke auf Vektorbuendeln (algebraische Zykel = Vakuumszustaende) | Gemini | Niedrig | Niedrig |
| 7 | p-adische Positivitaet via O-Minimalitaet (topologische Starrheit => Spezialisierung) | Deep Research | Hoch | Niedrig |
| 8 | Galois-Defekte: transzendente Klassen => negatives GHR-Spektrum => Diskriminante | Deep Research | Hoch | Mittel |
| 9 | Ricci-artiger Fluss auf (p,p)-Formen (algebraische Zykel = thermodynamische Minima) | Deep Research | Mittel | Niedrig |
| 10 | Voisin-Beispiele numerisch auf GHR-Negativitaet pruefen | Copilot Einzel | **Sofort** | **Hoch** |
| 11 | CM-Typen: Conjecture IX.1 via Mumford-Tate-Analyse | Copilot Einzel | Hoch | Mittel |
| 12 | Tannaka-Formalismus: Positivitaet in Mumford-Tate-Algebra uebersetzen | Copilot Einzel | Mittel | Mittel |
| 13 | VHS-Familien: monodromische Stabilitaet => algebraische Herkunft erzwingen | Copilot Einzel | Mittel | Mittel |
| 14 | Q_L-Algorithmus: numerische Evaluation unter verschiedenen Polarisierungen | Copilot Einzel | Mittel | Hoch |
| 15 | SDP-Modell fuer AP-Kegel in K3/Abelsche Varietaeten | Copilot/Gemini | Mittel | Mittel |
| 16 | Tannakische Rekonstruktion fuer CM-Faelle (Konsens 4/4 Reviews) | Copilot/Gemini | Mittel | Mittel |

### Empfohlener naechster Schritt
**KRITISCH: AP1 Vorzeichen-Fehler korrigieren** (HOD-5) -- muss vor jeder Weiterfuehrung geschehen. Dann: **Voisin-Beispiele numerisch pruefen** (GHR-Negativitaet) -- geringer Aufwand, klares Ergebnis.

### Schluessel-Literatur
- Bhatt-Scholze-Morrow-Tsuji (2019-2026), Binyamini-Schmidt-Thomas, Voisin

---

## 4. NS LogDistance (9.0/10)

### Kritische Luecke
- **TLL + LDI:** Nicht verifiziert fuer konkretes PDE-System (Framework ohne Proof-of-Life)

### Loesungsansaetze (konsolidiert aus allen Reviews)

| # | Ansatz | Quelle | Prioritaet | Machbarkeit |
|---|--------|--------|:-----------:|:-----------:|
| 1 | Reaktions-Diffusion als Pilotverifikation (2D, glatte Nichtlinearitaet) | Copilot | **Sofort** | **Hoch** |
| 2 | Hyperviskose NSE als Brueckensystem (beta >= 5/4, Inertialmannigfaltigkeiten bewiesen) | Copilot | Hoch | Mittel |
| 3 | Gedaempfte Wellengleichung | Copilot | Mittel | Mittel |
| 4 | Ginzburg-Landau + Liftung auf 2D NS | Gemini | Hoch | Mittel |
| 5 | Kuramoto-Sivashinsky (1D, deterministisches Chaos, keine Singularitaeten) | Gemini / Deep Research | **Sofort** | **Hoch** |
| 6 | Stochastische Stabilisierung: Feller-Rauschen erzwingt TLL "fast sicher" | Gemini | Mittel | Mittel |
| 7 | DFC aus Turbulenz als topologischer Trichter fuer TLL | Deep Research | Mittel | Mittel |
| 8 | Kanonisches Auswahlprinzip: "minimal metric derivative selection" oder a.e.-Eindeutigkeit | Copilot Einzel | Hoch | Mittel |
| 9 | TLL als Aubin-property / metric regularity der Projektionskorrespondenz Gamma | Copilot Einzel | Hoch | Mittel |
| 10 | Bridge Lemma: squeezing + determining modes => log-Modulus => TLL | Copilot Einzel | **Sofort** | **Mittel** |
| 11 | Doubling/Assouad: approximate tangent cone => TLL mit log(1/r)-Lipschitz | Copilot Einzel | Mittel | Mittel |
| 12 | LDI als Tail-Integral-Aequivalenz: int_0^R (1/s)(int_{0<d<=s} |u'|) ds < inf | Copilot Einzel | Hoch | Hoch |
| 13 | STC-Tracking-Lemma: d(t) <= e^{-lambda t} => log(1/d) in L^p | Copilot Einzel | Mittel | Mittel |
| 14 | Proposition 7.1 praezisieren: BV -> L^inf Einbettung + d(t)>0 a.e. | Copilot Einzel | Mittel | Hoch |
| 15 | Toy model: endlichdimensionales System mit fraktalem Attraktor (TLL+LDI komplett) | Copilot Einzel | **Sofort** | **Hoch** |
| 16 | "High-mode stability conjecture" als hinreichende Bedingung fuer TLL bei NS | Copilot Einzel | Hoch | Mittel |
| 17 | Variational analysis / metric regularity als Literatur-Achse einbauen | Copilot Einzel | Mittel | Hoch |

### Empfohlener naechster Schritt
**TLL+LDI auf Kuramoto-Sivashinsky oder Reaktions-Diffusion beweisen** -- "Proof of Life" fuer das Framework. Parallel: **Toy model** (endlichdimensional) als Minimal-Beispiel.

---

## 5. NS Skeleton (7.5/10)

### Kritische Luecken
- **AGC1:** Positiver Reach des Attraktors (unbewiesen)
- **Condition (D):** Moeglicherweise vakuos
- **Gronwall-Luecke:** Zeitabhaengigkeit des Attractor-Minimierers (als Assumption G gekennzeichnet)
- **Zirkularitaet:** Regularitaet <-> Attraktor (Standardliteratur setzt Glaette voraus)

### Loesungsansaetze (konsolidiert aus allen Reviews)

| # | Ansatz | Quelle | Prioritaet | Machbarkeit |
|---|--------|--------|:-----------:|:-----------:|
| 1 | Reach-Lemma: Rueckwaerts-Stabilitaet + Injektivitaet endlicher Projektion => positiver Reach | Copilot | Hoch | Mittel |
| 2 | Whitney-Stratifizierung statt glatte Mannigfaltigkeit | Copilot | Mittel | Mittel |
| 3 | No-cusp via Determining Modes (bi-Lipschitz Projektion) | Copilot | Hoch | Mittel |
| 4 | Multiplikative Oseledets-Kegel (invariante Konvexitaet im Tangentialraum) | Gemini | Mittel | Niedrig |
| 5 | Fraktionale Einbettung H^{1+eps}: Log-Lipschitz genuegt fuer glatte Normalenbuendel | Gemini | Mittel | Mittel |
| 6 | Defekte Masstheorie: Dissipation als Radon-Mass, Condition (D) = Positivitaet des singulaeren Teils | Gemini | Mittel | Mittel |
| 7 | Kompensierte Kompaktheit (Tartar): Oszillationen in Vortizitaet => integrale Divergenz | Gemini | Mittel | Niedrig |
| 8 | Ginzburg-Landau als TLL/LDI-Bruecke + Liftung auf 2D NS | Gemini | Hoch | Mittel |
| 9 | Approximate Prox-Regularity: lokaler Reach nach Entfernung magerer Menge | Deep Research | Hoch | Mittel |
| 10 | Schwache Auswahlfunktionen beschraenkter Hausdorff-Variation statt Lipschitz-Projektionen | Deep Research | Hoch | Mittel |
| 11 | Geometrisches Euler-System: Galerkin-Hierarchie als funktionales Reach-Surrogat | Deep Research | Mittel | Niedrig |
| 12 | H^1-Lift fuer Condition (D): Enstrophie divergiert unter Blow-up (nicht vakuos!) | Deep Research | **Sofort** | **Mittel** |
| 13 | Vortex-Stretching skaliert hoechstens C||Delta u||^{2-delta} => Viskositaet dominiert | Deep Research | Hoch | Mittel |
| 14 | Selection Regularity im trajectory-attractor Setting (Chepyzhov-Vishik) | Copilot Einzel | **Sofort** | **Mittel** |
| 15 | Ekeland-Variationsprinzip als Projektions-Ersatz ("almost minimisers") | Copilot Einzel | **Sofort** | **Mittel** |
| 16 | Condition D als Rate-Klasse: welche alpha machen D wahr/falsch? | Copilot Einzel | Hoch | Mittel |
| 17 | limsup-Dominanz auf Sequenz t_n statt volle Dominanz | Copilot Einzel | Hoch | Hoch |
| 18 | Set-valued BV-Selection Mini-Programm | Copilot Einzel | Mittel | Mittel |
| 19 | Quantitative Tracking-Lemmas: Massabschaetzung fuer {t: d(t) <= eps} | Copilot Einzel | Mittel | Mittel |
| 20 | Gronwall-Luecke via Lipschitz-Regularitaet (Barbu/Cannone Log-Lipschitz) | Ext. Eval. | Hoch | Mittel |
| 21 | Log-Dirichlet-Quotienten in KL-Joint-Minimierung integrieren | Ext. Eval. | Hoch | Mittel |
| 22 | Zirkularitaet durchbrechen: strikte Kontraktion von Q_N unter KL-Metrik | Ext. Eval. | Hoch | Mittel |
| 23 | Zelik (2022) Trajektorien-Attraktoren einarbeiten | Lit-Review | Hoch | Hoch |
| 24 | MFG-Wasserstein-Formulierung vertiefen (U Michigan 2026) | Lit-Review | Mittel | Mittel |

### Empfohlener naechster Schritt
**H^1-Lift fuer Condition (D)** (Enstrophie hat keine obere Schranke) + **Ekeland-Ersatz** fuer Projektion + **Chepyzhov-Vishik** fuer Selection Regularity. Drei parallele Angriffe auf die drei Hauptluecken.

### Schluessel-Literatur
- Zelik (2022), Barbu/Cannone (2016), Ilyin-Kalantarov-Zelik (2025), Lytchak (2023), Mallet-Paret/Sell (1988)

---

## 6. P vs NP (8.3/10)

### Kritische Luecke
- **Uniformitaets-Bruecke:** Uebergang von existenzieller zu uniformer Schranke ungeloest
- **Review-Finding PNP-1:** Witness Entropy Gap -- Luecke bei vielen Witnesses
- **Review-Finding PNP-2:** Barrier-Immunitaet gilt nur fuer unbeschraenktes K, nicht K^t

### Loesungsansaetze (konsolidiert aus allen Reviews)

| # | Ansatz | Quelle | Prioritaet | Machbarkeit |
|---|--------|--------|:-----------:|:-----------:|
| 1 | Compression-to-circuit Transfer-Lemma (AC^0 oder monotone circuits) | Copilot | Hoch | Mittel |
| 2 | Min-Entropy statt worst-case K^t (samplable distributions => average-case) | Copilot | Hoch | Mittel |
| 3 | Adversariale Zeugenfamilie (globaler Konsistenzcheck, nicht lokal extrahierbar) | Copilot | Mittel | Niedrig |
| 4 | Spieltheoretisches Minimax: Entropie-Maximierer vs Schaltkreis-Minimierer | Gemini | Hoch | Niedrig |
| 5 | Quanten-Zertifikate: ESC in MIP* (Verschraenkung -> lokale Greifbarkeit) | Gemini | Niedrig | Niedrig |
| 6 | Kolyvagin-Hierarchie (BSD-Analogie fuer "Raenge" der Ineffizienz) | Gemini | Mittel | Niedrig |
| 7 | Selbstreferentielle Diagonalisierung (Instanz beschreibt TM A) | Deep Research | Hoch | Niedrig |
| 8 | Extended Frege Refutationen + hohe K(pi|phi) => kein uniformer Suchalgorithmus | Deep Research | Hoch | Mittel |
| 9 | MCSP nicht in P/poly als minimaler Durchbruch (impliziert Entropy Hardness) | Deep Research | Hoch | Niedrig |
| 10 | Bridge Target 1: Instance compression => K^t lower bounds | Copilot Einzel | **Sofort** | **Mittel** |
| 11 | Bridge Target 2: Proof complexity -- Entropie ueber proofs statt witnesses | Copilot Einzel | **Sofort** | **Mittel** |
| 12 | Bridge Target 3: PRG-basierte Kandidatenfamilien (pseudorandom witness-Struktur) | Copilot Einzel | Hoch | Mittel |
| 13 | "Uniformity Bridge" Sektion als Lemma-Target-Liste schreiben | Copilot Einzel | **Sofort** | **Hoch** |
| 14 | Levins Universeller Suchalgorithmus als Uniformity-Collapse | Gemini | Mittel | Mittel |
| 15 | Valiant-Vazirani + K^t Few-Witness Extension | Copilot | Mittel | Mittel |
| 16 | Resolution-Breite ~ K(pi|phi) formal beweisen | Copilot | Mittel | Mittel |
| 17 | MCSP-Internalisierung (Allender => EH => P!=NP) | Copilot | Mittel | Niedrig |

### Empfohlener naechster Schritt
**Uniformity Bridge Sektion schreiben** (Lemma-Targets mit Statement-Form) -- verwandelt Re-Kodierung in konkretes Forschungsprogramm. Parallel: K vs K^t klar trennen (PNP-2 adressieren).

### Schluessel-Literatur
- Allender-Hirahara-Koucky (2018-2024), Li-Vitanyi

---

## 7. Turbulenz (8.5/10 -- "Clean Win", Prioritaet 1)

### Offene Punkte
- **NL/DFC:** Assumptions NL (entropy-compatible cascade) und ND nicht rigoros bewiesen
- **Spektrale Poincare-Ungleichung:** D_F^nu <-> epsilon_nu Verknuepfung fehlt
- **Intrinsische E*-Charakterisierung:** Ohne a priori K41

### Loesungsansaetze (konsolidiert aus allen Reviews)

| # | Ansatz | Quelle | Prioritaet | Machbarkeit |
|---|--------|--------|:-----------:|:-----------:|
| 1 | Referenz-Robustheit: K41 ist einzige selbstkonsistente Referenz (Fixed-Point) | Copilot Einzel | **Sofort** | **Mittel** |
| 2 | Minimax-Formulierung (Closure vs Energieprofil, K41 = Saddle-Point) | Copilot Einzel | Hoch | Mittel |
| 3 | DFC-with-slack: summierbare Fehler eta_j, delta_j | Copilot Einzel | **Sofort** | **Hoch** |
| 4 | "How to falsify" Sektion: Plots, Shell-Definitionen, Fensterbreiten fuer DNS | Copilot Einzel | **Sofort** | **Hoch** |
| 5 | Spektrale Poincare: D_F^nu <-> eps_nu Verknuepfung suchen | Copilot | Hoch | Mittel |
| 6 | Intrinsisches E* aus KHM-Relation (ohne a priori K41) | Copilot | Hoch | Mittel |
| 7 | Intermittenz-Exponenten aus Hessian-Fluktuationen rigoros ableiten | Copilot | Mittel | Mittel |
| 8 | Figuren erstellen (Beweiskette-Schema, F[E]-Landschaft) | Copilot | **Sofort** | **Hoch** |
| 9 | Energiekaskade als Minimax via Inertial Manifolds (Eden/Foias) | Ext. Eval. | Hoch | Mittel |
| 10 | Vortex-Tubes als knotted Quasipartikel (topologischer Strafterm) | Knoten-Input | Mittel | Niedrig |
| 11 | Entropy-Kompatibilitaet als "downhill in knot-space" (Remark in Discussion) | Knoten-Input | Niedrig | Mittel |

### Empfohlener naechster Schritt
**DFC-with-slack + How-to-falsify Sektion + Figuren** -- macht das Paper DNS-testbar und Phys. Rev. E ready. Das ist der naechste Schritt vor Einreichung.

---

## 8. Yang-Mills (7.5/10)

### Kritische Luecken
- **Kontinuumslimes (beta -> inf):** Holley-Stroock-Konstanten divergieren exponentiell
- **CL2:** Uniforme LSI in physikalischer Distanzskala unbewiesen
- **RP im Weak-Coupling:** Reflexionspositivitaet nur fuer Strong-Coupling bewiesen

### Loesungsansaetze (konsolidiert aus allen Reviews)

| # | Ansatz | Quelle | Prioritaet | Machbarkeit |
|---|--------|--------|:-----------:|:-----------:|
| 1 | FST-Minimierung => uniform konvexe LSI/PI (beta-uniforme Schranken) | Copilot | Hoch | Mittel |
| 2 | RG-Monotonie-Kandidat M(beta): bei beta=0 bekannt, Monotonie zeigen | Copilot | Hoch | Mittel |
| 3 | Zwei-Phasen-Interpolation: analytische Fortsetzbarkeit als Phasenuebergangs-Ausschluss | Copilot | Mittel | Niedrig |
| 4 | Kramers-Wannier-Dualitaet: schwache Kopplung -> stark gekoppeltes Monopol-Modell | Gemini | Hoch | Niedrig |
| 5 | Haar-Mass-Entropie: Confinement rein entropisch (Phasenraumentropie unterdrueckt) | Gemini | Mittel | Niedrig |
| 6 | Kirk Pro-Polymer auf SU(3) erweitern (Peter-Weyl aller irreduziblen Darstellungen) | Deep Research | Hoch | Mittel |
| 7 | Massenluecke als RG-Fixpunkt-Problem (IR-Fixpunkt mit positiver Kruemmung) | Deep Research | Hoch | Mittel |
| 8 | Zero-Free Strip fuer 4D SU(3): verbietet Phasenuebergaenge 1. Ordnung | Deep Research | Hoch | Niedrig |
| 9 | Gap-Transport-Lemma: PI-Konstante K(a) bleibt unter Block-RG stabil | Copilot Einzel | **Sofort** | **Mittel** |
| 10 | Effective influence matrix (Dobrushin-Style): Spektralradius < 1 => Gap | Copilot Einzel | **Sofort** | **Mittel** |
| 11 | "Hessian = Var" durch PI/LSI-Konstante ersetzen als primaeres Gap-Objekt | Copilot Einzel | **Sofort** | **Hoch** |
| 12 | Constants & Norms Anhang (Normierung, Metrik, Laplace-Beltrami-Gap, Dirichlet-Form) | Copilot Einzel | **Sofort** | **Hoch** |
| 13 | Topologische Sektoren: scharfe Definition oder in Outlook verschieben | Copilot Einzel | Mittel | Hoch |
| 14 | CL2 ersetzen durch finite-volume bootstrap (strukturell schwaecher) | Copilot Einzel | Hoch | Mittel |
| 15 | RFEP-Paper: Claim-Level trennen (Level 0/1/2) | Copilot Einzel | **Sofort** | **Hoch** |
| 16 | Typo "info(M) =4>0" in Einleitung entfernen | Copilot Einzel | **Sofort** | **Hoch** |
| 17 | Kontinuumslimes: hierarchische Block-Spin-RG INNERHALB des LSI-Frameworks | Ext. Eval. | Hoch | Mittel |
| 18 | Reflexionspositivitaet auf Weak-Coupling erweitern | Ext. Eval. | Hoch | Mittel |
| 19 | LSI-Konstanten unter Skalentransformation binden (Birkhoff-projektive Kontraktion) | Ext. Eval. | Hoch | Mittel |
| 20 | Kirk-Methodik auf FST-Beweis abbilden (Detailvergleich Pro-Polymer-Abschaetzungen) | Lit-Review | Hoch | Mittel |
| 21 | Vergleichstabelle FST vs. Kirk vs. TMT | Lit-Review | Mittel | Hoch |
| 22 | Kirk's zero-free strip fuer RG-Fluss nutzen | Lit-Review | Hoch | Mittel |
| 23 | Mixture/Metastable Decomposition statt Holley-Stroock (Konsens 4/4 Reviews) | Cross-Review | Hoch | Mittel |
| 24 | Kirk-Extension SU(3): Character Ratio Bounds + Heat-Kernel Domination | Cross-Review | Hoch | Mittel |

### Empfohlener naechster Schritt
**PI/LSI als primaeres Objekt** (statt Hessian=Var) + **Gap-Transport-Lemma** formulieren + **Constants Anhang** schreiben. Dann: Kirk-Detailvergleich fuer CL1-CL3.

### Schluessel-Literatur
- Kirk (2026, Zenodo), Garcia Baquero (2026), Zegarlinski (1992), "Emergence as Regime Formation" (2025)

---

## 9. Framework / FST-RH (9.4/10)

### Offene Punkte
- FST-RH Parts I-III unveroeffentlicht
- Numerische Zertifikate (29/30 abgeschlossen, lambda=640000 offen)
- Claim-Level-Trennung fehlt

### Loesungsansaetze (konsolidiert aus allen Reviews)

| # | Ansatz | Quelle | Prioritaet | Machbarkeit |
|---|--------|--------|:-----------:|:-----------:|
| 1 | Minimal Verifiable Core (MVC): nur Axiome, Schluessel-Lemmas, Proof-Graph | Copilot | **Sofort** | **Hoch** |
| 2 | Verification Appendix: jede Abschaetzung als Checkable Claim | Copilot | Hoch | Hoch |
| 3 | Independence Note: Framework bleibt korrekt auch wenn RH-Parts falsch sind | Copilot | **Sofort** | **Hoch** |
| 4 | Lean4/Coq Formalisierung (Teil-Lemma, z.B. Shift Parity Lemma) | Gemini | Hoch | Mittel |
| 5 | Part I als bedingter Beweis auf Zenodo: "Unter Annahme AP3 folgt..." | Gemini | **Sofort** | **Hoch** |
| 6 | Entropisches RH-Argument (Nullstellen => Kompressions-Algorithmus der FST-Schranken verletzt) | Gemini | Mittel | Niedrig |
| 7 | Resolventen-Operator auf diskretem Gitter + Lean4/Coq-Verifikation | Deep Research | Hoch | Mittel |
| 8 | L-Funktions-Nullstellen als Eigenwerte eines Hamiltonoperators | Deep Research | Mittel | Niedrig |
| 9 | Arithmetisches Euler-System ueber Spektralfluss | Deep Research | Mittel | Niedrig |
| 10 | Shift Parity Lemma in Lean4/Coq formalisieren | Deep Research | Hoch | Mittel |
| 11 | Claim-Level formal trennen (Level 0 Framework / Level 1 RH / Level 2 Domain) | Copilot Einzel | **Sofort** | **Hoch** |
| 12 | Universelles Proof-Schema formalisieren (6 Slots) als eigene Note | Copilot Einzel | Hoch | Hoch |
| 13 | "Second-Order Resolvent Dominance" als eigenstaendiges Resultat dokumentieren | Deep Research / Cross | Hoch | Hoch |

### Empfohlener naechster Schritt
**Part I als bedingter Beweis auf Zenodo** + **MVC + Independence Note** + **Claim-Level trennen**. Das macht das Framework fuer die Community zugaenglich, ohne ueberzogene Ansprueche.

---

## 10. Cross-Paper / Uebergreifend

| # | Ansatz | Quelle | Prioritaet | Machbarkeit |
|---|--------|--------|:-----------:|:-----------:|
| 1 | Arakelov-Bruecke: Hodge HR-Form + BSD Neron-Tate als Spezialfaelle | Cross-Review | Mittel | Niedrig |
| 2 | Meta-Theorem: Pattern A Master Stability Theorem formal formulieren | Cross-Review | Hoch | Mittel |
| 3 | YM->NS Transfer: Poincare-Dobrushin auf Galerkin-trunkierte NS | Cross-Review | Mittel | Mittel |
| 4 | SOS-Zertifikat als universelles Werkzeug (Convexity Trap + Gronwall + Kontinuumslimes) | Ext. Eval. | Hoch | Niedrig |
| 5 | Universelles Proof-Schema als Note/Framework-Kapitel | Copilot Einzel | Hoch | Hoch |
| 6 | Alle internen Referenzen (FST-Unified, FST-RH3) mit Zenodo-DOIs versehen | Review-Finding | **Sofort** | **Hoch** |
| 7 | Resolvent-Remarks: nur eigene Zeile + Kontext statt volle 5x4-Tabelle (Anti-Self-Plagiarism) | Review-Finding | Hoch | Hoch |

---

## 11. Review-Findings: Korrekturen vor Publikation

Diese Punkte betreffen keine neuen Forschung, sondern Korrekturen bestehender Papers.

### Hodge (KRITISCH)
- [ ] **HOD-5:** AP1 Vorzeichen: (-1)^p Q_L >= 0 statt Q_L >= 0 (Def 3.1, Eq 6, alle Folge-Propositionen)
- [ ] **HOD-1:** Proposition 2.3 Beweis umformulieren (zwei separate Annahmen klar machen)
- [ ] **HOD-4:** Theorem 5.1 als "Proof outline" oder "Strategy" labeln
- [ ] **HOD-2:** Lefschetz-Operator Notation klaeren (Fussnote oder Umbenennung)
- [ ] **HOD-3:** 7 nicht-zitierte Bibliographie-Eintraege mit \cite versehen
- [ ] **HOD-6:** Functoriality Proofs in Sec 5 praezisieren

### BSD
- [ ] **BSD-1:** Theorem 4.1 als "Conditional/Framework Theorem" labeln
- [ ] **BSD-2:** Step 2 praezise formulieren
- [ ] **BSD-3:** E2 als "under assumption of leading-term formula" umformulieren
- [ ] **BSD-4:** CRM-Referenz durch Zenodo-DOI ersetzen
- [ ] **BSD-5:** Iwasawa-Anwendung praezisieren (welche Faelle genau)

### P vs NP
- [ ] **PNP-2:** K vs K^t klar trennen (Barrier-Immunitaet nur fuer K)
- [ ] **PNP-3:** Abstract deutlicher als Reformulierung kennzeichnen
- [ ] **PNP-4:** Theorem 4.2 Proof Sketch praezisieren
- [ ] **PNP-5:** AP4 differenzierter formulieren (nur random k-SAT nahe Schwelle)

### Alle alten Papers
- [ ] Interne Referenzen durch Zenodo-DOIs ersetzen
- [ ] Resolvent-Remarks auf eigene Zeile kuerzen

---

## Prioritaeten-Matrix (Impact x Machbarkeit)

### SOFORT machbar + hoher Impact (Quick Wins) -- Sprint 4

| Prio | Aufgabe | Paper | Aufwand |
|:----:|---------|-------|:-------:|
| 1 | AP1 Vorzeichen-Fehler korrigieren | Hodge | 1-2h |
| 2 | DFC-with-slack + How-to-falsify + Figuren | Turbulenz | 1-2d |
| 3 | PI/LSI als primaeres Objekt + Constants Anhang | Yang-Mills | 1-2d |
| 4 | Uniformity Bridge Lemma-Targets schreiben | P vs NP | 1d |
| 5 | FST-RH Part I auf Zenodo + Independence Note | Framework | 1-2d |
| 6 | MVC + Claim-Level trennen | Framework | 1d |
| 7 | Numerische Experimente BSD Rang 2-4 | BSD | 2-3d |
| 8 | Voisin-Beispiele GHR-Negativitaet pruefen | Hodge | 1d |
| 9 | TLL+LDI Toy model (endlichdimensional) | NS LogDist | 2-3d |
| 10 | TLL+LDI auf Kuramoto-Sivashinsky | NS LogDist | 3-5d |
| 11 | Review-Korrekturen (BSD-1, PNP-2, Hodge, Referenzen) | Alle | 2-3d |

### Mittelfristig (1-2 Jahre, hoher Hebel) -- Phase 2

| Prio | Aufgabe | Paper |
|:----:|---------|-------|
| 12 | Zwei-Skalen gamma(R) fuer Screening | Dark Energy |
| 13 | Semiklassische Effektivwirkung | Dark Energy |
| 14 | Gap-Transport-Lemma entlang RG | Yang-Mills |
| 15 | Bridge Lemma squeezing => TLL | NS LogDist |
| 16 | H^1-Lift fuer Condition D | NS Skeleton |
| 17 | Ekeland-Variationsprinzip als Projektions-Ersatz | NS Skeleton |
| 18 | Kirk-Detailvergleich fuer CL1-CL3 | Yang-Mills |
| 19 | Mumford-Tate-Analyse fuer CM-Typen | Hodge |
| 20 | Instance compression => K^t + Proof complexity | P vs NP |
| 21 | Lean4/Coq Teil-Formalisierung | Framework |
| 22 | Pattern A Master Stability Theorem | Cross-Paper |
| 23 | Euclid/DESI/LISA Forecasts fuer eta=5/4 | Dark Energy |
| 24 | Mixture/Metastable Decomposition (Konsens 4/4) | Yang-Mills |

### Langfristig (Forschungsprogramm, 2+ Jahre) -- Phase 3

| Prio | Aufgabe | Paper |
|:----:|---------|-------|
| 25 | Higher Gross-Zagier fuer Rang >= 2 | BSD |
| 26 | O-minimale Strukturen fuer Sha | BSD |
| 27 | Prismatische Positivitaet AP' | Hodge |
| 28 | Kontinuumslimes alle beta | Yang-Mills |
| 29 | Uniformitaets-Beweis P != NP | P vs NP |
| 30 | SOS-Zertifikat auf Hilbertraeumen | Cross-Paper |
| 31 | Arakelov-motivische Volumenformel | BSD |
| 32 | Thermodynamische Holographie | Dark Energy |

---

## Publikationsreihenfolge (konsolidiert aus allen Reviews)

| Rang | Paper | Score | Status | Naechster Meilenstein |
|:----:|-------|:-----:|--------|----------------------|
| 1 | Turbulenz | 8.5 | "Clean Win" | DFC-slack + Figuren, dann Phys. Rev. E |
| 2 | Navier-Stokes (Skeleton) | 7.5 | Conditional Framework | H^1-Lift + Ekeland, dann Zenodo |
| 3 | Hodge | 7.0 | AP1-Fix noetig | Vorzeichen-Korrektur + Voisin-Test |
| 4 | BSD | 7.5 | Reformulierung | Conditional labeln + Numerik Rang 2-4 |
| 5 | Yang-Mills | 7.5 | Reframing noetig | PI/LSI + Kirk-Mapping + Constants |
| 6 | P vs NP | 8.3 | Ueberzeugend, Luecke gross | Uniformity Bridge schreiben |
| 7 | Dark Energy | 9.0 | Framework Note | Zwei-Skalen gamma + Effektivwirkung |
| 8 | NS LogDistance | 9.0 | Proof-of-Life fehlt | KS/RD Verifikation |
| 9 | Framework/FST-RH | 9.4 | Unveroeffentlicht | Part I auf Zenodo + MVC + IndNote |

> Hinweis: Hoehere Scores bedeuten NICHT hoehere Publikationsbereitschaft -- Score misst den Reifegrad des Beweises, nicht die Publikationsnaehe. Turbulenz (8.5) steht auf Rang 1, weil es am geschlossensten ist.

---

## Meta-Einsicht: "Second-Order Resolvent Dominance"

Alle drei LLM-Reviews konvergieren auf ein universelles Metamuster:
In jeder FST-Domaene ist die Loesung auf dem **linearen** Niveau entartet oder vakuos.
Der Durchbruch erfordert zwingend eine **Kopplung zweiter Ordnung**:

| Domaene | 0. Ordnung (neutral) | 1. Ordnung (degeneriert) | 2. Ordnung (entscheidet) |
|---------|----------------------|--------------------------|--------------------------|
| BSD | Funktionalgleichung L(s)=L(1-s) | L'(E,1)=0 bei hohem Rang | L^(r) + Regulator + Sha |
| Dark Energy | M_Pl^4 Vakuumenergie | Loop-Korrekturen (divergent) | Topologische Selektion rho_CC |
| NS | Energieerhaltung | Enstrophie-Dissipation (monoton) | Attraktor-Reach + Projektion |
| P vs NP | Existenz von Witnesses | SAT/UNSAT (NP-Entscheidung) | Uniforme Entropie-Schranke |
| Yang-Mills | Freie Gibbs-Masse | Poincare-Konstante (skalen-abhaengig) | Skalenuebergreifender Gap (RG) |

Dieses Muster sollte als eigenstaendiges Resultat im Framework-Paper dokumentiert werden (TODO bereits in Section 10 vermerkt).

---

## Literatur-Abgleich (ausstehend)

Aus `LITERATUR_DeepResearch_2026-03-16.md` -- gegen Paper-Bibliographien pruefen:

- [ ] Kim (2022-2024) -> BSD Paper
- [ ] Burns/Sano/Sakamoto -> BSD Paper
- [ ] Burungale/Tian (2026) -> BSD Paper
- [ ] Wetterich (2024) -> Dark Energy Paper
- [ ] Vasak/Garcia Baquero (2026) -> Dark Energy Paper
- [ ] Bhatt/Scholze/Morrow/Tsuji -> Hodge Paper
- [ ] Binyamini/Schmidt/Thomas -> Hodge Paper
- [ ] Ilyin/Kalantarov/Zelik (2025) -> NS Papers
- [ ] Lytchak (2023) -> NS LogDistance Paper
- [ ] Allender/Hirahara/Koucky -> P-vs-NP Paper
- [ ] Kirk (2026) -> Yang-Mills Paper (bereits referenziert)

---

*Generiert am 2026-03-16 aus 12 Review-Dokumenten. Naechstes Update nach Sprint 4.*
