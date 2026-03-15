# BEWEISNOTIZ -- Yang-Mills Massenluecke
# Stand: 2026-03-15 (nach Sprint 3 Review)
# Status: Strong-Coupling bewiesen, Kontinuumslimes OFFEN
# Review: 6C/9M/9L Issues identifiziert und korrigiert (Sprint 3)

===============================================================================
## Problemstellung
===============================================================================

**Millennium-Problem (Jaffe-Witten 2000):** Fuer jede kompakte einfache
Eichgruppe G existiert die Quantenfeldtheorie auf R^4 (Wightman-Axiome)
und besitzt eine Massenluecke Delta > 0:

  inf sigma(M)|_{H perp Omega} = Delta > 0

wobei M = sqrt(P_mu P^mu) der Massenoperator auf dem physikalischen
Hilbertraum ist.

**Ansatz (FST-YM):** Variationelle Free-Energy-Methode. Konstruktion eines
renormierten Free-Energy-Funktionals F[mu] ueber Gitter-Eichfeldmasse mu.
Strikte Konvexitaet von F liefert Spektralluecke im Transferoperator.


===============================================================================
## Beweiskette
===============================================================================

### Schritt 1: Wilson-Gitter-Regularisierung
**Status: BEWIESEN (Standard)**

Wilson-Wirkung auf Hyperkubus-Gitter Lambda = (aZ)^4 cap [-L,L]^4:
  S_W[U] = beta * Sigma_p (1 - Re tr U_p / N)
mit beta = 2N/g^2, Plaquettenvariablen U_p.
Gittermass mu_Lambda = Z^{-1} exp(-S_W) prod dU_ell (Haar-Mass).

Free-Energy-Funktional:
  F[mu] = <S_W>_mu + T * D_KL(mu || lambda)
mit KL-Divergenz. Eindeutiger Minimierer bei T=1: mu_Lambda.


### Schritt 2: Strikte Konvexitaet von F
**Status: BEWIESEN (Lemma 2.3)**

F[mu] ist strikt konvex auf P_ac(C_Lambda):
  F[theta mu_0 + (1-theta) mu_1] < theta F[mu_0] + (1-theta) F[mu_1]

Beweis: Energie-Term linear, Entropie-Term (KL) strikt konvex
(x log x strikt konvex auf (0,inf)). Summe strikt konvex.

ANMERKUNG: Konvexitaet allein gibt KEINEN Mass Gap -- jedes Gibbs-Mass
hat das. Der Motor ist Poincare/Dobrushin (Schritte 3-4).


### Schritt 3: Single-Link Poincare-Ungleichung via Holley-Stroock
**Status: BEWIESEN (Theorem 4.1, Step 1)**

Fuer gauge-invariante f in L^2(mu_Lambda) mit <f> = 0:
  Var_mu(f) <= (1/kappa) * E(f,f)
mit Dirichlet-Form E des Heat-Bath-Dynamics.

Bedingte Dichte am Link ell hat beschraenkte Oszillation:
  osc(W) <= 2 * d_ell * beta  (d_ell = 6 in d=4)

Holley-Stroock Bounded Perturbation Lemma liefert:
  kappa_link(beta, G) >= kappa_Haar(G) * exp(-2 * d_ell * beta)

UNIFORM in Randkonfiguration und Gittervolumen |Lambda|.


### Schritt 4: Volumenunabhaengigkeit via Dobrushin-Zegarlinski
**Status: BEWIESEN fuer beta < beta_0 (Strong-Coupling)**

Dobrushin-Einflusskonstanten c_{ell,ell'} <= C(G) * beta (Lipschitz).
Dobrushin-Konstante: eta(beta) = sup_ell Sigma_{ell'} c_{ell,ell'} <= C(d,G) * beta.

Fuer beta < beta_0 = 1/C(d,G):  eta < 1, und Zegarlinski liefert:
  kappa(Lambda) >= (1 - eta) * kappa_link =: kappa_* > 0

UNABHAENGIG von |Lambda|!

**EINSCHRAENKUNG:** Gilt NUR im Strong-Coupling-Regime. Fuer beta >= beta_0
(Weak-Coupling, physikalisch relevant) bricht das Argument zusammen, weil
die Holley-Stroock-Konstanten exponentiell mit beta divergieren.


### Schritt 5: Transferoperator-Spektralluecke
**Status: BEWIESEN (Theorem 4.1, Step 3)**

Transferkern K(U,W) durch Integration ueber temporale Links.
Symmetrisch (Zeitumkehr), positiv (Reflexionspositivitaet, Lemma 3.3).
Normierter Markov-Operator P ist reversibel bzgl. Slice-Marginal pi.

Poincare-Ungleichung uebersetzt zu:
  lambda_1(P) <= 1 - kappa

Exponentielles Clustering:
  |<O(t) O(0)>_conn| <= (1-kappa)^t * ||f||^2

Gitter-Massenluecke:
  Delta_Lambda = -log(1-kappa) >= kappa > 0


### Schritt 6: Reflexionspositivitaet
**Status: BEWIESEN (Lemma 3.3)**

3-Schritt-Beweis:
1. Reduktion auf Crossing-Kern-Positivitaet (fuer festes U_0)
2. Einzelne Plaquette: Peter-Weyl + nichtnegative Charakterkoeffizienten
   a_pi >= 0 (Osterwalder-Seiler 1978)
3. Produkt positiver Kerne: Schur-Produktsatz

Gilt fuer Wilson-Wirkung und Heat-Kernel-Wirkung.


### Schritt 7: Kontinuumslimes (CONDITIONAL)
**Status: OFFEN -- Kernproblem des Millennium-Problems**

Theorem 4.2 formuliert drei Annahmen:
  (CL1) OS-positiver Kontinuumslimes existiert (Schwinger-Funktionen)
  (CL2) Uniformes exponentielles Clustering in physikalischer Distanz
  (CL3) Uniforme Normierung der renormierten Observablen

Unter (CL1)-(CL3): Rekonstruiertes QFT hat Massenluecke Delta >= Delta_0.

**Problem:** Im Kontinuumslimes a -> 0 gilt beta(a) -> inf (asymptotische
Freiheit), d.h. das Strong-Coupling-Regime wird verlassen. Die Gitter-
Massenluecke Delta_lat >= kappa_* > 0 (Gittereinheiten) schrumpft,
waehrend Delta_phys = Delta_lat/a positiv bleiben soll.


===============================================================================
## Offene Luecken
===============================================================================

### L1: Kontinuumslimes (KRITISCH)
Strong-Coupling-Gap resultiert aus Gitter-Diskretheit. Im Limes beta -> inf
divergieren Holley-Stroock-Konstanten exponentiell:
  kappa_link ~ exp(-12 * beta)
Die Poincare-Ungleichung wird wertlos.

**Moegliche Loesung:** Hierarchische Block-Spin-RG innerhalb des LSI-Frameworks.
Skalenaufgeloeste Zerlegung statt globales Holley-Stroock.
Referenz: "Emergence as Regime Formation" (2025) -- polynomielle untere Schranken
via Dobrushin-Shlosman + LSI.

### L2: Extension auf alle beta
Dobrushin-Zegarlinski gilt nur fuer beta < beta_0. Persistenz eines uniformen
Gaps fuer alle beta > 0 im thermodynamischen Limes ist offene Vermutung.
Finite-Volume-Analytizitaet schliesst scharfe Phasenuebergaenge bei fixem
|Lambda| aus, aber nicht im Limes |Lambda| -> inf.

### L3: Kirk-Mapping (conditional)
Kirk (2026) hat SU(2)-Kontinuumslimes via Dobrushin-Shlosman Complete Analyticity
und Pro-Polymer Cluster Expansion konstruiert:
- CL1: Anisotropie O(a^2), volle O(4)-Restauration
- CL2: Zero-Free Strip -> uniformes Clustering
- CL3: Pro-Polymer Summability Bounds
Wenn auf allgemeine kompakte G erweiterbar, waeren (CL1)-(CL3) substantiiert.

### L4: Reflexionspositivitaet im Weak-Coupling
Lemma 3.3 gilt fuer alle beta, aber die NUTZUNG im Beweis (Schritt 5) erfordert
die Poincare-Konstante aus Schritt 4, die nur fuer beta < beta_0 gilt.


===============================================================================
## Naechste Schritte
===============================================================================

1. **LSI-Konstanten unter RG-Fluss binden:** Beweisen dass Holley-Stroock-Konstanten
   unter Skalentransformation invariant oder polynomiell beschraenkt bleiben
   (Birkhoff-projektive Kontraktion + KL-Minimierung zwischen RG-Skalen).

2. **Kirk-Ergebnisse auf SU(N) erweitern:** Peter-Weyl-Koeffizienten a_pi(beta)
   fuer alle Darstellungen kontrollieren.

3. **Topologische Sektoren als Stabilitaetsargument:** Instanton-Sektoren (nu in Z)
   und deren Holley-Stroock-Invarianz unter RG-Fluss formalisieren.

4. **Reframing:** Paper als "Strong-Coupling Gap via Poincare/Dobrushin" positionieren,
   Kontinuumslimes als conditional programme labeln.


===============================================================================
## Verbindung zum Meta-Framework (FST Positivity Pattern)
===============================================================================

### Pattern A Stability Logic

Yang-Mills realisiert das FST-Positivity-Pattern auf drei Ebenen:

| Ebene | Mechanismus | YM-Instanziierung |
|-------|-------------|-------------------|
| Leading neutral | Free Energy konvex | F[mu] strikt konvex (trivial) |
| 1st order degenerate | KL trivial | mu_Lambda als Minimierer (uninteressant) |
| 2nd order decides | Hessian/Transfer Gap | Delta = log(lambda_0/lambda_1) >= kappa |

**Kernidee:** Die Massenluecke ist ein genuines Second-Order-Resolvent-Phaenomen.
Sie emergiert aus der Entkopplung topologisch verschiedener Link-Konfigurationen
via Dobrushin-Zegarlinski, analog zur Even/Odd-Entkopplung bei RH.

### Topologische Sektoren als Shift-Parity-Analogon

Das RH Shift Parity Lemma zeigt: Sektor-Labels (even/odd) bleiben unter
spektralen Shifts invariant. In der Gitter-Eichtheorie operiert ein analoger
Mechanismus via Instanton-Sektoren: Holley-Stroock-Konstanten bleiben
sektorweise invariant unter dem RG-Fluss beta(a). Kirks Zero-Free-Strip ist
das YM-Analogon der Shift-Parity.

### Gemeinsame Programmlinie

```
F mit Entropieterm -> Strikte Konvexitaet -> Eindeutiger Minimierer -> Physik
```

Bei YM ist die "Physik" = Massenluecke. Die Nichttrivialitaet liegt NICHT in
der KL-Konvexitaet (die ist trivial), sondern in der Poincare/Dobrushin-
Maschinerie, die die Free-Energy-Kruemmung in eine Spektralluecke uebersetzt.


===============================================================================
## Referenzen
===============================================================================

- Jaffe & Witten (2000): Millennium-Formulierung
- Wilson (1974): Gitter-Eichtheorie
- Osterwalder & Schrader (1973, 1975): OS-Rekonstruktion
- Osterwalder & Seiler (1978): Nichtnegative PW-Koeffizienten
- Holley & Stroock (1987): Bounded Perturbation Lemma
- Zegarlinski (1992): Log-Sobolev auf dem Gitter
- Balaban (1984, 1985): UV-Stabilitaet via RG
- Kirk (2026): SU(2) Kontinuumslimes, Zero-Free Strip
- Garcia Baquero (2026): Three Matter Theory (alternatives Programm)
