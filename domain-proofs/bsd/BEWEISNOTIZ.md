# BEWEISNOTIZ -- BSD via Hoehen-Saettigung und Entropie-Normalform
# Stand: 2026-03-15 (nach Sprint 3 Review)
# Status: DRAFT -- Rang <= 1 bewiesen, Paper als "structural characterisation" repositioniert
# Review: 4C/8M/10L Issues identifiziert und korrigiert (Sprint 3)
# WICHTIG: Zirkularitaet in Prop 4.2 aufgeloest, Axiome als teilweise BSD-aequivalent deklariert

===============================================================================
## Problemstellung
===============================================================================

**Millennium-Problem (Birch & Swinnerton-Dyer 1965):**

**Schwache Form:** Fuer eine elliptische Kurve E/Q:
  ord_{s=1} L(E,s) = rank E(Q)

**Starke Form:** Die Leitkoeffizient-Formel:
  L^{(r)}(E,1) / r! = Omega_E * R_E * prod c_p * |Sha| / |E(Q)_tors|^2

wobei r = rank E(Q), Omega_E = reale Periode, R_E = Regulator
(Determinante der Neron-Tate-Paarung), c_p = Tamagawa-Zahlen,
Sha = Tate-Shafarevich-Gruppe.

**Bekannt fuer:**
- Rang 0: L(E,1) != 0 => rank = 0, |Sha| < inf (Kolyvagin 1990)
- Rang 1: ord_{s=1} L = 1 => rank = 1, |Sha| < inf (Gross-Zagier + Kolyvagin)
- Rang >= 2: OFFEN

**Ansatz (FST-BSD):** BSD als Positivity Normal-Form Theorem: Die BSD-Formel
ist das EINZIGE "Saettigungsprofil" das mit vier axiomatischen Constraints
kompatibel ist. Rang-Diskrepanz wuerde Positivitaet verletzen.


===============================================================================
## Beweiskette
===============================================================================

### Schritt 1: Neron-Tate-Hoehe als kanonische Positivitaets-Struktur
**Status: BEWIESEN (Theorem, nicht Vermutung)**

Neron-Tate-Hoehe h_hat: E(Q) -> R_>=0:
- h_hat(P) >= 0 fuer alle P, mit Gleichheit gdw P Torsion
- h_hat(nP) = n^2 h_hat(P) (quadratische Skalierung)
- Neron-Tate-Paarung <P,Q> = (1/2)(h_hat(P+Q) - h_hat(P) - h_hat(Q))
  ist positiv semi-definit auf E(Q), positiv definit auf E(Q)/E(Q)_tors

Regulator: R_E = det(<P_i, P_j>) > 0 fuer r >= 1 (Gram-Determinante
einer positiv definiten Form).


### Schritt 2: Die vier Axiome
**Status: FORMULIERT (Section 3)**

**(Axiom I) Endliche Arithmetische Kapazitaet:**
|Sha(E/Q)| < inf und |Sel_n(E/Q)| < inf fuer alle n.
Status: Bewiesen fuer Rang <= 1 (Kolyvagin). Offen fuer Rang >= 2.

**(Axiom II) Kooperative Hoehen-Positivitaet:**
<P,P> > 0 fuer P nicht-Torsion. Hoehen kombinieren quadratisch: R_E > 0.
Status: UNBEDINGT WAHR (Theorem, nicht Vermutung).

**(Axiom III) Monotone Iwasawa-Kontrolle:**
Im zyklotomischen Z_p-Turm Q_inf/Q wachsen Selmer-Gruppen kontrolliert:
  rank_{Z_p} Sel_{p^inf}(E/Q_n) = mu*p^n + lambda*n + nu + O(1)
Iwasawa Main Conjecture liefert p-adische L-Funktion L_p(E,s).
Status: Bewiesen in vielen Faellen (Skinner-Urban 2014, Kato 2004).

**(Axiom IV) Euler-System-Kompositionskohaerenz:**
Kohaerente Familie {c_m} von Kohomologie-Klassen mit Norm-Kompatibilitaet:
  cores_{m*ell/m}(c_{m*ell}) = P_ell(Frob_ell^{-1}) * c_m
Obere Schranken fuer Selmer-Gruppen via Kolyvagin-Derivat-Maschinerie.
Status: Euler-Systeme existieren (Heegner-Punkte, Kato, Beilinson-Flach).


### Schritt 3: Normal-Form Theorem
**Status: FORMULIERT (Section 4)**

Die BSD-Leitkoeffizienten-Formel ist die EINZIGE Struktur, die alle
vier Axiome gleichzeitig erfuellt.

**Argument:** Betrachte Alternativen:
- Analytische Nullstellen OHNE geometrische Richtungen:
  Verletzt Axiom IV (Euler-System-Kohaerenz erzwingt Punkte)
- Geometrischer Rang OHNE analytische Ableitungen:
  Verletzt Axiom III (Iwasawa-Kontrolle koppelt analytisch/algebraisch)
- Falscher Leitkoeffizient (z.B. R_E durch andere Invariante ersetzt):
  Verletzt Axiom II (nur Neron-Tate ist die kanonische positiv-definite Form)
- Unendliches Sha:
  Verletzt Axiom I (absorbiert Rang-Diskrepanz)


### Schritt 4: Rang <= 1 -- Gross-Zagier als Positivitaets-Identitaet
**Status: BEWIESEN (Gross-Zagier 1986 + Kolyvagin 1990)**

**Rang 0:** L(E,1) != 0 => (via Kolyvagin Euler-System)
  |Sel(E/Q)| beschraenkt => rank E(Q) = 0, |Sha| < inf.
  Alle vier Axiome erfuellt.

**Rang 1:** Gross-Zagier-Theorem:
  L'(E,1) = c * h_hat(P_K)
mit Heegner-Punkt P_K. Dies ist eine POSITIVITAETS-IDENTITAET:
  - Linke Seite: Analytisches Objekt (Ableitung der L-Funktion)
  - Rechte Seite: Nichtnegativer geometrischer Wert (Hoehe)
  - L'(E,1) != 0 <=> h_hat(P_K) > 0 <=> P_K nicht-Torsion

Kolyvagin vervollstaendigt: P_K nicht-Torsion => rank = 1, |Sha| < inf.


### Schritt 5: Height Saturation Conjecture
**Status: OFFEN (Section 5)**

Fuer allgemeinen Rang r:
  L^{(r)}(E,1) / r! ~ det(<P_i, P_j>) = R_E

wobei P_1, ..., P_r unabhaengige Punkte sind.

Nichtnegativitaet von det(<P_i, P_j>) (aus positiver Definitheit von h_hat)
erzwingt: ord_{s=1} L(E,s) = rank E(Q).

**Dies ist das "Higher Gross-Zagier": eine r-dimensionale Verallgemeinerung
der Rang-1-Positivitaets-Identitaet.**


===============================================================================
## Offene Luecken
===============================================================================

### L1: Higher Gross-Zagier (KERNLUECKE)
L^{(r)}(E,1) ~ R_E fuer r >= 2 ist NICHT bewiesen.
Es fehlt eine universelle Formel, die den r-ten Ableitungswert der
L-Funktion an eine Gram-Determinante von Hoehen koppelt.

Partielle Resultate:
- Darmon-Rotger: Diagonal-Restriction fuer Rang 2 (Spezialfaelle)
- Loeffler-Skinner-Zerbes: Euler-Systeme fuer hoehere Motive

### L2: Vollstaendige Kolyvagin-Systeme fuer hoeheren Rang
Fuer Rang 1 liefern Heegner-Punkte ein vollstaendiges Euler-System.
Fuer Rang >= 2 existiert kein analoges System.

### L3: Sha-Endlichkeit fuer Rang >= 2
|Sha(E/Q)| < inf (Axiom I) ist nur fuer Rang <= 1 bewiesen.
Ohne Sha-Endlichkeit koennte "versteckter Rang" in Sha die
Rang-Gleichheit unterlaufen.

### L4: Uniformity Bridge
Die Verbindung zwischen dem p-adischen Main Conjecture (Axiom III)
und der globalen Rang-Gleichheit (schwaches BSD) ist nicht geschlossen.
Der Uebergang von p-adischer zu archimedischer Information ist
der fundamentale Engpass.


===============================================================================
## Naechste Schritte
===============================================================================

1. **Darmon-Rotger Diagonal-Restriction studieren:** Rang-2-Spezialfaelle
   als Testfall fuer die Height Saturation Conjecture.

2. **Loeffler-Skinner-Zerbes:** Euler-Systeme fuer hoehere Motive
   auf Kompatibilitaet mit den vier Axiomen pruefen.

3. **Numerische Verifikation:** BSD-Formel fuer Kurven mit Rang 2-4
   numerisch verifizieren (LMFDB-Datenbank).

4. **Verbindung zu Hodge formalisieren:** Gleiche Positivitaets-Architektur
   (Easy Direction bewiesen, Hard Direction offen) explizit machen.


===============================================================================
## Verbindung zum Meta-Framework (FST Positivity Pattern)
===============================================================================

### BSD als Positivity Normal-Form

Die BSD-Formel ist das arithmetische Analogon der Variationsprinzipien
in den physikalischen Folgebeweisen:

| Physik-Problem | Positiv-Form | BSD-Analogon |
|----------------|-------------|--------------|
| YM: Free Energy | F[mu] strikt konvex | h_hat positiv definit |
| NS: H(u|A) >= 0 | Attractor-Distanz | Rang-Distanz ord - rank |
| Turbulenz: F[E] >= 0 | KL-Divergenz | BSD-Leitkoeffizient >= 0 |
| DE: Phi[rho] -> min | Vakuum-Minimum | Regulator R_E > 0 |

### Die Neron-Tate-Hoehe als kanonische Positivitaet

Die Hoehe spielt fuer BSD dieselbe Rolle wie:
- Die Hodge-Riemann-Form fuer die Hodge-Vermutung
- Die Frobenius-Positivitaet fuer die Weil-Vermutungen
- Die spektrale Positivitaet fuer RH

### Easy Direction / Hard Direction Pattern

```
Positivitaet (Theorem) -> Normal-Form (Rang <= 1 bewiesen) -> Hard Direction (Rang >= 2 offen)
```

Dieses Pattern teilt BSD mit Hodge:
- Easy: "algebraisch/geometrisch => positiv" (bewiesen)
- Hard: "positiv => algebraisch/geometrisch" (offen)

### Gemeinsame Programmlinie mit Hodge

BSD und Hodge sind die beiden "arithmetischen" Folgebeweise.
Beide basieren auf positiv definiten Formen (Neron-Tate bzw. Hodge-Riemann),
beide haben die Easy Direction als Theorem und die Hard Direction offen.
Ein Durchbruch bei einem koennte den anderen katalysieren.


===============================================================================
## Referenzen
===============================================================================

- Birch & Swinnerton-Dyer (1965): Original-Vermutung
- Gross & Zagier (1986): Heegner-Punkte und L-Funktions-Ableitungen
- Kolyvagin (1990): Euler-Systeme und Rang-Schranken
- Kato (2004): p-adische Hodge-Theorie und Euler-Systeme
- Skinner & Urban (2014): Iwasawa Main Conjecture fuer GL(2)
- Loeffler, Lei & Zerbes (2014): Euler-Systeme fuer Rankin-Selberg
- Darmon & Rotger: Diagonal-Restriction (Rang-2-Spezialisierungen)
- Silverman (1986, 1994): Arithmetic of Elliptic Curves
