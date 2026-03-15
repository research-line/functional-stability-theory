# Naming Hierarchy -- Natur&Technik Forschungsprogramm

> Stand: 2026-03-15
> Zweck: Interne Referenz fuer die Beziehungen zwischen den Projekten

## Chronologische Entstehung

```
2025-2026   RH-Trilogie          CRM I-IV
            (eigenstaendig)       (eigenstaendig)
                  |                     |
                  v                     v
            +-----------+         +-----------+
            | Riemann   |         | Cosmic    |
            | Hypothesis|         | Recursion |
            | Part I-III|         | Model I-V |
            +-----------+         +-----------+
                  |                     |
                  +----------+----------+
                             |
                    Abstraktion / Generalisierung
                             |
                             v
                  +---------------------+
                  | RFEP Framework      |
                  | (Renormalized Free- |
                  | Energy Principle)   |
                  | = Verbindungsglied  |
                  +---------------------+
                             |
                    Instantiierung / Anwendung
                             |
                  +----------+----------+
                  |                     |
                  v                     v
       +------------------+   +------------------+
       | Domain Proofs    |   | Applications     |
       | (Folgebeweise)   |   | (Anwendungsfaelle)|
       +------------------+   +------------------+
       | NS, YM, TU, DE  |   | FST-I  Thermo    |
       | Hodge, BSD, PNP  |   | FST-II Chemical  |
       +------------------+   | FST-III Biology  |
                              +------------------+
                                      |
                              unter dem Dachnamen
                                      |
                                      v
                        +---------------------------+
                        | FST = Functional          |
                        | Stability Theory          |
                        | (Programmname, kam ZULETZT)|
                        +---------------------------+
```

## Abhaengigkeitskette (Beweisrichtung)

```
RH  ──────────> RFEP Framework ──────────> FST Domain Proofs
(bewiesen)      (abstrahiert aus RH)       (instantiieren RFEP)
                       |
CRM ──────────> RFEP Framework ──────────> FST Applications
(Modell)        (bestaetigt durch CRM)     (bestaetigen RFEP empirisch)
```

Wichtig: Die Pfeile zeigen die LOGISCHE Abhaengigkeit, nicht die zeitliche.
- RH steht OHNE RFEP/FST auf eigenen Beinen
- RFEP referenziert RH als "Reference Instantiation"
- FST referenziert RFEP als theoretische Grundlage
- Anwendungsfaelle BESTAETIGEN RFEP-Vorhersagen (nicht umgekehrt)

## 4-Repo-Struktur (GitHub)

```
github.com/lukisch/
    |
    +-- rh-even-dominance/            Repo 1: RH (eigenstaendig)
    |   +-- RH_Part_I_*.tex (EN+DE)
    |   +-- RH_Part_II_*.tex (EN+DE)
    |   +-- RH_Part_III_*.tex (EN+DE)
    |   +-- Zenodo DOI: 10.5281/zenodo.19035845
    |
    +-- crm-cosmology/                Repo 2: CRM (eigenstaendig)
    |   +-- CRM_Paper_I-V (EN+DE)
    |   +-- Zenodo DOI: 10.5281/zenodo.18728936
    |
    +-- rfep-framework/               Repo 3: RFEP (Verbindungsglied) NEU
    |   +-- FST_Unified_*.tex (EN+DE)
    |   +-- Zenodo DOI: 10.5281/zenodo.19039190
    |
    +-- functional-stability-theory/  Repo 4: FST (Programm)
        +-- domain-proofs/
        |   +-- navier-stokes/
        |   +-- yang-mills/
        |   +-- turbulence/
        |   +-- dark-energy/
        |   +-- hodge/
        |   +-- bsd/
        |   +-- p-vs-np/
        +-- applications/
        |   +-- fst-i-particles/
        |   +-- fst-ii-chemistry/
        |   +-- fst-iii-biology/
        +-- fst_references.bib
        +-- NAMING_HIERARCHY.md       (diese Datei, oeffentlich)
```

## Namenskonventionen

| Ebene | Name | Akronym | Bedeutung |
|-------|------|---------|-----------|
| Prinzip | Renormalized Free-Energy Principle | RFEP | Das mathematische Kernprinzip |
| Muster | Pattern A: Second-Order Resolvent Dominance | Pattern A | Das universelle Stabilitaetsmuster |
| Programm | Functional Stability Theory | FST | Der Programmname (Dach ueber alles) |
| Fundament | Riemann Hypothesis Proof | RH | Eigenstaendiger Beweis, Referenz-Instantiierung |
| Fundament | Cosmic Recursion Model | CRM | Eigenstaendiges Modell |

## Was FST NICHT ist

- FST ist NICHT "Fractal Game Theory" (FGT) -- Fraktale spielen keine zentrale Rolle
- FST ist NICHT "Free-Energy Spectral Theory" -- das war der historische RH-Name
- FST ist NICHT das RFEP selbst -- RFEP ist das Prinzip, FST ist das Programm
- Die Anwendungsfaelle (FST-I/II/III) BEGRUENDEN nicht das RFEP, sie BESTAETIGEN es

## Zenodo-Publikationen (Stand 2026-03-15)

| Record | Titel | Version | DOI |
|--------|-------|---------|-----|
| RH | The Riemann Hypothesis: A Gauge-Theoretic Proof via Even Dominance and Spectral Stability | v1.1 | 19035845 |
| CRM | Cosmic Recursion Model I-IV | v6 | 18728936 |
| CRM5 | The Saturation Theorem | v1.0 | 19036189 |
| RFEP | The Renormalized Free Energy Framework | v1.2 | 19039190 |
| DE | Dark Energy as Residual Vacuum Free Energy | v1.2 | 19038011 |
