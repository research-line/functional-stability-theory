# Functional Stability Theory (FST)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19036190.svg)](https://doi.org/10.5281/zenodo.19036190)

**Functional Stability Theory** is a unified mathematical framework that connects thermodynamic stability, game-theoretic equilibrium selection, and renormalized energy functionals to derive rigorous results across multiple domains of mathematical physics and complexity theory.

The central object is a **renormalized free-energy functional** whose critical-point structure, combined with selection principles from evolutionary game theory, yields domain-specific theorems when instantiated with the appropriate state spaces and energy densities. The framework is governed by the **Dissipative Selection Principle** (DSP): among all critical points of the renormalized functional, physical evolution selects the unique dissipation-minimising attractor -- a meta-theorem that unifies the domain-specific results below.

## Repository Structure

- **masters/** -- Foundation papers (three masters)
  - **zookeeper/** -- The Spectral Zookeeper: RH via microcluster closure ([DOI](https://doi.org/10.5281/zenodo.19673127))
  - **zeta-zoo/** -- The Zeta Zoo: HP-BL classification, SGE taxonomy ([DOI](https://doi.org/10.5281/zenodo.19673227))
  - **spectrum-duality/** -- FST Spectrum Duality / RFEP ([DOI](https://doi.org/10.5281/zenodo.19162705))
- **domain-proofs/** -- Domain-specific instantiations and proofs
- **scripts/** -- Numerical validation scripts
- **fst_references.bib** -- Shared bibliography

## Domain Proof Papers -- Current Status

| Paper | Version | Status | Open Problem | DOI |
|-------|---------|--------|--------------|-----|
| [**Turbulence**](https://doi.org/10.5281/zenodo.19056813) | v1.3 | Journal-ready | DFC1 empirical (only input) | [10.5281/zenodo.19056813](https://doi.org/10.5281/zenodo.19056813) |
| [**Dark Energy**](https://doi.org/10.5281/zenodo.19036235) | v1.6 | Framework Note | Hu–Sawicki parameters quantitatively open | [10.5281/zenodo.19036235](https://doi.org/10.5281/zenodo.19036235) |
| [**Yang–Mills**](https://doi.org/10.5281/zenodo.19087433) | v2.1 | Conditional | Analytical proof of λ < 0 | [10.5281/zenodo.19087433](https://doi.org/10.5281/zenodo.19087433) |
| [**Navier–Stokes**](https://doi.org/10.5281/zenodo.19087449) | v2.1 | Conditional | Assumption G2 (projection regularity) | [10.5281/zenodo.19087449](https://doi.org/10.5281/zenodo.19087449) |
| [**NS Log-Distance**](https://doi.org/10.5281/zenodo.19056807) | v1.3 | Proof of Life ✓ | TLL for 3D NS analytically open | [10.5281/zenodo.19056807](https://doi.org/10.5281/zenodo.19056807) |
| [**BSD**](https://doi.org/10.5281/zenodo.19087443) | v1.1 | Reformulation | Higher Gross–Zagier (rank ≥ 2) | [10.5281/zenodo.19087443](https://doi.org/10.5281/zenodo.19087443) |
| [**Hodge**](https://doi.org/10.5281/zenodo.19087439) | v1.1 | No-Go Theorem | = Deligne's question (1982) | [10.5281/zenodo.19087439](https://doi.org/10.5281/zenodo.19087439) |
| [**P vs NP**](https://doi.org/10.5281/zenodo.19056809) | v1.2 | Reformulation | Uniformity Bridge | [10.5281/zenodo.19056809](https://doi.org/10.5281/zenodo.19056809) |
| [**Framework (RFEP)**](https://doi.org/10.5281/zenodo.19036190) | v1.6 | Meta-Theorem | Pattern A falsifiability clarified | [10.5281/zenodo.19036190](https://doi.org/10.5281/zenodo.19036190) |

## Proof Architecture

```
                    Framework (Pattern A)
                    +-- DS1-DS3 (Dissipative Selection)
                    +-- Second-Order Resolvent Dominance
                         |
          +--------------+--------------+
          |              |              |
     PROVEN        CONDITIONAL        OPEN
     (rigorous)   (reduced to       (= open
                  threshold axiom)   research)
          |              |              |
         TU:            YM:           BSD:
      DFC => NL'    Doeblin a > 0   Rank >= 2
      (journal)     (Kingman l < 0)  Gross-Zagier
                        |
         DE:           NS:          Hodge:
      Screening    G2 (projection)  Deligne's
      (validated)  G3 (Gronwall)    question
                        |
        NS-LDI:       PvNP:
      TLL + LDI     Uniformity
      (Lorenz OK)    Bridge
```

## Numerical Validation Scripts

The `scripts/` directory contains computational validation scripts:

| Script | Paper | Description |
|--------|-------|-------------|
| `scripts/turbulence/compute_F_spectrum.py` | Turbulence | Verifies K41 as unique minimiser of F[E]; strict convexity test |
| `scripts/navier-stokes/compute_ds3_lorenz.py` | Navier-Stokes | DS3 stress test on Lorenz attractor; TV saturation |
| `scripts/navier-stokes/compute_tll_ldi_lorenz.py` | NS-LDI | **Proof of Life**: TLL+LDI on Lorenz attractor (5/5 tests passed) |
| `scripts/dark-energy/compute_w_vs_desi.py` | Dark Energy | w_eff(z) comparison with DESI constraints |
| `scripts/dark-energy/compute_w_mapping.py` | Dark Energy | Correct w_eff -> w_DE mapping + DESI grid scan |
| `scripts/dark-energy/compute_husawicki_mcmc.py` | Dark Energy | Hu-Sawicki f(R) MCMC fit against DESI+Planck+Cassini |
| `scripts/bsd/compute_height_saturation.py` | BSD | Height saturation test for quadratic twists |
| `scripts/hodge/compute_ghr_spectrum.py` | Hodge | GHR spectrum numerical verification |

## Genesis and Architecture

### Chronological Development

```
2025-2026   RH Trilogy             CRM I-IV
            (self-contained)       (self-contained)
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
                    Abstraction / Generalisation
                             |
                             v
                  +---------------------+
                  | RFEP Framework      |
                  | (Renormalized Free- |
                  | Energy Principle)   |
                  | = Connecting Link   |
                  +---------------------+
                             |
                    Instantiation / Application
                             |
                  +----------+----------+
                  |                     |
                  v                     v
       +------------------+   +------------------+
       | Domain Proofs    |   | Applications     |
       +------------------+   +------------------+
       | NS, YM, TU, DE  |   | FST-I  Thermo    |
       | Hodge, BSD, PNP  |   | FST-II Chemical  |
       +------------------+   | FST-III Biology  |
                              +------------------+
                                      |
                              under the umbrella name
                                      |
                                      v
                        +---------------------------+
                        | FST = Functional          |
                        | Stability Theory          |
                        | (programme name, came LAST)|
                        +---------------------------+
```

### Dependency Chain (Proof Direction)

```
RH  ──────────> RFEP Framework ──────────> FST Domain Proofs
(proven)        (abstracted from RH)       (instantiate RFEP)
                       |
CRM ──────────> RFEP Framework ──────────> FST Applications
(model)         (confirmed by CRM)         (validate RFEP empirically)
```

Arrows indicate **logical** dependency, not chronological order:
- RH stands on its own without RFEP/FST
- RFEP references RH as its "Reference Instantiation"
- FST references RFEP as its theoretical foundation
- Applications **confirm** RFEP predictions (not the other way around)

### Naming Conventions

| Level | Name | Acronym | Meaning | DOI |
|-------|------|---------|---------|-----|
| Principle | Renormalized Free-Energy Principle | RFEP | The core mathematical principle | [10.5281/zenodo.19036190](https://doi.org/10.5281/zenodo.19036190) |
| Pattern | Pattern A: Second-Order Resolvent Dominance | Pattern A | The universal stability pattern | [10.5281/zenodo.19036190](https://doi.org/10.5281/zenodo.19036190) |
| Programme | Functional Stability Theory | FST | The programme name (umbrella over all) | -- |
| Master | The Spectral Zookeeper | Zookeeper | RH proof via CCM microcluster closure | [10.5281/zenodo.19673127](https://doi.org/10.5281/zenodo.19673127) |
| Master | The Zeta Zoo | Zeta Zoo | HP-BL classification, SGE taxonomy | [10.5281/zenodo.19673227](https://doi.org/10.5281/zenodo.19673227) |
| Foundation | RH Even Dominance | RH v2.1 | Independent RH proof, even dominance route | [10.5281/zenodo.19546593](https://doi.org/10.5281/zenodo.19546593) |
| Foundation | Cosmic Recursion Model | CRM | Self-contained model | [10.5281/zenodo.18728935](https://doi.org/10.5281/zenodo.18728935) |

## Related Repositories

- [rh-even-dominance](https://github.com/research-line/rh-even-dominance) -- Riemann Hypothesis: Even-dominance proof (foundation)
- [crm-cosmology](https://github.com/research-line/crm-cosmology) -- Curvature Relaxation Model (foundation)
- [rfep-framework](https://github.com/research-line/rfep-framework) -- archived, now integrated here under `framework/`

## Author

Lukas Geiger -- ORCID: [0009-0005-7296-1534](https://orcid.org/0009-0005-7296-1534)

## License

This work is licensed under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).

---

## Haftung / Liability

Dieses Projekt ist eine **unentgeltliche Open-Source-Schenkung** im Sinne der §§ 516 ff. BGB. Die Haftung des Urhebers ist gemäß **§ 521 BGB** auf **Vorsatz und grobe Fahrlässigkeit** beschränkt. Ergänzend gelten die Haftungsausschlüsse aus GPL-3.0 / MIT / Apache-2.0 §§ 15–16 (je nach gewählter Lizenz).

Nutzung auf eigenes Risiko. Keine Wartungszusage, keine Verfügbarkeitsgarantie, keine Gewähr für Fehlerfreiheit oder Eignung für einen bestimmten Zweck.

This project is an unpaid open-source donation. Liability is limited to intent and gross negligence (§ 521 German Civil Code). Use at your own risk. No warranty, no maintenance guarantee, no fitness-for-purpose assumed.

