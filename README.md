# Functional Stability Theory (FST)

**Functional Stability Theory** is a unified mathematical programme that identifies a single structural challenge — *Functional Positivity under Gauge Constraint* (Pattern A) — as the common substrate of open problems in number theory, mathematical physics, and cosmology.

## The Three Masters

The programme rests on three foundation papers:

| Master | Title | Role | DOI |
|--------|-------|------|-----|
| [**Zookeeper**](masters/zookeeper/) | The Spectral Zookeeper | RH proof via CCM microcluster closure | [10.5281/zenodo.19673127](https://doi.org/10.5281/zenodo.19673127) |
| [**Zeta Zoo**](masters/zeta-zoo/) | The Zeta Zoo | Mathematical classification (SGE taxonomy, Boundary Theorem) | [10.5281/zenodo.19673227](https://doi.org/10.5281/zenodo.19673227) |
| [**Spectrum Duality**](masters/spectrum-duality/) | FST Spectrum Duality / RFEP | Physical instantiation (Pattern A, DS1-DS3) | [10.5281/zenodo.19162705](https://doi.org/10.5281/zenodo.19162705) |

## Domain Supplements — Mathematics (Zeta Zoo)

Classified by the SGE taxonomy from the Zeta Zoo. These instantiate Pattern A on number-theoretic / algebraic structures.

| Paper | Version | Status | Open Problem | DOI |
|-------|---------|--------|--------------|-----|
| [**BSD**](https://doi.org/10.5281/zenodo.19087443) | v1.1 | Reformulation | Higher Gross–Zagier (rank ≥ 2) | [10.5281/zenodo.19087443](https://doi.org/10.5281/zenodo.19087443) |
| [**Hodge**](https://doi.org/10.5281/zenodo.19087439) | v1.1 | No-Go Theorem | = Deligne's question (1982) | [10.5281/zenodo.19087439](https://doi.org/10.5281/zenodo.19087439) |
| [**P vs NP**](https://doi.org/10.5281/zenodo.19056809) | v1.2 | Reformulation | Uniformity Bridge | [10.5281/zenodo.19056809](https://doi.org/10.5281/zenodo.19056809) |

## Domain Supplements — Physics (Spectrum Duality)

Derive Pattern A + DS1–DS3 from Spectrum Duality. These instantiate the Dissipative Selection Principle on physical systems.

| Paper | Version | Status | Open Problem | DOI |
|-------|---------|--------|--------------|-----|
| [**Turbulence**](https://doi.org/10.5281/zenodo.19056813) | v1.3 | Journal-ready | DFC1 empirical (only input) | [10.5281/zenodo.19056813](https://doi.org/10.5281/zenodo.19056813) |
| [**Dark Energy**](https://doi.org/10.5281/zenodo.19036235) | v1.6 | Framework Note | Hu–Sawicki parameters quantitatively open | [10.5281/zenodo.19036235](https://doi.org/10.5281/zenodo.19036235) |
| [**Yang–Mills**](https://doi.org/10.5281/zenodo.19087433) | v2.1 | Conditional | Analytical proof of λ < 0 | [10.5281/zenodo.19087433](https://doi.org/10.5281/zenodo.19087433) |
| [**Navier–Stokes**](https://doi.org/10.5281/zenodo.19087449) | v2.1 | Conditional | Assumption G2 (projection regularity) | [10.5281/zenodo.19087449](https://doi.org/10.5281/zenodo.19087449) |
| [**NS Log-Distance**](https://doi.org/10.5281/zenodo.19056807) | v1.3 | Proof of Life | TLL for 3D NS analytically open | [10.5281/zenodo.19056807](https://doi.org/10.5281/zenodo.19056807) |

## Proof Architecture

```
                     THREE MASTERS
    ┌────────────────────┼────────────────────┐
    │                    │                    │
 Zookeeper          Zeta Zoo          Spectrum Duality
 (RH proof)      (Classification)    (Pattern A, DS1-DS3)
    │                    │                    │
    │              Math Supplements     Physics Supplements
    │                    │                    │
    │              ┌─────┴─────┐        ┌────┴────┐
    │              │     │     │        │    │    │
    │             BSD  Hodge  PvNP     TU   YM   NS
    │           (rank≥2)(No-Go)(Unif.) (✓) (cond)(cond)
    │                                       │
    │                                      DE    NS-LDI
    │                                   (valid.) (PoL)
    │
 Status: PROVEN (unconditional, CCM route)
```

### Detailed Proof Status

```
  PROVEN (rigorous)     CONDITIONAL (threshold)     OPEN (= research)
       │                       │                         │
      TU:                     YM:                      BSD:
   DFC⇒NL'                Doeblin                   Rank ≥ 2
   (journal-ready)        α(β)>0                    Gross-Zagier
                              │
      DE:                    NS:                     Hodge:
   Screening              G₂ (proj.)               Deligne's
   (validated)            G₃ (Gronwall)            question
                              │
     NS-LDI:                PvNP:
   TLL+LDI                Uniformity
   (Lorenz ✓)              Bridge
```

## Genesis and Architecture

### Chronological Development

```
2025–2026   RH Trilogy            CRM I–IV
            (independent)         (independent)
                  |                     |
                  v                     v
            +-----------+         +-----------+
            | Riemann   |         | Cosmic    |
            | Hypothesis|         | Recursion |
            | Parts I–III|        | Model I–V |
            +-----------+         +-----------+
                  |                     |
                  +----------+----------+
                             |
                    Abstraction / Generalization
                             |
                             v
            +────────────────────────────────────+
            │       THREE MASTERS (FST)          │
            ├────────────────────────────────────┤
            │ Zookeeper    │ RH proof (CCM)      │
            │ Zeta Zoo     │ Classification (SGE)│
            │ Spectrum Duality │ Physics (RFEP)  │
            +────────────────────────────────────+
                             |
                    Instantiation / Application
                             |
                  +----------+----------+
                  |                     |
                  v                     v
       +------------------------+   +------------------------+
       | Domain Proofs          |   | Applications           |
       | (follow-up proofs)     |   | (use cases)            |
       +------------------------+   +------------------------+
       | NS, YM, TU, DE        |   | FST-I  Thermo          |
       | Hodge, BSD, P vs NP   |   | FST-II Chemical        |
       +------------------------+   | FST-III Biology        |
                                    +------------------------+
                                             |
                                    under the umbrella name
                                             |
                                             v
                              +-----------------------------+
                              | FST = Functional            |
                              | Stability Theory            |
                              | (program name, introduced   |
                              | last)                       |
                              +-----------------------------+
```

### Dependency Chain (Proof Direction)

```
RH  ──────────> Zookeeper (Master)  ──────> Zeta Zoo (classification)
(proven)        (CCM proof route)           (SGE taxonomy of all zetas)
                       │
CRM ──────────> Spectrum Duality (Master) ──> Domain Proofs
(model)         (Pattern A, RFEP)             (instantiate Pattern A)
```

Arrows indicate **logical** dependency, not chronological order:
- RH and CRM stand on their own (independent foundations)
- The Zookeeper proves RH via CCM microcluster closure
- The Zeta Zoo classifies all zeta-type families (including Riemann zeta)
- Spectrum Duality provides the physical axiom package (Pattern A, DS1–DS3)
- Domain supplements instantiate Pattern A from Spectrum Duality

### Core Foundations — Principles and Naming

| Level | Name | Acronym | Meaning | DOI |
|-------|------|---------|---------|-----|
| Principle | Renormalized Free-Energy Principle | RFEP | The mathematical core principle | [10.5281/zenodo.19162705](https://doi.org/10.5281/zenodo.19162705) |
| Pattern | Pattern A: Functional Positivity under Gauge Constraint | Pattern A | The universal stability pattern | [10.5281/zenodo.19162705](https://doi.org/10.5281/zenodo.19162705) |
| Programme | Functional Stability Theory | FST | The programme name (umbrella over all) | — |
| Master | The Spectral Zookeeper | Zookeeper | RH proof via CCM microcluster closure | [10.5281/zenodo.19673127](https://doi.org/10.5281/zenodo.19673127) |
| Master | The Zeta Zoo | Zeta Zoo | HP-BL classification, SGE taxonomy | [10.5281/zenodo.19673227](https://doi.org/10.5281/zenodo.19673227) |
| Master | FST Spectrum Duality | RFEP | Physical instantiation (Pattern A, DS1–DS3) | [10.5281/zenodo.19162705](https://doi.org/10.5281/zenodo.19162705) |
| Foundation | RH Even Dominance | RH v2.1 | Independent RH proof, second route | [10.5281/zenodo.19546593](https://doi.org/10.5281/zenodo.19546593) |
| Foundation | Cosmic Recursion Model | CRM | Independent dark energy model | [10.5281/zenodo.18728935](https://doi.org/10.5281/zenodo.18728935) |

## Independent Foundations

| Name | Role | DOI |
|------|------|-----|
| RH Even Dominance (v2.1) | Independent RH proof, second route | [10.5281/zenodo.19546593](https://doi.org/10.5281/zenodo.19546593) |
| CRM Cosmology (I-V) | Independent dark energy model | [10.5281/zenodo.18728935](https://doi.org/10.5281/zenodo.18728935) |

These stand independently of FST. The RFEP was abstracted from them; they are not derived from it.

## Numerical Validation Scripts

| Script | Paper | Description |
|--------|-------|-------------|
| `scripts/turbulence/compute_F_spectrum.py` | Turbulence | Verifies K41 as unique minimiser of F[E]; strict convexity test |
| `scripts/navier-stokes/compute_ds3_lorenz.py` | Navier-Stokes | DS3 stress test on Lorenz attractor; TV saturation |
| `scripts/navier-stokes/compute_tll_ldi_lorenz.py` | NS-LDI | **Proof of Life**: TLL+LDI on Lorenz attractor (5/5 tests passed) |
| `scripts/dark-energy/compute_w_vs_desi.py` | Dark Energy | w_eff(z) comparison with DESI constraints |
| `scripts/dark-energy/compute_w_mapping.py` | Dark Energy | Correct w_eff → w_DE mapping + DESI grid scan |
| `scripts/dark-energy/compute_husawicki_mcmc.py` | Dark Energy | Hu-Sawicki f(R) MCMC fit against DESI+Planck+Cassini |
| `scripts/bsd/compute_height_saturation.py` | BSD | Height saturation test for quadratic twists |
| `scripts/hodge/compute_ghr_spectrum.py` | Hodge | GHR spectrum numerical verification |

## Repository Structure

```
functional-stability-theory/
├── masters/                Three foundation papers
│   ├── zookeeper/          RH proof (microcluster closure)
│   ├── zeta-zoo/           Classification (SGE taxonomy)
│   └── spectrum-duality/   Physical axioms (RFEP, Pattern A)
├── domain-proofs/          Domain supplements (8 papers)
│   ├── turbulence/         K41 spectrum (unconditional)
│   ├── yang-mills/         Mass gap (conditional)
│   ├── navier-stokes/      Regularity (conditional)
│   ├── dark-energy/        CRM screening (validated)
│   ├── bsd/                Rank-1 positivity
│   ├── hodge/              No-go + easy direction
│   └── p-vs-np/            Witness entropy gap
└── scripts/                Numerical validation
```

## Related Repositories

- [rh-even-dominance](https://github.com/research-line/rh-even-dominance) — Riemann Hypothesis: Even-dominance proof (foundation)
- [crm-cosmology](https://github.com/research-line/crm-cosmology) — Curvature Relaxion Model (foundation)

## Author

Lukas Geiger -- ORCID: [0009-0005-7296-1534](https://orcid.org/0009-0005-7296-1534)

## License

[CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/)
