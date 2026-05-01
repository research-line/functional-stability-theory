# Functional Stability Theory (FST)

**Functional Stability Theory** is a unified mathematical programme that identifies a single structural challenge — *Functional Positivity under Gauge Constraint* (Pattern A) — as the common substrate of open problems in number theory, mathematical physics, and cosmology.

## The Three Masters

The programme rests on three foundation papers:

| Master | Title | Role | DOI |
|--------|-------|------|-----|
| [**Zookeeper**](masters/zookeeper/) | The Spectral Zookeeper | RH proof via CCM microcluster closure | [10.5281/zenodo.19673127](https://doi.org/10.5281/zenodo.19673127) |
| [**Zeta Zoo**](masters/zeta-zoo/) | The Zeta Zoo — The Mathematical Side of FST | Classification (SGE taxonomy, Boundary Theorem) | [10.5281/zenodo.19673227](https://doi.org/10.5281/zenodo.19673227) |
| [**Spectrum Duality**](masters/spectrum-duality/) | FST Spectrum Duality / RFEP | Physical instantiation (Pattern A, DS1–DS3) | [10.5281/zenodo.19162705](https://doi.org/10.5281/zenodo.19162705) |

## Domain Supplements

### FST-Mathematics

Classified by the SGE taxonomy from the Zeta Zoo. These instantiate Pattern A on number-theoretic and algebraic structures. BSD, Hodge, and P vs NP are *bridge species* — they appear in both the mathematical and physical branches of the programme.

| Paper | Version | Status | Open Problem | DOI |
|-------|---------|--------|--------------|-----|
| [**BSD**](fst-mathematics/bsd/) | v1.1 | Reformulation | Higher Gross–Zagier (rank ≥ 2) | [10.5281/zenodo.19087443](https://doi.org/10.5281/zenodo.19087443) |
| [**Hodge**](fst-mathematics/hodge/) | v1.1 | No-Go Theorem | = Deligne's question (1982) | [10.5281/zenodo.19087439](https://doi.org/10.5281/zenodo.19087439) |
| [**P vs NP**](fst-mathematics/p-vs-np/) | v1.2 | Reformulation | Uniformity Bridge | [10.5281/zenodo.19056809](https://doi.org/10.5281/zenodo.19056809) |

### FST-Physics

Derive Pattern A + DS1–DS3 from Spectrum Duality. These instantiate the Dissipative Selection Principle on physical systems.

| Paper | Version | Status | Open Problem | DOI |
|-------|---------|--------|--------------|-----|
| [**Turbulence**](fst-physics/turbulence/) | v1.3 | Journal-ready | DFC1 empirical (only input) | [10.5281/zenodo.19056813](https://doi.org/10.5281/zenodo.19056813) |
| [**Yang–Mills**](fst-physics/yang-mills/) | v2.1 | Conditional | Analytical proof of λ < 0 | [10.5281/zenodo.19087433](https://doi.org/10.5281/zenodo.19087433) |
| [**Navier–Stokes**](fst-physics/navier-stokes/) | v2.1 | Conditional | Assumption G2 (projection regularity) | [10.5281/zenodo.19087449](https://doi.org/10.5281/zenodo.19087449) |
| [**NS Log-Distance**](fst-physics/navier-stokes/) | v1.3 | Proof of Life | TLL for 3D NS analytically open | [10.5281/zenodo.19056807](https://doi.org/10.5281/zenodo.19056807) |

### FST-Cosmology

The cosmological branch of FST. The Dark Energy paper instantiates Pattern A on cosmological screening mechanisms (Hu–Sawicki f(R) gravity). A dedicated CRM paper re-examining the Cooperative Renormalization Model from the FST perspective is planned.

| Paper | Version | Status | Open Problem | DOI |
|-------|---------|--------|--------------|-----|
| [**Dark Energy**](fst-cosmology/dark-energy/) | v1.6 | Framework Note | Hu–Sawicki parameters quantitatively open | [10.5281/zenodo.19036235](https://doi.org/10.5281/zenodo.19036235) |

### FST-Biology

In development. See [`fst-biology/`](fst-biology/).

### FST-Chemistry

Planned. See [`fst-chemistry/`](fst-chemistry/).

## Proof Architecture

```
                        THREE MASTERS
     ┌──────────────────────┼──────────────────────┐
     │                      │                      │
  Zookeeper            Zeta Zoo            Spectrum Duality
  (RH proof)        (Classification)      (Pattern A, DS1-DS3)
     │                      │                      │
     │              FST-Mathematics          FST-Physics
     │                      │                      │
     │              ┌───────┼───────┐        ┌─────┼─────┐
     │              │       │       │        │     │     │
     │            BSD†    Hodge†  PvNP†      TU    YM    NS
     │                                                   │
     │                                      FST-Cosmology│
     │                                            │    NS-LDI
     │                                           DE
     │
  Status: PROVEN (unconditional, CCM route)

  † = bridge species (math + physics)
```

### Detailed Proof Status

```
  PROVEN (rigorous)     CONDITIONAL (threshold)     OPEN (= research)
       │                       │                         │
      TU:                     YM:                      BSD†:
   DFC⇒NL'                Doeblin                   Rank ≥ 2
   (journal-ready)        α(β)>0                    Gross-Zagier
                              │
      DE:                    NS:                     Hodge†:
   Screening              G₂ (proj.)               Deligne's
   (validated)            G₃ (Gronwall)            question
                              │
     NS-LDI:                PvNP†:
   TLL+LDI                Uniformity
   (Lorenz ✓)              Bridge

  † = bridge species
```

## Hierarchy

```
FST (Functional Stability Theory)
│
├── Masters
│   ├── Zookeeper         RH proof (CCM microcluster closure)
│   ├── Zeta Zoo          Mathematical classification (SGE taxonomy)
│   └── Spectrum Duality  Physical instantiation (RFEP, Pattern A)
│
├── FST-Mathematics       BSD, Hodge, P vs NP
├── FST-Physics           Turbulence, Yang–Mills, Navier–Stokes, NS-LDI
├── FST-Cosmology         Dark Energy
├── FST-Biology           (in development)
└── FST-Chemistry         (planned)
```

## Chronological Development

```
2025/2026  CRM I–IV (dark energy)     RH "light" proof (even dominance)
           developed independently    developed independently
                 \                       /
                  \                     /
                   +---------+---------+
                             |
                   Recognition: both share the same
                   structural pattern (Pattern A)
                             |
                             v
                   RFEP formulated (general principle)
                             |
                   Several dead ends
                             |
                             v
                   Idea: classify zeta-type families
                   using techniques from the RH proof
                             |
                   Not enough — need deeper tools
                             |
                             v
        RH via Connes framework (CCM)
        microcluster closure → unconditional proof
                             |
                             v
                   Zeta Zoo opens: SGE taxonomy
                   classifies all zeta families
                             |
                             v
                   Reorganization under FST
                   (Three Masters + five domain branches)
```

### Core Foundations — Principles and Naming

| Level | Name | Acronym | Meaning | DOI |
|-------|------|---------|---------|-----|
| Programme | Functional Stability Theory | FST | The programme name (umbrella over all) | — |
| Principle | Renormalized Free-Energy Principle | RFEP | The mathematical core principle | [10.5281/zenodo.19162705](https://doi.org/10.5281/zenodo.19162705) |
| Pattern | Pattern A: Functional Positivity under Gauge Constraint | Pattern A | The universal stability pattern | [10.5281/zenodo.19162705](https://doi.org/10.5281/zenodo.19162705) |

## Independent Foundations

| Name | Role | DOI |
|------|------|-----|
| RH Even Dominance (v2.1) | Independent RH proof, second route | [10.5281/zenodo.19546593](https://doi.org/10.5281/zenodo.19546593) |
| CRM Cosmology (I-V) | Independent dark energy model | [10.5281/zenodo.18728935](https://doi.org/10.5281/zenodo.18728935) |

These stand independently of FST. The RFEP was abstracted from them; they are not derived from it.

## Numerical Validation Scripts (Domain Supplements)

| Script | Paper | Description |
|--------|-------|-------------|
| `scripts/turbulence/compute_F_spectrum.py` | Turbulence | Verifies K41 as unique minimiser of F[E]; strict convexity test |
| `scripts/navier-stokes/compute_ds3_lorenz.py` | FST-Physics: Navier-Stokes | DS3 stress test on Lorenz attractor; TV saturation |
| `scripts/navier-stokes/compute_tll_ldi_lorenz.py` | NS-LDI | **Proof of Life**: TLL+LDI on Lorenz attractor (5/5 tests passed) |
| `scripts/dark-energy/compute_w_vs_desi.py` | Dark Energy | w_eff(z) comparison with DESI constraints |
| `scripts/dark-energy/compute_w_mapping.py` | Dark Energy | Correct w_eff → w_DE mapping + DESI grid scan |
| `scripts/dark-energy/compute_husawicki_mcmc.py` | Dark Energy | Hu-Sawicki f(R) MCMC fit against DESI+Planck+Cassini |
| `scripts/bsd/compute_height_saturation.py` | BSD | Height saturation test for quadratic twists |
| `scripts/hodge/compute_ghr_spectrum.py` | Hodge | GHR spectrum numerical verification |
| `scripts/hodge/compute_voisin_test.py` | Hodge | Voisin-style negative-control stress test |

## Repository Structure

```
functional-stability-theory/
├── masters/                    Three foundation papers
│   ├── zookeeper/              RH proof (microcluster closure)
│   ├── zeta-zoo/               Classification (SGE taxonomy)
│   └── spectrum-duality/       Physical axioms (RFEP, Pattern A)
├── fst-mathematics/            Domain supplements — Mathematics
│   ├── bsd/                    Rank-1 positivity (reformulation)
│   ├── hodge/                  No-go + easy direction
│   └── p-vs-np/                Witness entropy gap (reformulation)
├── fst-physics/                Domain supplements — Physics
│   ├── turbulence/             K41 spectrum (unconditional)
│   ├── yang-mills/             Mass gap (conditional)
│   └── navier-stokes/          Regularity + NS-LDI (conditional)
├── fst-cosmology/              Domain supplements — Cosmology
│   └── dark-energy/            CRM screening (validated)
├── fst-biology/                Domain supplements — Biology (in development)
├── fst-chemistry/              Domain supplements — Chemistry (planned)
└── scripts/                    Numerical validation
```

## Related Repositories

- [rh-even-dominance](https://github.com/research-line/rh-even-dominance) — Riemann Hypothesis: Even-dominance proof (foundation)
- [crm-cosmology](https://github.com/research-line/crm-cosmology) — Cooperative Renormalization Model (foundation)

## Author

Lukas Geiger -- ORCID: [0009-0005-7296-1534](https://orcid.org/0009-0005-7296-1534)

## License

[CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/)
