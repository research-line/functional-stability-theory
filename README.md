# Functional Stability Theory (FST)

**Functional Stability Theory** is a unified mathematical programme that identifies a single structural challenge — *Functional Positivity under Gauge Constraint* (Pattern A) — as the common substrate of open problems in number theory, mathematical physics, and cosmology.

## The Three Masters

The programme rests on three foundation papers:

| Master | Title | Role | DOI |
|--------|-------|------|-----|
| [**Zookeeper**](masters/zookeeper/) | The Spectral Zookeeper | RH proof via CCM microcluster closure | [10.5281/zenodo.19673127](https://doi.org/10.5281/zenodo.19673127) |
| [**Zeta Zoo**](masters/zeta-zoo/) | The Zeta Zoo | Mathematical classification (SGE taxonomy, Boundary Theorem) | [10.5281/zenodo.19673227](https://doi.org/10.5281/zenodo.19673227) |
| [**Spectrum Duality**](masters/spectrum-duality/) | FST Spectrum Duality / RFEP | Physical instantiation (Pattern A, DS1-DS3) | [10.5281/zenodo.19162705](https://doi.org/10.5281/zenodo.19162705) |

## Domain Supplements

Domain-specific instantiations of Pattern A, each proving or reformulating a major open problem:

| Paper | Version | Status | Open Problem | DOI |
|-------|---------|--------|--------------|-----|
| [**Turbulence**](domain-proofs/turbulence/) | v1.3 | Unconditional | DFC1 empirical (only input) | [10.5281/zenodo.19056813](https://doi.org/10.5281/zenodo.19056813) |
| [**Dark Energy**](domain-proofs/dark-energy/) | v1.6 | Framework Note | Hu-Sawicki parameters open | [10.5281/zenodo.19036235](https://doi.org/10.5281/zenodo.19036235) |
| [**Yang-Mills**](domain-proofs/yang-mills/) | v2.1 | Conditional | Analytical proof of lambda < 0 | [10.5281/zenodo.19087433](https://doi.org/10.5281/zenodo.19087433) |
| [**Navier-Stokes**](domain-proofs/navier-stokes/) | v2.1 | Conditional | Assumption G2 (projection) | [10.5281/zenodo.19087449](https://doi.org/10.5281/zenodo.19087449) |
| [**NS Log-Distance**](domain-proofs/navier-stokes/) | v1.3 | Proof of Life | TLL for 3D NS open | [10.5281/zenodo.19056807](https://doi.org/10.5281/zenodo.19056807) |
| [**BSD**](domain-proofs/bsd/) | v1.1 | Reformulation | Higher Gross-Zagier (rank >= 2) | [10.5281/zenodo.19087443](https://doi.org/10.5281/zenodo.19087443) |
| [**Hodge**](domain-proofs/hodge/) | v1.1 | No-Go Theorem | = Deligne's question (1982) | [10.5281/zenodo.19087439](https://doi.org/10.5281/zenodo.19087439) |
| [**P vs NP**](domain-proofs/p-vs-np/) | v1.2 | Reformulation | Uniformity Bridge | [10.5281/zenodo.19056809](https://doi.org/10.5281/zenodo.19056809) |

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

## Architecture

```
         Zeta Zoo                    Spectrum Duality (RFEP)
    (Classification)                (Physical Instantiation)
           |                                  |
    SGE Taxonomy                      Pattern A + DS1-DS3
    Boundary Theorem                  Gauge Constraint
           |                                  |
           +----------------+-----------------+
                            |
                     Domain Supplements
                            |
        +--------+--------+--------+--------+
        |        |        |        |        |
       TU      YM       NS      DE     BSD/Hodge/PNP
   (proven) (cond.)  (cond.)  (valid.)  (reformul.)
```

The **Zookeeper** supplies the solved Riemann Hypothesis via the CCM route (independent of Pattern A).
The **Zeta Zoo** classifies all zeta-type families and provides the mathematical framework.
**Spectrum Duality** provides the physical axiom package from which domain supplements derive.

## Independent Foundations

| Name | Role | DOI |
|------|------|-----|
| RH Even Dominance (v2.1) | Independent RH proof, second route | [10.5281/zenodo.19546593](https://doi.org/10.5281/zenodo.19546593) |
| CRM Cosmology (I-V) | Independent dark energy model | [10.5281/zenodo.18728935](https://doi.org/10.5281/zenodo.18728935) |

These stand independently of FST. The RFEP was abstracted from them; they are not derived from it.

## Author

Lukas Geiger -- ORCID: [0009-0005-7296-1534](https://orcid.org/0009-0005-7296-1534)

## License

[CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/)
