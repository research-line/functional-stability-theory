# FST Applications

Application-scale companions of the **Functional Stability Theory** programme, structured by physical scale (Level 1 = particles, 1b = chemistry, 1c = biology, Level 2 = cosmology).

| Slot | Scale | Title | Concept-DOI | Latest |
|---|---|---|---|---|
| **FST-I** | Particles / Thermodynamics | Game-Theoretic Framework for the Thermodynamic Stability of Fundamental Parameters | [10.5281/zenodo.20130544](https://doi.org/10.5281/zenodo.20130544) | v1.6 live |
| **FST-II** | Chemistry / Autocatalysis | Chemical Stability and Autocatalytic Selection | [10.5281/zenodo.20130563](https://doi.org/10.5281/zenodo.20130563) | v1.3 live / local post-v1.3 guardrail |
| **FST-III** | Biology / Protein folding | Biological Stability and Nash Frustration | [10.5281/zenodo.20130573](https://doi.org/10.5281/zenodo.20130573) | v1.2 live / local v1.3 guardrail |
| **FST-IV** | Cosmology (collector slot) | Cosmological Stability — Collector for CRM, Saturation, QG-CRM and Dark Energy | *(in preparation, v0 skeleton)* | v0.2 (local) |
| ↳ FST-DE | sub: Dark Energy | Dark Energy as Residual Vacuum Free Energy | [10.5281/zenodo.19036235](https://doi.org/10.5281/zenodo.19036235) | v1.9 |
| **Hub** | Programme umbrella | FST — A Programmatic Hub: Universal Convexity-Uniqueness Across Scales | [10.5281/zenodo.20130499](https://doi.org/10.5281/zenodo.20130499) | v1.4 |

All papers are published on Zenodo with persistent DOIs (except the FST-IV Cosmology collector which is currently a v0 skeleton and not yet on Zenodo). This repository holds the LaTeX sources and compiled PDFs.

The **FST-DE** paper lives under `../fst-cosmology/dark-energy/` (cosmology domain) and is conceptually part of two programme lines:
- the FST application series (sub-component of the FST-IV Cosmology slot), and
- the CRM follow-up series (uses the R + γR² Lagrangian; CRM-V provides the universal tanh normal form; CRM-VI provides the UV completion).

## Honest status

- **Hub:** v1.4, programmatic survey with FST Series + FST Foundations tables.
- **FST-I:** v1.6 is live on Zenodo ([record DOI 10.5281/zenodo.21386047](https://doi.org/10.5281/zenodo.21386047)); the public repo package mirrors the July 2026 alpha-tolerance guardrail paper set. No further Zenodo upload is open without a new file or metadata change.
- **FST-II:** v1.3 is live on Zenodo; the public repo mirrors the local post-v1.3 proof/paper-math plus citation guardrail package, while an optional v1.4 maintenance upload remains a separate decision.
- **FST-III:** v1.2 is live on Zenodo; the repo contains the tightened local v1.3 guardrail package, while the protein-folding code/data expansion continues in the private companion workflow.
- **FST-IV Cosmology collector:** KONZEPT + v0.2 skeleton; LaTeX body to be filled with the CRM-series / FST-DE / CRM-VI synthesis and the Flow-Zeta math reformulation.

## Scripts &amp; results

- FST-I public artefacts are mirrored under `fst-i-thermodynamic/papers/`, `fst-i-thermodynamic/scripts/` and `fst-i-thermodynamic/results/`, including the entropy scan outputs and the 2D log-Hessian pilot JSON. A dedicated scripted reproducer for the Hessian pilot is still future work.
- FST-II reproducibility artefacts are mirrored under `fst-ii-chemical/scripts/` and `fst-ii-chemical/results/`.
- FST-III public artefacts are mirrored under `fst-iii-biological/papers/`, `scripts/fst_iii/` and `results/fst_iii/`. Curated JSON/plot outputs are public; raw per-protein work directories and server logs remain local-only. The chaperone game-theory companion paper **FST-Nash** is published separately ([DOI: 10.5281/zenodo.20402751](https://doi.org/10.5281/zenodo.20402751)) with its own public repository [`research-line/fst-nash`](https://github.com/research-line/fst-nash).
- FST-IV computational artefacts are still primarily tracked locally and will be mirrored as the collector paper stabilizes.
