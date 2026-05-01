# RFEP Spectrum Duality

**Status:** v1.8 draft synchronized from the local RFEP master.
**Live Zenodo:** v1.7, concept DOI `10.5281/zenodo.19036190`.
**Next publication step:** upload the prepared v1.8 EN/DE PDFs as the next Zenodo version.

## Overview

This directory contains the public RFEP Spectrum Duality draft, the physics-master layer of the Functional Stability Theory hierarchy. RFEP translates stability and operator gains from the Zeta-Zoo core into a reusable physical normal form for Yang-Mills, Navier-Stokes, NS-LDI, turbulence, and cosmology.

The v1.8 draft incorporates the Post-Zookeeper transfer as a reusable pattern:

- rank-one or finite-rank defect language for physical obstructions
- mass-based cluster gaps instead of isolated eigenvalue checks
- coercive-complement control outside the target cluster
- Galerkin-faithfulness as an explicit proof obligation
- Prime-Hub boundary transfer as frontier or bouquet-style boundary data

## Public Files

| File | Purpose |
|---|---|
| `paper/FST_SPECTRUM_DUALITY_v0-8_en.tex` | English LaTeX source |
| `paper/FST_SPECTRUM_DUALITY_v0-8_en.pdf` | English compiled draft |
| `paper/FST_SPECTRUM_DUALITY_v0-8_de.tex` | German LaTeX source |
| `paper/FST_SPECTRUM_DUALITY_v0-8_de.pdf` | German compiled draft |

No private proof notes, review material, Zenodo credentials, or internal planning documents are published in this directory. Those remain local until the publication gate permits release.

## Relation to FST Hierarchy

```text
Level 0: RH even dominance and CRM reference systems
Level 1: Zookeeper and Zeta Zoo math masters
Level 2: RFEP Spectrum Duality physics master
Level 3: Domain supplements and applications
```

## Publication Gate

The repository is public, but the v1.8 Zenodo upload is still pending. Before uploading v1.8:

- verify that Zenodo metadata uses the concept DOI `10.5281/zenodo.19036190`
- set the language metadata for the bilingual EN/DE file set
- place the draft/preprint notice in the Zenodo description, not only in notes
- use concept DOIs for Zookeeper and Zeta Zoo related identifiers

## Reproducibility

This RFEP master is a theory manuscript and does not add a standalone numerical script bundle in this directory. Numerical Zookeeper diagnostics live under `masters/zookeeper/scripts/`; domain-level generated figures live under the repository-level `scripts/` tree.
