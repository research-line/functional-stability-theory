---
name: functional-stability-theory
version: "1.0.6"
type: research-toolkit
standalone: true
visibility: public
provenance: "research-line / open-bricks"
description: "Unified mathematical programme identifying Functional Positivity under Gauge Constraint (Pattern A) as common substrate across number theory, mathematical physics, and cosmology. Provides deterministic, zero-egress numerical verification ledgers."
---

# Functional Stability Theory (FST) Skill

Use this skill when interacting with the `research-line/functional-stability-theory` research software, conducting scholarly literature analysis, auditing mathematical proofs across the Five Masters, or executing deterministic numerical verification scripts.

## Core Conceptual Model

Functional Stability Theory establishes that open problems across distinct mathematical and physical domains share a common mathematical substrate: **Functional Positivity under Gauge Constraint** (Pattern A), governed by the **Renormalized Free-Energy Principle (RFEP)**.

### The Five Masters Foundation Hierarchy

1. **Master I (Zookeeper / Analytic Number Theory):** Shift-parity symmetry, spectral determinant decomposition, Weil quadratic-form positivity, and SGE (Spectral Gap Equilibrium).
2. **Master II (Soliton & Integrability):** Inverse scattering transform, Lax pairs, finite-band potentials, and non-perturbative stability.
3. **Master III (Statistical Mechanics & Criticality):** RG flows, partition function zeros, conformal invariance, and thermodynamic limit bounds.
4. **Master IV (Differential Geometry & Gauge Theory):** Curvature bounds, Yang-Mills Casimir corridor reduction, orbit rigidity, and geometric analysis.
5. **Master V (Metric & Measure Geometry):** Optimal transport, curvature-dimension conditions, Wasserstein geometry, and measure concentration.

## Agent Execution & Verification Runbooks

All numerical scripts in this repository follow strict research integrity invariants (`INV-LOCAL-01` to `INV-SLA-10`):
- **100% Offline / Zero-Egress (`INV-NET-03`):** Scripts run strictly offline; no network sockets or external telemetry.
- **Deterministic Bit-for-Bit Output (`INV-DET-04`):** Fixed pseudo-random seeds and stable numeric representations ensure identical output.
- **Unprivileged Execution (`INV-PRIV-05`):** Standard `RunAsInvoker` user mode without administrative or elevated rights.

### 1. Running Automated Quality Gates & Contract Tests

To execute the full test suite and verify contract invariants:

```bash
# Execute full test suite (pytest)
pytest -ra -v

# Run metadata and governance contract tests specifically
pytest tests/test_metadata.py -v
```

### 2. Domain Validation Scripts

Standalone numerical scripts reside in `scripts/<domain>/`:

- **Yang-Mills Continuum Gate Audit:**
  ```bash
  python scripts/yang-mills/compute_continuum_gate_audit.py
  ```
- **Yang-Mills Casimir Corridor & Orbit Rigidity:**
  ```bash
  python scripts/yang-mills/compute_casimir_corridor_ledger.py
  python scripts/yang-mills/compute_orbit_rigidity_ledger.py
  ```
- **Navier-Stokes BV Selection & Attractors:**
  ```bash
  python scripts/navier-stokes/compute_bv_selection.py
  python scripts/navier-stokes/compute_bv_multi_attractor.py
  ```
- **Turbulence Shell DFC Waterline Ledger:**
  ```bash
  python scripts/turbulence/compute_shell_dfc_waterline_ledger.py
  ```
- **Birch and Swinnerton-Dyer (BSD) Regulator Verification:**
  ```bash
  python scripts/bsd/compute_bsd_verification.py
  ```

## Reference Documentation & Machine-Readable Context

- **LLM Context & Search Index:** [`llms.txt`](llms.txt) and [`llm.txt`](llm.txt)
- **Contributing & Governance:** [`CONTRIBUTING.md`](CONTRIBUTING.md)
- **Code of Conduct:** [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md)
- **Security & Vulnerability SLA:** [`SECURITY.md`](SECURITY.md)
- **Supply Chain & Plain-Text SBOM:** [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md) and [`THIRD_PARTY_LICENSES.txt`](THIRD_PARTY_LICENSES.txt)
- **Legal Attribution & Disclaimer:** [`NOTICE`](NOTICE) and [§ 521 BGB Statutory Notice](README.md#security-policy--statutory-notice)
