# Third-Party Licenses & Software Inventory

- **Repository:** `research-line/functional-stability-theory`
- **Version:** `1.0.6`
- **Audit Date:** `2026-09-26`
- **License Status:** `100% Permissive Open Source (0 AGPL, 0 Copyleft, 0 Cloud Telemetry)`
- **Umbrella:** [`open-bricks`](https://github.com/open-bricks) | **Parent Organization:** [`research-line`](https://github.com/research-line)
- **Primary Work License:** [Creative Commons Attribution 4.0 International (CC-BY-4.0)](https://creativecommons.org/licenses/by/4.0/)
- **Attribution & Notice:** Canonical author copyright and open-science ecosystem attribution are declared in [`NOTICE`](NOTICE).

---

## 1. Overview & Compliance Summary

`functional-stability-theory` is an open-science mathematical research programme and computational reproducibility repository under the `research-line` initiative and `open-bricks` federation. In strict compliance with open-science reproducibility standards and open-bricks governance policies, all direct runtime, optional scientific, and development-time dependencies undergo rigorous licensing audits.

- **Primary Repository License:** Creative Commons Attribution 4.0 International (`CC-BY-4.0`), granting unrestricted scholarly use, adaptation, sharing, and commercial application with appropriate attribution (see [`NOTICE`](NOTICE) and [`LICENSE`](LICENSE)).
- **Copyleft / Viral License Risk:** 0% (Zero AGPL, Zero GPL, Zero LGPL, Zero proprietary components).
- **Zero-Egress & Air-Gap Compliance:** 100% Offline (No network sockets, no phone-home telemetry, no external cloud API requests during execution).
- **Supply Chain Permissiveness:** All bundled, optional, and development dependencies are distributed under permissive open-source licenses (CC-BY-4.0, BSD-3-Clause, MIT, PSFL 2.0).

---

## 2. Direct Runtime Dependencies

The core verification ledgers, contract assertion suites, and metadata parsers execute entirely using standard Python without mandatory third-party runtime package dependencies:

| Component | Version Constraint | License | SPDX Identifier | Role in Project | Upstream Repository |
|---|---|---|---|---|---|
| **Python Standard Library** | `>=3.10` | Python Software Foundation License 2.0 | `PSF-2.0` | Core arithmetic, hashing, math/cmath operations, JSON/CSV ledger export, CLI parsing (`math`, `cmath`, `hashlib`, `json`, `csv`, `pathlib`, `typing`, `argparse`, `sys`, `os`, `re`, `subprocess`, `itertools`, `collections`) | `https://github.com/python/cpython` |

---

## 3. Optional Scientific & Diagnostic Dependencies

Specific numerical simulation scripts (e.g., in `scripts/`, `masters/atlas/scripts/`, `masters/zookeeper/scripts/`) utilize established scientific computing libraries for high-precision arithmetic, symbolic validation, and visualization:

| Library | License | SPDX Identifier | Purpose in Numerical Simulation | Upstream Repository |
|---|---|---|---|---|
| [`mpmath`](https://github.com/mpmath/mpmath) | BSD 3-Clause | `BSD-3-Clause` | Arbitrary-precision floating-point arithmetic for zeta zeros, Galerkin coefficients, and high-precision evaluation | `https://github.com/mpmath/mpmath` |
| [`sympy`](https://github.com/sympy/sympy) | BSD 3-Clause | `BSD-3-Clause` | Computer algebra and symbolic verification of algebraic identities and transfer relations | `https://github.com/sympy/sympy` |
| [`numpy`](https://github.com/numpy/numpy) | BSD 3-Clause | `BSD-3-Clause` | Vectorized numerical linear algebra, array operations, and grid computations | `https://github.com/numpy/numpy` |
| [`scipy`](https://github.com/scipy/scipy) | BSD 3-Clause | `BSD-3-Clause` | Numerical integration, ODE/PDE solving, optimization routines, and special functions | `https://github.com/scipy/scipy` |
| [`matplotlib`](https://github.com/matplotlib/matplotlib) | PSF-compatible (BSD-style) | `Matplotlib` | Publication-quality figure generation for spectral densities and phase diagrams | `https://github.com/matplotlib/matplotlib` |

*Note: All optional scientific dependencies are non-copyleft, permissive packages. Core contract verification and unit test suites execute cleanly without external scientific packages.*

---

## 4. Development, Linting & CI/CD Dependencies

Development tools are used exclusively during local testing, static analysis, packaging, and automated CI matrix verification. They are never bundled into distributed source snapshots or required by end-user research evaluators:

| Tool | Version Range | License | SPDX Identifier | Role in Development Workflow | Upstream Repository |
|---|---|---|---|---|---|
| [`pytest`](https://github.com/pytest-dev/pytest) | `>=8.0.0` | MIT License | `MIT` | Automated test runner, contract tests, and parameter sweep assertions | `https://github.com/pytest-dev/pytest` |
| [`ruff`](https://github.com/astral-sh/ruff) | `>=0.6.0` | MIT / Apache-2.0 | `MIT OR Apache-2.0` | High-speed static analysis, linting, and Python syntax hygiene | `https://github.com/astral-sh/ruff` |
| [`setuptools`](https://github.com/pypa/setuptools) | `>=61.0.0` | MIT License | `MIT` | PEP 517 / PEP 621 build backend and metadata specification | `https://github.com/pypa/setuptools` |

---

## 5. Governance & Research Invariant Compliance

The dependency and software architecture directly supports the 10 Governance and Research Invariants of Functional Stability Theory:

1. **100% Local-First & Zero-Egress (INV-LOCAL-01):** Numerical validation scripts run entirely on local computation cores; zero network sockets, external cloud calls, or phone-home telemetry.
2. **Unprivileged Mode / RunAsInvoker (INV-UNPRIV-02):** Execution strictly occurs within user-space privileges without administrative or root escalation.
3. **Deterministic Reproducibility (INV-DETERM-03):** Calculations produce bit-for-bit identical results and deterministic machine-readable CSV/JSON ledgers across all supported platforms.
4. **Claim-Level Disambiguation (INV-DISAMBIG-04):** Proven theorems, conditional transfer hypotheses, and open problems are explicitly separated in documentation and test assertions.
5. **Fail-Closed Ledger Gatekeeping (INV-FAILCLOSED-05):** Verification scripts immediately fail-closed and reject post-hoc envelopes, circular certificates, or degenerate parameters.
6. **Immutable Zenodo Anchors (INV-ZENODO-06):** Long-term scholarly provenance is secured through versioned Concept-DOIs deposited on CERN/Zenodo.
7. **Multi-OS Platform Parity (INV-PLATFORM-07):** Continuous verification ensures identical results on Windows, Linux, and macOS.
8. **Multi-Host & Cloud-Sync Hardening (INV-SYNCHARD-08):** Strict `.gitignore` rules prevent repository pollution from cloud-sync conflicts and multi-agent lock contention.
9. **Transparent Diagnostic Floor (INV-FLOOR-09):** Positive controls remain strictly separated from conditional transfer hypotheses.
10. **48h Security & 5-Day Triage SLA (INV-SLA-10):** Coordinated security disclosures through defined channels with guaranteed response times.

---

## 6. Authoritative License Texts

### Creative Commons Attribution 4.0 International (CC-BY-4.0)
```text
Creative Commons Corporation ("Creative Commons") is not a law firm and does
not provide legal services or legal advice. Distribution of Creative Commons
public licenses does not create a lawyer-client or other relationship.

Creative Commons makes its licenses and related information available on an
"as-is" basis. Creative Commons gives no warranties regarding its licenses, any
material licensed under their terms and conditions, or any related information.
Creative Commons disclaims all liability for damages resulting from their use
to the fullest extent possible.

Using Creative Commons Public Licenses
Creative Commons public licenses provide a standard set of terms and conditions
that creators and other rights holders may use to share original works of
authorship and other material subject to copyright and certain other rights
specified in the licensing terms.
```

### The 3-Clause BSD License
```text
Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

### The MIT License
```text
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

### Python Software Foundation License Version 2
```text
1. This LICENSE AGREEMENT is between the Python Software Foundation ("PSF"), and
the Individual or Organization ("Licensee") accessing and otherwise using this
software ("Python") in source or binary form and its associated documentation.

2. Subject to the terms and conditions of this License Agreement, PSF hereby
grants Licensee a nonexclusive, royalty-free, world-wide license to reproduce,
analyze, test, perform and/or display publicly, prepare derivative works, distribute,
and otherwise use Python alone or in any derivative version, provided, however, that
PSF's License Agreement and PSF's notice of copyright, i.e., "Copyright (c) 2001-2026
Python Software Foundation; All Rights Reserved" are retained in Python alone or
in any derivative version prepared by Licensee.
```
