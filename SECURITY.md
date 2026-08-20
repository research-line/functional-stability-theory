# Security & Research Integrity Policy / Sicherheits- und Integritätsrichtlinie

English | [Deutsch](#deutsche-fassung)

---

## English

### 1. Scope & Research Integrity

`functional-stability-theory` is an open-science theoretical and mathematical research repository. It provides analytical papers, structural proofs, LaTeX sources, and reproducible Python numerical validation scripts.

Security in this repository centers on:
1. **Mathematical & Research Integrity**: Ensuring numerical diagnostics and simulation scripts produce verifiable, deterministic, and untampered results.
2. **Zero-Egress & Offline Execution**: All numerical scripts (`scripts/`, `masters/*/scripts/`) run 100% locally with zero external network calls, telemetry, or remote telemetry egress.
3. **Deterministic Dependencies**: Minimal and well-audited Python dependencies (`numpy`, `scipy`, `matplotlib`, `pytest`, `ruff`), with no unvetted dynamic code loading.

---

### 2. Threat Model & Invariants

| Dimension | Policy / Invariant |
|---|---|
| **Network Egress** | Zero network egress during script execution and test runs. |
| **Execution Environment** | Unprivileged local execution; never requires administrator or root privileges. |
| **Data Privacy** | No personally identifiable information (PII) or proprietary telemetry is collected or stored. |
| **Reproducibility** | All numerical validation scripts are self-contained and reproducible. |
| **Provenance** | Key paper milestones and artifacts are anchored to immutable Zenodo Concept-DOIs. |

---

### 3. Vulnerability & Anomaly Reporting

If you identify a security vulnerability, anomalous script behavior, or supply-chain issue:

- **Email**: Contact [`security@ellmos.ai`](mailto:security@ellmos.ai) or [`support@lukasgeiger.com`](mailto:support@lukasgeiger.com).
- **GitHub Security Advisory**: Open a private advisory via [GitHub Security Advisories](https://github.com/research-line/functional-stability-theory/security/advisories/new).

Please do not open public issues for sensitive security or integrity vulnerabilities before coordination.

---

<a name="deutsche-fassung"></a>
## Deutsche Fassung

### 1. Geltungsbereich & Forschungsintegrität

`functional-stability-theory` ist ein Open-Science-Repository für theoretische und mathematische Forschung. Es umfasst analytische Arbeiten, mathematische Beweisführungen, LaTeX-Quellen und reproduzierbare Python-Validierungsskripte.

Sicherheit und Integrität umfassen in diesem Kontext:
1. **Mathematische & Wissenschaftliche Integrität**: Gewährleistung, dass numerische Diagnostiken und Simulationsskripte deterministische, überprüfbare und unverfälschte Ergebnisse liefern.
2. **Zero-Egress & Offline-Ausführung**: Alle numerischen Skripte (`scripts/`, `masters/*/scripts/`) laufen zu 100% lokal ohne externe Netzwerkaufrufe, Tracking oder Telemetrie.
3. **Minimale Abhängigkeiten**: Schlanke, geprüfte Python-Standard- und Wissenschaftsbibliotheken (`numpy`, `scipy`, `matplotlib`, `pytest`, `ruff`) ohne dynamisches Nachladen von Fremdcode.

---

### 2. Bedrohungsmodell & Invarianten

| Dimension | Richtlinie / Invariante |
|---|---|
| **Netzwerk-Egress** | 0% Netzwerk-Egress während der Skriptausführung und Testläufe. |
| **Ausführungskontext** | Standard-Benutzerkontext; benötigt niemals Administrator- oder Root-Rechte. |
| **Datenschutz** | Keine Erfassung oder Übertragung von personenbezogenen Daten oder Telemetrie. |
| **Reproduzierbarkeit** | Alle Validierungsskripte sind in sich geschlossen und deterministisch reproduzierbar. |
| **Provenienz** | Wichtige Meilensteine und Arbeiten sind über unveränderliche Zenodo-Konzept-DOIs verankert. |

---

### 3. Meldung von Schwachstellen & Anomalien

Wenn Sie eine Sicherheitslücke, unerwartetes Skriptverhalten oder Integritätsprobleme entdecken:

- **E-Mail**: Senden Sie einen Bericht an [`security@ellmos.ai`](mailto:security@ellmos.ai) oder [`support@lukasgeiger.com`](mailto:support@lukasgeiger.com).
- **GitHub Security Advisory**: Erstellen Sie einen vertraulichen Hinweis via [GitHub Security Advisories](https://github.com/research-line/functional-stability-theory/security/advisories/new).

Bitte erstellen Sie keine öffentlichen GitHub-Issues für sicherheitsrelevante Schwachstellen vor Abschluss der Prüfung.
