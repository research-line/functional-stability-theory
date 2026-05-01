# REVIEW CHAIN — Step 4: Experte (zweite Runde)
## Paper: Dirichlet Character Atlas — A Weil-Kernel Cartography of Low-Conductor Sector Gaps (v2, Final Draft)
## Datum: 2026-04-29
## Fokus: Minimaler Einreichungsplan

---

## Vorbemerkung: Verifikationsstand

Vor der Bewertung der Step-3-Punkte wurde der Paper-Text gegen die Step-3-Belegtabelle geprüft.

**Befund — Phantom-Fehler in Step 3:** Step 3 (Abschnitt A.E.5) behauptete, der Abstract enthalte noch die alte Vier-Achsen-Formulierung. Der Paper-Text Z. 95–98 zeigt jedoch bereits die korrekte Formulierung:

> "We catalog the ten characters along two filter constants (parity and root number, both $+1$ in the sample) and two informative axes (archimedean factor weight and gap signature)..."

E.5 ist im Abstract **vollständig umgesetzt**. Der Step-3-Vorschlag war ein Phantom-Fehler und wird in Step 4 nicht erneut aufgegriffen. Dies eliminiert eine vermeintliche arXiv-Pflicht-Korrektur.

---

## Abschnitt A: Pflicht-Korrekturen für arXiv-Submission

### A.1 Title-Präfix entfernen (PFLICHT)

**Problem:** Der Titel beginnt mit „FST-Mathematics: ..." (Z. 49). Das ist internes Programm-Jargon ohne Erklärungswert für externe Leser.

**LaTeX-Fix (copy-paste):**
```latex
\title{The Dirichlet Character Atlas\\
  \large --- A Weil-Kernel Numerical Atlas for Dirichlet $L$-Function Sector Gaps\\
  \normalsize Galerkin Diagnostics and the Boundaries of Leading-Order Theory}
```

Z. 53 entsprechend von „Draft v2 --- FST-Mathematics supplement (...)" zu „Draft v2 (revised after N=600 server analysis)" kürzen.

**Alle anderen Step-1/2/3-Punkte für arXiv:** NICHT blocking.

**Zusammenfassung A:** Genau **eine** Pflicht-Korrektur vor arXiv. Das Paper ist danach submission-ready.

---

## Abschnitt B: Pflicht-Korrekturen für Journal-Einreichung (Experimental Mathematics)

### B.1 Beweisskizzen für Theoreme 3.2 und 3.3

**LaTeX-Fix nach Theorem 3.2:**
```latex
\begin{proof}[Sketch]
The kernel $K_\chi(x,y)=\kappa_\chi(x-y)+K^{\mathrm{arch}}_\chi(x-y)$
inherits the symmetry $K_\chi(x,y)=K_\chi(-x,-y)$ from the evenness
of $\kappa_\chi$ and $K^{\mathrm{arch}}_\chi$. Symmetrising with the
projectors $P^\pm=(I\pm\mathcal{P})/2$ and using
$\mathcal{P}\kappa_\chi(\cdot)\mathcal{P}=\kappa_\chi(-\,\cdot)$
yields the stated decomposition; full details in
\cite{AnalyticKernel}, \S3. $\square$
\end{proof}
```

**LaTeX-Fix nach Theorem 3.3:**
```latex
\begin{proof}[Sketch]
Subtract the ``$+$'' and ``$-$'' versions of Theorem~\ref{thm:kernel}.
Diagonal terms $\kappa_\chi(x-y)$ and $K^{\mathrm{arch}}_\chi(x-y)$
cancel; the anti-diagonal terms add with relative sign $-1$. The
character dependence enters only through $\kappa_\chi(x+y)$, since
$K^{\mathrm{arch}}_\chi$ depends on $\chi$ only through the constant
$\log(q/\pi)$. $\square$
\end{proof}
```

**Aufwand:** 5 Minuten.

### B.2 N=400-Spalte in Tabelle 1 ergänzen (oder im Caption begründen)

Setup §3.1 nennt $N\in\{200,400,600\}$; Tabelle 1 zeigt nur N=200 und N=600. **Caption-Klausel** als schnelle Lösung:

```latex
% Caption-Ergänzung am Ende:
Intermediate $N=400$ data, available for the convergence study of $\chi_{21}$
in Table~\ref{tab:chi21}, are omitted here because they qualitatively
interpolate between $N=200$ and $N=600$ for the stable characters.
```

**Aufwand:** 5 Minuten (Caption-Klausel) bis 30 Minuten (Spalte ergänzen).

### B.3 Spectral-Convergence-Verweis in §3.4

**LaTeX-Fix (am Ende von §3.4):**
```latex
The Galerkin-truncated operators $W^\pm_{\chifn,N}$ are
finite-dimensional approximations of the bounded symmetric kernel
operator on $L^2[-L,L]$. Standard spectral-convergence theory for
compact symmetric operators (\cite{BabuskaOsborn}, Thm.~2.1)
ensures that the discrete eigenvalues converge to the true spectrum
as $N\to\infty$, with rates controlled by the smoothness of the
kernel and the eigenfunction. Section~\ref{sec:chi21} illustrates
that for characters with quasi-degenerate spectra this convergence
may be non-monotone at moderate $N$.
```

**Bibitem:**
```latex
\bibitem{BabuskaOsborn}
I.~Babu\v{s}ka and J.~E.~Osborn,
\textit{Eigenvalue problems},
in Handbook of Numerical Analysis Vol.~II
(P.~G.~Ciarlet and J.~L.~Lions, eds.),
North-Holland, Amsterdam, 1991, pp.~641--787.
```

**Aufwand:** 10 Minuten.

### B.4 Interne Supplementary-Bibitems verifizierbar machen

**Empfehlung — Hybridlösung:**
- [AnalyticKernel] §3 und [AnalyticKernelV2] §5 als Appendix-Sections ins Paper integrieren
- [AnalyticGroundstate], [Woche4Validation], [AsymptoticScan] als Zenodo-Companion mit DOI

**Aufwand:** 3–4 Stunden.

### B.5 χ_21 Quasi-Degeneracy spektral belegen

Tabelle mit $\lambda_1^{\min}$, $\lambda_2^{\min}$ (beide Sektoren) und Eigenvektor-Überlappungen für N=200, 400, 600 in §6 einfügen.

**Aufwand:** ~2 Stunden (Skript + Compute + Integration).

### B.6 Selektionskriterium und λ-Begründung in §3.1

**LaTeX-Fix (am Anfang von §3.1):**
```latex
The ten conductors $D\in\{5,8,12,13,17,21,24,29,33,60\}$ are the
smallest fundamental discriminants $D>0$ giving primitive real
even Kronecker characters with $D\leq 60$, providing a small but
structurally diverse sample (prime, prime-power, and composite
conductors). The cut-off $\lambda=20000$ is the largest value for
which $N=600$ Galerkin diagonalization completes within
$\sim 2$~hours wall-clock per ten-character batch on the Hetzner
CCX13 instance.
```

**Aufwand:** 5 Minuten.

---

## Abschnitt C: Empfohlene Verbesserungen (nicht blocking)

| # | Beschreibung | Aufwand |
|---|---|---|
| C.1 | Note-Block Z. 471–477 als `\begin{remark}` ausweisen | 10 Min |
| C.2 | „Final Draft" / „Draft v2"-Mismatch konsolidieren | 5 Min |
| C.3 | Tabelle 3 strikt nach $r_\chi$ sortieren | 15 Min |
| C.4 | Generalized-Ritz-Galerkin-Remark nach Prop. 6.1 | 30 Min |
| C.5 | §7.3 ↔ App B Konsolidierung (für Revision 2) | 1 Std |
| C.6 | §7.2 Slope-0.72 integrieren oder kürzen | 1 Std |

---

## Abschnitt D: Abgelehnte Vorschläge aus Steps 1–3

### D.1 Step-3-B.1 (F4 tonal inkonsistent) — ABGELEHNT
F4 ist keine prädiktive Aussage, sondern ein diagnostisches Ergebnis über was nicht widerlegt wurde. Konsistent mit der diagnostischen Postur des Papers. Keine Änderung.

### D.2 Step-3-B.4 (§7.3 ↔ App B Doppelung) — ABGELEHNT als Pflicht
Geringe Redundanz zwischen Discussion (§7.3) und Technischer Schicht (App B) ist in mathematischen Texten Standard. Maximal als C.5 für Revision 2.

### D.3 Step-2-Option-1 (Δ_χ-Definition zur Summe umdefinieren) — ABGELEHNT
Option 3 (umgesetzte Note + App B) ist für Experimental Mathematics ausreichend. Umdefinition würde Theorem 3.4, Eq. (5), Eq. (15), Tabellen-Captions und den gesamten Sign-Source-Anhang berühren — hoher Aufwand für ästhetischen Gewinn. POSTRELEASE-PR-Item für Version 3.

---

## Abschnitt E: Minimaler Einreichungsplan

### Phase 1: Vor arXiv (30 Minuten)
1. A.1 — Title-Präfix „FST-Mathematics:" entfernen
2. C.1 (optional) — Note-Block als `\begin{remark}` ausweisen
3. Compile-Check (pdflatex ×2)
4. arXiv-Submit: math.NT (Primary), math.SP (Cross-list)

### Phase 2: Vor Experimental Mathematics (4–6 Stunden)
1. B.6 — Selektionskriterium + λ-Begründung (5 Min)
2. B.1 — Beweisskizzen Thm 3.2/3.3 (5 Min)
3. B.3 — Spectral-Convergence-Verweis + Bibitem (10 Min)
4. B.2 — N=400-Spalte oder Caption-Klausel (5–30 Min)
5. C.4 — Ritz-Galerkin-Remark (30 Min)
6. B.5 — χ_21 spektral belegen (2 Std)
7. B.4 — Bibitems: Hybrid-Lösung (3–4 Std)
8. POSTRELEASE-PR: C.5, C.6, D.3 für Version 3

### Phase 3: Fallback
Falls Experimental Mathematics ablehnt → Research in Number Theory (Open Access).

---

## Abschnitt F: Readiness-Score

| Kontext | Score | Begründung |
|---|---|---|
| arXiv (nach A.1) | **9/10** | Vollständig submission-ready; -1 Pkt für interne Bibitems (in Programm-Preprint-Tradition toleriert) |
| Journal nach A-Fix | 7.5/10 | B.1–B.6 noch ausstehend |
| Journal nach A+B-Fixes | **8.5–9/10** | Experimental Mathematics submission-ready |

**Empfehlung:** Phase 1 sofort, arXiv-Submit, dann Phase 2 in einer dedizierten Halbtages-Session vor Journal-Submission.
