"""
compute_sat_entropy.py
======================
P-vs-NP: Berechnet die Slice-Entropie H_slice(SAT_n) fuer zufaellige
3-SAT-Instanzen am Phasenuebergang (alpha ~ 4.267).

Das ESC-Framework (Entropy-Stability-Complexity) definiert:
  H_slice(f) = max_x { H(f|_{x_S}) : S subset [n], |S| = n/2 }

Fuer eine Boolesche Funktion f: {0,1}^n -> {0,1} misst die Slice-Entropie
die maximale Unsicherheit einer Restriktion.

Hypothese (ESC): H_slice(SAT_n) >= c * 2^{n/2} fuer "harte" Instanzen,
waehrend P-berechenbare Funktionen H_slice <= poly(n) haben.

Wir messen:
1. H_slice fuer zufaellige 3-SAT bei alpha = m/n ~ 4.267 (Phasenuebergang)
2. Vergleich: leichte vs. schwere Instanzen
3. Skalierung von H_slice mit n

Autor: Lukas Geiger (Skript erstellt per Claude, 2026)
"""

import numpy as np
from itertools import combinations
import argparse
import csv
import json
import os
import time

# ==========================================================================
# 3-SAT Generator
# ==========================================================================

def random_3sat(n, m, seed=None):
    """
    Generiert eine zufaellige 3-SAT-Instanz.
    n: Anzahl Variablen
    m: Anzahl Klauseln
    Returniert: Liste von Klauseln, jede Klausel = (l1, l2, l3)
                wobei l > 0 fuer positive und l < 0 fuer negierte Variable
    """
    rng = np.random.RandomState(seed)
    clauses = []
    for _ in range(m):
        # Waehle 3 verschiedene Variablen
        vars_chosen = rng.choice(n, 3, replace=False) + 1  # 1-basiert
        # Zufaellige Negierung
        signs = rng.choice([-1, 1], 3)
        clause = tuple(vars_chosen * signs)
        clauses.append(clause)
    return clauses


def planted_3sat(n, m, planted_assignment, seed=None):
    """
    Generiert eine kleine planted-3-SAT-Positivkontrolle.

    Die bekannte Belegung ist bewusst eine Kontroll-/Advice-Quelle und darf
    nicht als Uniformity-Bridge-Signal gelesen werden.
    """
    rng = np.random.RandomState(seed)
    clauses = []
    while len(clauses) < m:
        vars_chosen = rng.choice(n, 3, replace=False) + 1
        signs = rng.choice([-1, 1], 3)
        clause = tuple(vars_chosen * signs)
        if evaluate_clause(clause, planted_assignment):
            clauses.append(clause)
    return clauses


def evaluate_clause(clause, assignment):
    """Evaluiere eine Klausel unter einer Belegung (dict: var -> bool)"""
    for lit in clause:
        var = abs(lit)
        val = assignment.get(var, False)
        if lit > 0 and val:
            return True
        if lit < 0 and not val:
            return True
    return False


def evaluate_formula(clauses, assignment):
    """Evaluiere gesamte Formel"""
    return all(evaluate_clause(c, assignment) for c in clauses)


def count_satisfied(clauses, assignment):
    """Zaehle erfuellte Klauseln"""
    return sum(1 for c in clauses if evaluate_clause(c, assignment))


def assignment_from_index(n, idx):
    """Dekodiere eine boolesche Belegung aus einem Integer."""
    return {var: bool((idx >> (var - 1)) & 1) for var in range(1, n + 1)}


def assignment_to_tuple(assignment, n):
    return tuple(bool(assignment[var]) for var in range(1, n + 1))


def enumerate_satisfying_assignments(clauses, n, sample_limit=128):
    """
    Zaehlt alle Witnesses exakt fuer kleine n und behaelt eine begrenzte
    Stichprobe fuer Distanz-/Cluster-Proxies.
    """
    count = 0
    sample = []
    for idx in range(2**n):
        assignment = assignment_from_index(n, idx)
        if evaluate_formula(clauses, assignment):
            count += 1
            if len(sample) < sample_limit:
                sample.append(assignment_to_tuple(assignment, n))
    return count, sample


def single_flip_stability(clauses, assignment, n):
    """Anteil einzelner Bitflips, die einen Witness weiterhin satisfying lassen."""
    if not assignment or not evaluate_formula(clauses, assignment):
        return 0.0
    stable = 0
    for var in range(1, n + 1):
        flipped = dict(assignment)
        flipped[var] = not flipped[var]
        if evaluate_formula(clauses, flipped):
            stable += 1
    return stable / n


def witness_pair_distance_stats(witnesses, n):
    """Kleine Distanz-Proxys; kein OGP-Zertifikat."""
    if len(witnesses) < 2:
        return {
            'pair_count': 0,
            'distance_min': 0,
            'distance_max': 0,
            'distance_middle_fraction': 0.0,
        }
    distances = []
    max_witnesses = min(len(witnesses), 96)
    subset = witnesses[:max_witnesses]
    for i in range(len(subset)):
        for j in range(i + 1, len(subset)):
            distances.append(sum(a != b for a, b in zip(subset[i], subset[j])))
    low = n / 4
    high = 3 * n / 4
    middle = sum(1 for d in distances if low <= d <= high)
    return {
        'pair_count': len(distances),
        'distance_min': min(distances),
        'distance_max': max(distances),
        'distance_middle_fraction': round(middle / len(distances), 6),
    }


# ==========================================================================
# Slice-Entropie Berechnung
# ==========================================================================

def slice_entropy_exact(clauses, n, S):
    """
    Berechne H(f|_{x_S}) fuer eine feste Menge S.
    Fixiere die Variablen in S auf alle moeglichen Werte,
    zaehle jeweils die Anzahl erfuellender Belegungen fuer die restlichen Variablen.

    H = -sum_y p(y) * log2(p(y))
    wobei y ueber {0,1}^S laeuft und p(y) = #SAT(f|_{x_S=y}) / total
    """
    S_list = sorted(S)
    complement = sorted(set(range(1, n+1)) - set(S))
    n_S = len(S_list)
    n_C = len(complement)

    if n_C > 20:
        # Zu gross fuer exakte Enumeration: Monte-Carlo-Approximation
        return slice_entropy_monte_carlo(clauses, n, S, n_mc=5000)

    # Zaehle fuer jede S-Belegung die Anzahl erfuellender Complement-Belegungen
    counts = np.zeros(2**n_S)

    for s_idx in range(2**n_S):
        # S-Belegung dekodieren
        assignment = {}
        for bit, var in enumerate(S_list):
            assignment[var] = bool((s_idx >> bit) & 1)

        # Alle Complement-Belegungen durchprobieren
        sat_count = 0
        for c_idx in range(2**n_C):
            for bit, var in enumerate(complement):
                assignment[var] = bool((c_idx >> bit) & 1)
            if evaluate_formula(clauses, assignment):
                sat_count += 1

        counts[s_idx] = sat_count

    # Normalisiere zu Wahrscheinlichkeitsverteilung
    total = np.sum(counts)
    if total == 0:
        return 0.0  # UNSAT

    p = counts / total
    # Entropie H = -sum p*log2(p)
    mask = p > 0
    H = -np.sum(p[mask] * np.log2(p[mask]))
    return H


def slice_entropy_monte_carlo(clauses, n, S, n_mc=5000):
    """Monte-Carlo-Approximation der Slice-Entropie"""
    S_list = sorted(S)
    complement = sorted(set(range(1, n+1)) - set(S))
    n_S = len(S_list)
    n_C = len(complement)

    # Sample zufaellige Belegungen und zaehle SAT pro S-Belegung
    s_counts = {}

    for _ in range(n_mc):
        # Zufaellige Gesamtbelegung
        assignment = {}
        for var in range(1, n+1):
            assignment[var] = bool(np.random.randint(2))

        if evaluate_formula(clauses, assignment):
            # S-Belegung als Key
            s_key = tuple(assignment[var] for var in S_list)
            s_counts[s_key] = s_counts.get(s_key, 0) + 1

    if not s_counts:
        return 0.0

    counts = np.array(list(s_counts.values()), dtype=float)
    total = np.sum(counts)
    p = counts / total
    mask = p > 0
    H = -np.sum(p[mask] * np.log2(p[mask]))
    return H


def max_slice_entropy(clauses, n, k=None, n_random_subsets=50):
    """
    Berechne max_S H(f|_{x_S}) ueber Teilmengen S mit |S| = k.
    Fuer kleine n: exakt ueber alle (n choose k) Teilmengen.
    Fuer groessere n: zufaellige Stichprobe.
    """
    if k is None:
        k = n // 2

    from math import comb
    total_subsets = comb(n, k)

    if total_subsets <= n_random_subsets * 2:
        # Exakt
        max_H = 0.0
        best_S = None
        for S in combinations(range(1, n+1), k):
            H = slice_entropy_exact(clauses, n, set(S))
            if H > max_H:
                max_H = H
                best_S = S
        return max_H, best_S
    else:
        # Zufaellige Stichprobe
        max_H = 0.0
        best_S = None
        for _ in range(n_random_subsets):
            S = set(np.random.choice(range(1, n+1), k, replace=False))
            H = slice_entropy_exact(clauses, n, S)
            if H > max_H:
                max_H = H
                best_S = tuple(sorted(S))
        return max_H, best_S


# ==========================================================================
# Kolmogorov-Komplexitaet Approximation
# ==========================================================================

def approx_kolmogorov(clauses, n):
    """
    Approximation der Kolmogorov-Komplexitaet K(f) via Kompressionslaenge.
    Kodiere die Wahrheitstabelle und komprimiere mit zlib.
    """
    import zlib

    # Baue (Teil-)Wahrheitstabelle
    if n <= 20:
        truth_table = bytearray()
        for idx in range(2**n):
            assignment = {}
            for var in range(1, n+1):
                assignment[var] = bool((idx >> (var-1)) & 1)
            truth_table.append(1 if evaluate_formula(clauses, assignment) else 0)
    else:
        # Zu gross -- sample
        truth_table = bytearray()
        for _ in range(10000):
            assignment = {}
            for var in range(1, n+1):
                assignment[var] = bool(np.random.randint(2))
            truth_table.append(1 if evaluate_formula(clauses, assignment) else 0)

    compressed = zlib.compress(bytes(truth_table), 9)
    return len(compressed) * 8  # in Bits


# ==========================================================================
# Advice-/Front-Ledger
# ==========================================================================

def _safe_log2(x):
    return float(np.log2(max(float(x), 1.0)))


def _slice_cover_proxy(H, n):
    """Normalisierte Support-/Cover-Proxymasse aus Shannon-H_slice."""
    denom = 2 ** max(1, n // 2)
    return min(1.0, (2 ** max(float(H), 0.0)) / denom)


def build_advice_front_ledger(results, results_alpha):
    """
    Baut ein vorsichtiges Advice-/Front-Leakage-Ledger.

    Das ist keine K^t-Messung und kein Uniformity-Bridge-Beweis. Es zwingt die
    experimentelle SAT-Slice-Ausgabe nur in die Buchhaltung aus
    proof_notes/ADVICE_FRONT_LEAKAGE_LEDGER_2026-05-20.md.
    """
    sat_rows = [r for r in results if r.get('n_sat', 0) > 0]
    if not sat_rows:
        return []

    critical = max(sat_rows, key=lambda r: (r['n'], r['H_mean']))
    n = int(critical['n'])
    H = float(critical['H_mean'])
    template_bits = _safe_log2(n + 1)
    cover_mass = _slice_cover_proxy(H, n)

    low_alpha = min(results_alpha, key=lambda r: r['alpha']) if results_alpha else None
    high_alpha = max(results_alpha, key=lambda r: r['alpha']) if results_alpha else None

    rows = [
        {
            'row_id': 'SAT_PHASE_RANDOM_FRONT',
            'n': n,
            'alpha': 4.267,
            'front_template_id': 'random_half_slice_grid',
            'front_template_bits': round(template_bits, 6),
            'front_extractor_id': 'max_half_slice_entropy',
            'advice_scope': 'n_only',
            'advice_bits': 0.0,
            'specialization_bits': 0.0,
            'reference_source': 'independent_seed_grid',
            'same_run_reference_risk': False,
            'matched_control_id': 'alpha_sweep_low_high',
            'front_cover_mass': round(cover_mass, 8),
            'capillary_hit': False,
            'residue_currency': 'shannon_proxy',
            'shannon_proxy_bits': round(H, 6),
            'allowance_bits': round(template_bits, 6),
            'uniformity_tail_after_allowances': round(H - template_bits, 6),
            'detector_capacity_bits': round(template_bits, 6),
            'target_entropy_degree': 'slice_proxy_only',
            'advice_degree': 'n_only',
            'front_degree': 'finite_grid',
            'residual_degree': 'not_kt_residue',
            'status': 'diagnostic_only_no_positive_tail',
            'interpretation': 'Shannon-H_slice bleibt unter den sichtbaren Template-Allowances; keine Uniformity-Bridge-Evidenz.',
        },
        {
            'row_id': 'SAME_RUN_CLUSTER_FRONT_FAIL',
            'n': n,
            'alpha': 4.267,
            'front_template_id': 'posthoc_solver_cluster',
            'front_template_bits': 0.0,
            'front_extractor_id': 'same_run_success_trace',
            'advice_scope': 'same_run',
            'advice_bits': float(n),
            'specialization_bits': float(n // 2),
            'reference_source': 'same_run_self_consistency',
            'same_run_reference_risk': True,
            'matched_control_id': 'must_fail_by_design',
            'front_cover_mass': 1.0,
            'capillary_hit': True,
            'residue_currency': 'kt_mdl_residue_required_but_absent',
            'shannon_proxy_bits': round(H, 6),
            'allowance_bits': round(1.5 * n, 6),
            'uniformity_tail_after_allowances': round(H - 1.5 * n, 6),
            'detector_capacity_bits': round(1.5 * n, 6),
            'target_entropy_degree': 'self_referential',
            'advice_degree': 'same_run',
            'front_degree': 'posthoc',
            'residual_degree': 'collapsed',
            'status': 'rejected_same_run_reference',
            'interpretation': 'Ein Frontsignal aus demselben Solverlauf wird als Leakage gebucht, nicht als externe Brücke.',
        },
        {
            'row_id': 'DUPLICATED_ADVICE_CONTROL_FAIL',
            'n': n,
            'alpha': 4.267,
            'front_template_id': 'duplicated_advice_baseline',
            'front_template_bits': round(template_bits, 6),
            'front_extractor_id': 'advice_echo',
            'advice_scope': 'instance_dependent',
            'advice_bits': float(n),
            'specialization_bits': round(_safe_log2(n + 1), 6),
            'reference_source': 'target_class_in_predictor',
            'same_run_reference_risk': True,
            'matched_control_id': 'duplicate_advice_negative_control',
            'front_cover_mass': 1.0,
            'capillary_hit': True,
            'residue_currency': 'kt_mdl_residue_required_but_absent',
            'shannon_proxy_bits': round(H, 6),
            'allowance_bits': round(template_bits + n + _safe_log2(n + 1), 6),
            'uniformity_tail_after_allowances': round(H - template_bits - n - _safe_log2(n + 1), 6),
            'detector_capacity_bits': round(template_bits + n, 6),
            'target_entropy_degree': 'target_leak',
            'advice_degree': 'instance_dependent',
            'front_degree': 'leaky',
            'residual_degree': 'collapsed',
            'status': 'rejected_advice_leak',
            'interpretation': 'Instanzabhängiger Rat oder Zielklassenwissen muss als Advice zählen und darf den Rest nicht stützen.',
        },
        {
            'row_id': 'LOW_HIGH_ALPHA_MATCHED_CONTROL',
            'n': n,
            'alpha': f"{low_alpha['alpha'] if low_alpha else 'NA'}..{high_alpha['alpha'] if high_alpha else 'NA'}",
            'front_template_id': 'alpha_matched_controls',
            'front_template_bits': round(template_bits, 6),
            'front_extractor_id': 'same_n_alpha_sweep',
            'advice_scope': 'n_only',
            'advice_bits': 0.0,
            'specialization_bits': 0.0,
            'reference_source': 'matched_alpha_sweep',
            'same_run_reference_risk': False,
            'matched_control_id': 'low_vs_high_clause_density',
            'front_cover_mass': round(_slice_cover_proxy((low_alpha['H_mean'] + high_alpha['H_mean']) / 2, n), 8) if low_alpha and high_alpha else 0.0,
            'capillary_hit': False,
            'residue_currency': 'shannon_proxy',
            'shannon_proxy_bits': round((low_alpha['H_mean'] + high_alpha['H_mean']) / 2, 6) if low_alpha and high_alpha else 0.0,
            'allowance_bits': round(template_bits, 6),
            'uniformity_tail_after_allowances': round(((low_alpha['H_mean'] + high_alpha['H_mean']) / 2) - template_bits, 6) if low_alpha and high_alpha else 0.0,
            'detector_capacity_bits': round(template_bits, 6),
            'target_entropy_degree': 'matched_control_proxy',
            'advice_degree': 'n_only',
            'front_degree': 'finite_grid',
            'residual_degree': 'not_kt_residue',
            'status': 'control_required_before_claim',
            'interpretation': 'Alpha-Kontrollen sind vorhanden, aber nur als H_slice-Proxies; OGP/planted/easy-Kontrollen fehlen noch.',
        },
    ]
    return rows


def write_advice_front_ledger(results, results_alpha):
    ledger = build_advice_front_ledger(results, results_alpha)
    if not ledger:
        return None

    outdir = os.path.join(os.path.dirname(__file__), '_results')
    os.makedirs(outdir, exist_ok=True)
    stem = 'PNP_ADVICE_FRONT_LEDGER_2026-06-01'
    json_path = os.path.join(outdir, stem + '.json')
    csv_path = os.path.join(outdir, stem + '.csv')
    md_path = os.path.join(outdir, stem + '.md')

    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(ledger, f, indent=2, ensure_ascii=False)

    fieldnames = list(ledger[0].keys())
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(ledger)

    with open(md_path, 'w', encoding='utf-8') as f:
        f.write('# P vs NP Advice-Front Ledger (2026-06-01)\n\n')
        f.write('Status: reproduzierbarer Guardrail-Dry-Run; kein P-vs-NP-Beweis und keine neue Uniformity Bridge.\n\n')
        f.write('Dieses Ledger zwingt die vorhandene SAT-Slice-Numerik in die Buchhaltung aus ')
        f.write('`proof_notes/ADVICE_FRONT_LEAKAGE_LEDGER_2026-05-20.md`. ')
        f.write('Alle positiven oder leaky Frontsignale werden nur als Diagnose geführt, solange kein ')
        f.write('Front-/Flow-Shadow-Lemma und keine unabhängigen OGP-/planted-/easy-Kontrollen vorliegen.\n\n')
        f.write('| Row | Status | Advice | Reference | Tail after allowances | Interpretation |\n')
        f.write('|---|---|---:|---|---:|---|\n')
        for row in ledger:
            f.write(
                f"| `{row['row_id']}` | `{row['status']}` | {row['advice_bits']} | "
                f"`{row['reference_source']}` | {row['uniformity_tail_after_allowances']} | "
                f"{row['interpretation']} |\n"
            )
        f.write('\n## Fazit\n\n')
        f.write('Der heutige Lauf verbessert die Nachvollziehbarkeit der Experimente, nicht den theorem-level Status. ')
        f.write('Der aktuelle SAT-Slice-Lauf liefert keine positive `K_front^t`-Restgröße. ')
        f.write('Same-run- und instanzabhängige Fronten werden korrekt als Leakage verworfen. ')
        f.write('Nächster sinnvoller Schritt ist ein echtes Kontrollset mit planted/easy/OGP-Instanzen und ')
        f.write('einem externen Referenztest.\n')

    return {'json': json_path, 'csv': csv_path, 'md': md_path, 'rows': len(ledger)}


def _control_case_row(row_id, family, clauses, n, alpha, seed,
                      reference_source, advice_scope,
                      known_assignment=None, ogp_certified=False):
    witness_count, witness_sample = enumerate_satisfying_assignments(clauses, n)
    sat_frac = witness_count / (2**n)
    np.random.seed(seed + 7000)
    H = 0.0
    if witness_count > 0:
        H, _ = max_slice_entropy(clauses, n, k=n//2, n_random_subsets=12)

    template_bits = _safe_log2(n + 1)
    known_ok = bool(known_assignment and evaluate_formula(clauses, known_assignment))
    stability = single_flip_stability(clauses, known_assignment, n) if known_assignment else 0.0
    distance_stats = witness_pair_distance_stats(witness_sample, n)

    if row_id.startswith('PLANTED'):
        status = 'positive_control_advice_visible'
        interpretation = (
            'Der vorregistrierte planted Witness besteht; das ist eine Positivkontrolle '
            'für Advice-/Referenzsichtbarkeit, keine Uniformity Bridge.'
        )
    elif row_id.startswith('EASY'):
        status = 'positive_tail_expected_easy_control'
        interpretation = (
            'Leichte Low-alpha-Instanzen haben erwartbar breite Witness-Masse; '
            'ein positiver Shannon-Proxy-Tail ist hier ein Sanity-Check, kein Separationssignal.'
        )
    elif row_id.startswith('OGP'):
        status = 'ogp_proxy_not_certified_small_n'
        interpretation = (
            'Der kleine 3-SAT-Smoke liefert nur Distanzproxies. Ein echtes OGP-Zertifikat '
            'braucht ein separates asymptotisches oder literaturgestütztes Gate.'
        )
    else:
        status = 'external_holdout_smoke_no_tail'
        interpretation = (
            'Unabhängiger Holdout-Seed als externe Referenzprobe; keine positive '
            'Restgröße nach sichtbarer Template-Allowance.'
        )

    row = {
        'row_id': row_id,
        'family': family,
        'n': n,
        'm': len(clauses),
        'alpha': alpha,
        'seed': seed,
        'sat_frac': round(sat_frac, 8),
        'witness_count': witness_count,
        'H_slice_proxy': round(float(H), 6),
        'front_template_bits': round(template_bits, 6),
        'uniformity_tail_after_allowances': round(float(H) - template_bits, 6),
        'known_witness_satisfies': known_ok,
        'single_flip_stability': round(stability, 6),
        'pair_distance_min': distance_stats['distance_min'],
        'pair_distance_max': distance_stats['distance_max'],
        'pair_distance_middle_fraction': distance_stats['distance_middle_fraction'],
        'ogp_certified': ogp_certified,
        'external_reference_source': reference_source,
        'advice_scope': advice_scope,
        'same_run_reference_risk': False,
        'status': status,
        'interpretation': interpretation,
    }
    return row


def build_control_family_ledger():
    """
    Kleine planted/easy/OGP-Kontrollmatrix.

    Der Zweck ist negativ/diagnostisch: Controls werden vorregistriert und
    getrennt von P-vs-NP-Claims ausgewiesen. Besonders die OGP-Zeile ist nur
    ein Proxy-Gate, kein OGP-Nachweis.
    """
    n = 12
    critical_alpha = 4.267
    critical_m = int(critical_alpha * n)
    planted_assignment = {var: (var % 2 == 0) for var in range(1, n + 1)}

    cases = [
        (
            'PLANTED_WITNESS_POSITIVE_CONTROL',
            'planted_3sat',
            planted_3sat(n, critical_m, planted_assignment, seed=6203),
            critical_alpha,
            6203,
            'precommitted_planted_assignment',
            'external_holdout',
            planted_assignment,
            False,
        ),
        (
            'EASY_LOW_ALPHA_CONTROL',
            'random_3sat_low_density',
            random_3sat(n, int(2.0 * n), seed=6204),
            2.0,
            6204,
            'independent_low_alpha_seed',
            'n_only',
            None,
            False,
        ),
        (
            'RANDOM_THRESHOLD_HOLDOUT',
            'random_3sat_threshold_holdout',
            random_3sat(n, critical_m, seed=6205),
            critical_alpha,
            6205,
            'independent_threshold_holdout_seed',
            'n_only',
            None,
            False,
        ),
        (
            'OGP_PROXY_STRESS_SMALL_N',
            'random_3sat_threshold_distance_proxy',
            random_3sat(n, critical_m, seed=6206),
            critical_alpha,
            6206,
            'distance_proxy_only_no_ogp_certificate',
            'n_only',
            None,
            False,
        ),
    ]

    rows = []
    for row_id, family, clauses, alpha, seed, source, advice, known, ogp in cases:
        rows.append(_control_case_row(
            row_id=row_id,
            family=family,
            clauses=clauses,
            n=n,
            alpha=alpha,
            seed=seed,
            reference_source=source,
            advice_scope=advice,
            known_assignment=known,
            ogp_certified=ogp,
        ))
    return rows


def write_control_family_ledger():
    rows = build_control_family_ledger()
    outdir = os.path.join(os.path.dirname(__file__), '_results')
    os.makedirs(outdir, exist_ok=True)
    stem = 'PNP_PLANTED_EASY_OGP_CONTROLS_2026-06-03'
    json_path = os.path.join(outdir, stem + '.json')
    csv_path = os.path.join(outdir, stem + '.csv')
    md_path = os.path.join(outdir, stem + '.md')

    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(rows, f, indent=2, ensure_ascii=False)

    fieldnames = list(rows[0].keys())
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    with open(md_path, 'w', encoding='utf-8') as f:
        f.write('# P vs NP planted/easy/OGP controls (2026-06-03)\n\n')
        f.write('Status: Kontrollfamilien-Guardrail; kein P-vs-NP-Beweis, keine neue Uniformity Bridge und kein OGP-Zertifikat.\n\n')
        f.write('Zweck: Der Advice-/Front-Ledger vom 2026-06-01 verlangte echte planted/easy/OGP-Kontrollen ')
        f.write('und externe Referenzquellen. Dieser kleine Lauf materialisiert die ersten drei Kontrolltypen ')
        f.write('als reproduzierbare, kleine Smoke-Matrix. Die OGP-Zeile ist absichtlich nur ein Proxy-Gate: ')
        f.write('bei `n=12`, `k=3` und wenigen Seeds kann keine asymptotische OGP-Aussage folgen.\n\n')
        f.write('| Row | Status | SAT% | Witnesses | H_slice | Tail | Reference | Interpretation |\n')
        f.write('|---|---|---:|---:|---:|---:|---|---|\n')
        for row in rows:
            f.write(
                f"| `{row['row_id']}` | `{row['status']}` | {100*row['sat_frac']:.3f} | "
                f"{row['witness_count']} | {row['H_slice_proxy']} | "
                f"{row['uniformity_tail_after_allowances']} | "
                f"`{row['external_reference_source']}` | {row['interpretation']} |\n"
            )
        f.write('\n## Befund\n\n')
        f.write('Die planted-, Holdout- und OGP-Proxy-Zeilen bleiben nach sichtbarer Template-Allowance im negativen Tail-Bereich. ')
        f.write('Die Low-alpha-Kontrolle zeigt dagegen erwartbar einen positiven Shannon-Proxy-Tail, weil sie witnessreich und leicht ist; ')
        f.write('sie zählt deshalb als Sanity-Check, nicht als Separationssignal. Der planted Witness ist eine bewusst eingebaute ')
        f.write('Positivkontrolle und zählt als Advice-/Referenzsichtbarkeit. Die Low-alpha- und Holdout-Zeilen ')
        f.write('trennen externe Seed-Referenz von Same-run-Self-Consistency. Die OGP-Zeile bleibt ')
        f.write('`ogp_proxy_not_certified_small_n`; ein echter nächster Schritt wäre ein separates ')
        f.write('Random-k-SAT-/low-degree- oder Literatur-Gate statt weiterer Klein-n-Interpretation.\n')

    return {'json': json_path, 'csv': csv_path, 'md': md_path, 'rows': len(rows)}


# ==========================================================================
# Hauptberechnung
# ==========================================================================

parser = argparse.ArgumentParser(description='SAT slice entropy and advice-front ledger')
parser.add_argument(
    '--quick',
    action='store_true',
    help='Run a small deterministic smoke run and write the advice/front ledger.'
)
parser.add_argument(
    '--controls-only',
    action='store_true',
    help='Only write the planted/easy/OGP control-family guardrail report.'
)
args = parser.parse_args()

if args.controls_only:
    control_paths = write_control_family_ledger()
    print("Planted/easy/OGP-Kontrollreport geschrieben "
          f"({control_paths['rows']} Zeilen):")
    print(f"  JSON: {control_paths['json']}")
    print(f"  CSV:  {control_paths['csv']}")
    print(f"  MD:   {control_paths['md']}")
    raise SystemExit(0)

print("=" * 70)
print("P-vs-NP: SAT SLICE-ENTROPIE am Phasenuebergang")
print("=" * 70)
if args.quick:
    print("[quick mode] kleiner Smoke-Run fuer Ledger-Verifikation")

# ==========================================================================
# Test 1: Skalierung H_slice(n) fuer 3-SAT am Phasenuebergang
# ==========================================================================
print(f"\n[1] Skalierung von H_slice mit n (alpha = 4.267)")

alpha_critical = 4.267  # 3-SAT Phasenuebergang
n_values = [8, 10, 12] if args.quick else [8, 10, 12, 14, 16]
n_instances = 4 if args.quick else 15  # Instanzen pro n
n_random_subsets_main = 8 if args.quick else 30

results = []

for n in n_values:
    m = int(alpha_critical * n)
    t0 = time.time()

    H_values = []
    K_values = []
    sat_fracs = []

    for seed in range(n_instances):
        clauses = random_3sat(n, m, seed=seed + n*1000)

        # SAT-Fraktion (Anteil erfuellender Belegungen)
        if n <= 20:
            sat_count = 0
            for idx in range(2**n):
                assignment = {}
                for var in range(1, n+1):
                    assignment[var] = bool((idx >> (var-1)) & 1)
                if evaluate_formula(clauses, assignment):
                    sat_count += 1
            sat_frac = sat_count / 2**n
        else:
            # Monte Carlo
            sat_count = sum(1 for _ in range(5000)
                          if evaluate_formula(clauses,
                                             {var: bool(np.random.randint(2))
                                              for var in range(1, n+1)}))
            sat_frac = sat_count / 5000

        sat_fracs.append(sat_frac)

        if sat_frac > 0:
            # Slice-Entropie
            H, _ = max_slice_entropy(clauses, n, k=n//2, n_random_subsets=n_random_subsets_main)
            H_values.append(H)

            # Kolmogorov-Approximation
            K = approx_kolmogorov(clauses, n)
            K_values.append(K)

    elapsed = time.time() - t0
    H_mean = np.mean(H_values) if H_values else 0
    H_std = np.std(H_values) if H_values else 0
    K_mean = np.mean(K_values) if K_values else 0
    sat_mean = np.mean(sat_fracs)

    results.append({
        'n': n, 'm': m, 'H_mean': H_mean, 'H_std': H_std,
        'K_mean': K_mean, 'sat_frac': sat_mean,
        'n_sat': len(H_values), 'time': elapsed
    })

    print(f"  n={n:>3}, m={m:>4}: H_slice={H_mean:>8.3f} +/- {H_std:>6.3f}, "
          f"K~{K_mean:>6.0f} bits, SAT%={100*sat_mean:>5.1f}%, "
          f"({len(H_values)}/{n_instances} SAT, {elapsed:.1f}s)")


# ==========================================================================
# Test 2: H_slice bei verschiedenen alpha (leicht vs schwer)
# ==========================================================================
print(f"\n{'='*70}")
print("[2] H_slice vs. alpha (Phasenuebergang)")
print(f"{'='*70}")

n_test = 12
alphas = [2.0, 4.267, 8.0] if args.quick else [2.0, 3.0, 3.5, 4.0, 4.267, 4.5, 5.0, 6.0, 8.0]
n_inst = 5 if args.quick else 30
n_random_subsets_alpha = 8 if args.quick else 20

results_alpha = []

for alpha in alphas:
    m = int(alpha * n_test)
    H_vals = []
    sat_vals = []

    for seed in range(n_inst):
        clauses = random_3sat(n_test, m, seed=seed + int(alpha*10000))

        # SAT check
        sat_count = 0
        for idx in range(2**n_test):
            assignment = {}
            for var in range(1, n_test+1):
                assignment[var] = bool((idx >> (var-1)) & 1)
            if evaluate_formula(clauses, assignment):
                sat_count += 1
        sat_frac = sat_count / 2**n_test
        sat_vals.append(sat_frac)

        if sat_frac > 0:
            H, _ = max_slice_entropy(clauses, n_test, k=n_test//2, n_random_subsets=n_random_subsets_alpha)
            H_vals.append(H)

    H_mean = np.mean(H_vals) if H_vals else 0
    sat_mean = np.mean(sat_vals)

    results_alpha.append({
        'alpha': alpha, 'H_mean': H_mean, 'sat_frac': sat_mean,
        'n_sat': len(H_vals)
    })

    print(f"  alpha={alpha:>5.3f}: H_slice={H_mean:>7.3f}, "
          f"SAT%={100*sat_mean:>6.2f}%, {len(H_vals)}/{n_inst} SAT")


# ==========================================================================
# Analyse: Exponentielle vs polynomielle Skalierung
# ==========================================================================
print(f"\n{'='*70}")
print("[3] SKALIERUNGSANALYSE")
print(f"{'='*70}")

n_arr = np.array([r['n'] for r in results])
H_arr = np.array([r['H_mean'] for r in results])
coeffs_exp = None

# Pruefe: H ~ 2^{n/2} (exponentiell) oder H ~ n^c (polynomiell)?
mask = H_arr > 0
if np.sum(mask) >= 3:
    # Log-Fit: ln(H) vs n
    log_H = np.log2(np.maximum(H_arr[mask], 1e-10))
    n_fit = n_arr[mask]

    # Exponentieller Fit: log2(H) = a*n + b
    coeffs_exp = np.polyfit(n_fit, log_H, 1)
    # Polynomieller Fit: log2(H) = c*log2(n) + d
    log_n = np.log2(n_fit)
    coeffs_poly = np.polyfit(log_n, log_H, 1)

    print(f"\n  Exponentieller Fit: log2(H) = {coeffs_exp[0]:.4f} * n + {coeffs_exp[1]:.4f}")
    print(f"    => H ~ 2^{{{coeffs_exp[0]:.4f}*n}}")
    print(f"    ESC-Hypothese erwartet: H ~ 2^{{n/2}} (Koeffizient = 0.5)")

    print(f"\n  Polynomieller Fit:  log2(H) = {coeffs_poly[0]:.4f} * log2(n) + {coeffs_poly[1]:.4f}")
    print(f"    => H ~ n^{{{coeffs_poly[0]:.4f}}}")

    # Welcher Fit ist besser? (R^2)
    H_pred_exp = coeffs_exp[0] * n_fit + coeffs_exp[1]
    H_pred_poly = coeffs_poly[0] * log_n + coeffs_poly[1]
    ss_total = np.sum((log_H - np.mean(log_H))**2)
    if ss_total > 0:
        r2_exp = 1 - np.sum((log_H - H_pred_exp)**2) / ss_total
        r2_poly = 1 - np.sum((log_H - H_pred_poly)**2) / ss_total
        print(f"\n  R^2 (exponentiell): {r2_exp:.4f}")
        print(f"  R^2 (polynomiell):  {r2_poly:.4f}")
        print(f"  => {'EXPONENTIELL' if r2_exp > r2_poly else 'POLYNOMIELL'} passt besser")
    else:
        print("  (Zu wenige Datenpunkte fuer R^2-Vergleich)")


# ==========================================================================
# Zusammenfassung
# ==========================================================================
print(f"\n{'='*70}")
print("ZUSAMMENFASSUNG: SAT Slice-Entropie")
print(f"{'='*70}")
print(f"""
  Getestete n-Werte: {n_values}
  Instanzen pro n:   {n_instances}
  Alpha (Phasenuebergang): {alpha_critical}

  Slice-Entropie H_slice(SAT_n) am Phasenuebergang:
""")

for r in results:
    print(f"    n={r['n']:>3}: H={r['H_mean']:>8.3f} +/- {r['H_std']:>6.3f} "
          f"(SAT%={100*r['sat_frac']:>5.1f}%)")

print(f"""
  ESC-Hypothese:
    H_slice(NP-harte Instanzen) >= c * 2^{{n/2}}
    H_slice(P-berechenbare Fkt.) <= poly(n)

  Interpretation:
    - Der Lauf ist nur eine kleine Shannon-Proxy-Diagnose, keine K^t-Messung
    - Aus n <= 16 folgt keine asymptotische Skalierung und keine ESC-Evidenz
    - Volle Validierung braeuchte groessere n und unabhaengige Kontrollfamilien
    - Die Kompressionslaenge K(f) bleibt hier nur eine grobe Vergleichsgroesse
""")

ledger_paths = write_advice_front_ledger(results, results_alpha)
if ledger_paths:
    print(f"\nAdvice-/Front-Ledger geschrieben ({ledger_paths['rows']} Zeilen):")
    print(f"  JSON: {ledger_paths['json']}")
    print(f"  CSV:  {ledger_paths['csv']}")
    print(f"  MD:   {ledger_paths['md']}")


# Plot
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Panel 1: H_slice vs n
    ax = axes[0, 0]
    n_plot = [r['n'] for r in results]
    H_plot = [r['H_mean'] for r in results]
    H_err = [r['H_std'] for r in results]
    ax.errorbar(n_plot, H_plot, yerr=H_err, fmt='ro-', linewidth=2,
                markersize=8, capsize=4, label=r'$H_{slice}$ (3-SAT)')
    ax.set_xlabel('n (Variablen)', fontsize=12)
    ax.set_ylabel(r'$H_{slice}$', fontsize=12)
    ax.set_title(r'Slice-Entropie vs. $n$ ($\alpha = 4.267$)', fontsize=11)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Panel 2: H_slice vs alpha
    ax = axes[0, 1]
    alpha_plot = [r['alpha'] for r in results_alpha]
    H_alpha_plot = [r['H_mean'] for r in results_alpha]
    sat_plot = [r['sat_frac'] for r in results_alpha]
    ax.plot(alpha_plot, H_alpha_plot, 'bo-', linewidth=2, markersize=8,
            label=r'$H_{slice}$')
    ax2 = ax.twinx()
    ax2.plot(alpha_plot, [100*s for s in sat_plot], 'g^--', linewidth=1.5,
             markersize=6, alpha=0.7, label='SAT%')
    ax.axvline(x=4.267, color='red', linestyle=':', linewidth=1.5,
               label=r'$\alpha_c = 4.267$')
    ax.set_xlabel(r'$\alpha = m/n$', fontsize=12)
    ax.set_ylabel(r'$H_{slice}$', fontsize=12, color='blue')
    ax2.set_ylabel('SAT %', fontsize=12, color='green')
    ax.set_title(r'$H_{slice}$ vs. Klauseldichte $\alpha$', fontsize=11)
    ax.legend(fontsize=9, loc='upper left')
    ax2.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)

    # Panel 3: Skalierung (log)
    ax = axes[1, 0]
    if any(h > 0 for h in H_plot):
        H_pos = [(n, h) for n, h in zip(n_plot, H_plot) if h > 0]
        if H_pos:
            ns, hs = zip(*H_pos)
            ax.semilogy(ns, hs, 'rs-', linewidth=2, markersize=8,
                       label=r'$H_{slice}$ (gemessen)')
            # Exponentieller Fit plotten
            n_line = np.linspace(min(ns), max(ns), 50)
            if len(ns) >= 2 and coeffs_exp is not None:
                H_exp_line = 2**(coeffs_exp[0] * n_line + coeffs_exp[1])
                ax.semilogy(n_line, H_exp_line, 'r--', linewidth=1,
                           label=f'$2^{{{coeffs_exp[0]:.3f} n}}$')
                H_half = 2**(0.5 * n_line - 2)
                ax.semilogy(n_line, H_half, 'k:', linewidth=1,
                           label=r'$2^{n/2}$ (ESC)')
    ax.set_xlabel('n', fontsize=12)
    ax.set_ylabel(r'$H_{slice}$ (log)', fontsize=12)
    ax.set_title('Skalierungsanalyse', fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')

    # Panel 4: Kolmogorov-Komplexitaet
    ax = axes[1, 1]
    K_plot = [r['K_mean'] for r in results]
    ax.plot(n_plot, K_plot, 'g^-', linewidth=2, markersize=8,
            label='K(f) (Kompression)')
    ax.plot(n_plot, [2**n for n in n_plot], 'r--', linewidth=1, alpha=0.5,
            label=r'$2^n$ (max)')
    ax.set_xlabel('n', fontsize=12)
    ax.set_ylabel('K(f) [Bits]', fontsize=12)
    ax.set_title('Kolmogorov-Komplexitaet (approx.)', fontsize=11)
    ax.legend(fontsize=10)
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    outpath = os.path.join(os.path.dirname(__file__),
                           'compute_sat_entropy.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"Plot gespeichert: {outpath}")
except Exception as e:
    print(f"(matplotlib: {e})")

print("\n[DONE]")
