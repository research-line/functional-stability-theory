"""Finite coercive-complement ledger for the Yang--Mills RG programme.

The ledger evaluates the three quantities introduced by the conditional
three-lemma schema in the bilingual Yang--Mills paper::

    s_lambda = |<r, (H-E)q>|,
    p_lambda = ||(I-|r><r|)(H-E)q||,
    g_*      = dist(E, spectrum(H on Ran(I-P))).

It then checks the finite-dimensional leakage bound

    ||(I-P)q|| <= (s_lambda + p_lambda) / g_*.

Only self-adjoint matrices and orthogonal projections that reduce the matrix
are accepted.  This is deliberate: residual bounds in spectral gaps depend on
self-adjointness and associated subspaces (Seelmann, arXiv:2309.07032), while
projected finite-dimensional spectra can otherwise suffer spectral pollution
(Levitin--Shargorodsky, arXiv:math/0212087).

The bundled cases are synthetic matched controls.  They calibrate the gate but
do not identify a physical Yang--Mills microcluster, prove a scale-uniform
outer gap, or upgrade any mass-gap claim.

Numerical gate status and transfer eligibility are intentionally separate.  A
post-hoc cluster can have good residual/gap numbers and still be rejected as
circular; review eligibility additionally requires an independently verified
physical predefinition certificate for both target and complement.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np


ROOT = Path(__file__).resolve().parent
DATE = "2026-08-26"
STEM = f"YM_COERCIVE_COMPLEMENT_LEDGER_{DATE}"


@dataclass(frozen=True)
class AuditCase:
    """One finite, explicitly labelled coercive-complement control."""

    case_id: str
    role: str
    complement_kind: str
    hamiltonian: np.ndarray
    candidate: np.ndarray
    reference_energy: float
    relevant_direction: np.ndarray
    cluster_projection: np.ndarray
    cluster_predefined: bool
    projection_predefined: bool
    cluster_fixed_before_leakage: bool
    complement_fixed_before_leakage: bool
    cluster_definition_source: str
    predefinition_scope: str
    predefinition_certificate_id: str
    predefinition_certificate_status: str
    mean_log_tau_B: float
    notes: str


@dataclass(frozen=True)
class AuditResult:
    """Serializable measurements and fail-closed diagnostic status."""

    case_id: str
    role: str
    complement_kind: str
    cluster_predefined: bool
    projection_predefined: bool
    cluster_fixed_before_leakage: bool
    complement_fixed_before_leakage: bool
    cluster_definition_source: str
    predefinition_scope: str
    predefinition_certificate_id: str
    predefinition_certificate_status: str
    mean_log_tau_B: float
    reference_energy: float
    self_adjoint_error: float
    projection_error: float
    reducing_error: float
    s_lambda: float
    p_lambda: float
    residual_norm: float
    residual_sum: float
    g_star: float
    residual_over_gap: float
    actual_complement_leakage: float
    residual_split_bound_holds: bool
    leakage_bound_holds: bool
    outer_gap_floor: float
    residual_ratio_ceiling: float
    diagnostic_status: str
    transfer_decision: str
    claim_status: str
    notes: str


FIELDS = tuple(AuditResult.__dataclass_fields__)


def _finite_array(value: np.ndarray, *, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=np.complex128)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} enthält nicht-endliche Werte")
    return array


def _validate_thresholds(outer_gap_floor: float, residual_ratio_ceiling: float, tolerance: float) -> None:
    if not np.isfinite(outer_gap_floor) or outer_gap_floor <= 0.0:
        raise ValueError("outer_gap_floor muss endlich und strikt positiv sein")
    if not np.isfinite(residual_ratio_ceiling) or residual_ratio_ceiling < 0.0:
        raise ValueError("residual_ratio_ceiling muss endlich und nichtnegativ sein")
    if not np.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("tolerance muss endlich und strikt positiv sein")


def decide_transfer(case: AuditCase, diagnostic_status: str) -> str:
    """Separate numerical gate quality from non-circular transfer eligibility."""

    fixed_before_measurement = (
        case.cluster_predefined
        and case.projection_predefined
        and case.cluster_fixed_before_leakage
        and case.complement_fixed_before_leakage
    )
    if not fixed_before_measurement:
        return "reject_circular_target_or_complement"
    if diagnostic_status != "finite_control_pass":
        return "reject_numeric_coercive_complement_gate"
    if (
        case.predefinition_certificate_status != "independently_verified"
        or not case.predefinition_certificate_id
    ):
        return "blocked_predefinition_not_independently_verified"
    if case.predefinition_scope != "physical_yang_mills":
        return "blocked_nonphysical_predefinition_scope"
    return "eligible_for_analytic_review_no_claim"


def audit_case(
    case: AuditCase,
    *,
    outer_gap_floor: float = 0.1,
    residual_ratio_ceiling: float = 0.25,
    tolerance: float = 1e-10,
) -> AuditResult:
    """Compute ``s_lambda``, ``p_lambda``, ``g_*`` and the leakage ratio.

    The function fails closed if the operator is not self-adjoint, the supplied
    cluster is not an orthogonal projection, or the projection does not reduce
    the operator.  No compressed-spectrum result is reported in those cases.
    """

    _validate_thresholds(outer_gap_floor, residual_ratio_ceiling, tolerance)

    hamiltonian = _finite_array(case.hamiltonian, name="hamiltonian")
    candidate = _finite_array(case.candidate, name="candidate")
    relevant = _finite_array(case.relevant_direction, name="relevant_direction")
    projection = _finite_array(case.cluster_projection, name="cluster_projection")

    if hamiltonian.ndim != 2 or hamiltonian.shape[0] != hamiltonian.shape[1]:
        raise ValueError("hamiltonian muss quadratisch sein")
    dimension = hamiltonian.shape[0]
    if candidate.shape != (dimension,) or relevant.shape != (dimension,):
        raise ValueError("candidate und relevant_direction müssen zur Operatordimension passen")
    if projection.shape != (dimension, dimension):
        raise ValueError("cluster_projection muss zur Operatordimension passen")
    if not np.isfinite(case.reference_energy):
        raise ValueError("reference_energy muss endlich sein")
    if not np.isfinite(case.mean_log_tau_B):
        raise ValueError("mean_log_tau_B muss endlich sein")
    allowed_scopes = {"synthetic_fixture", "physical_yang_mills"}
    if case.predefinition_scope not in allowed_scopes:
        raise ValueError(f"unbekannter predefinition_scope: {case.predefinition_scope}")
    allowed_certificate_statuses = {"missing", "synthetic_fixture_only", "independently_verified"}
    if case.predefinition_certificate_status not in allowed_certificate_statuses:
        raise ValueError(
            f"unbekannter predefinition_certificate_status: {case.predefinition_certificate_status}"
        )
    if not case.cluster_definition_source:
        raise ValueError("cluster_definition_source darf nicht leer sein")

    self_adjoint_error = float(np.linalg.norm(hamiltonian - hamiltonian.conj().T, ord=2))
    if self_adjoint_error > tolerance:
        raise ValueError("hamiltonian ist nicht selbstadjungiert")

    projection_adjoint_error = np.linalg.norm(projection - projection.conj().T, ord=2)
    projection_idempotence_error = np.linalg.norm(projection @ projection - projection, ord=2)
    projection_error = float(max(projection_adjoint_error, projection_idempotence_error))
    if projection_error > tolerance:
        raise ValueError("cluster_projection ist kein Orthogonalprojektor")

    candidate_norm = float(np.linalg.norm(candidate))
    relevant_norm = float(np.linalg.norm(relevant))
    if abs(candidate_norm - 1.0) > tolerance:
        raise ValueError("candidate muss normiert sein")
    if abs(relevant_norm - 1.0) > tolerance:
        raise ValueError("relevant_direction muss normiert sein")
    if np.linalg.norm(projection @ relevant - relevant) > tolerance:
        raise ValueError("relevant_direction muss im Zielcluster liegen")

    reducing_error = float(np.linalg.norm(hamiltonian @ projection - projection @ hamiltonian, ord=2))
    if reducing_error > tolerance:
        raise ValueError("cluster_projection reduziert hamiltonian nicht")

    projection_eigenvalues, projection_eigenvectors = np.linalg.eigh(projection)
    complement_basis = projection_eigenvectors[:, projection_eigenvalues < 0.5]
    if complement_basis.shape[1] == 0:
        raise ValueError("das äußere Komplement darf nicht leer sein")

    identity = np.eye(dimension, dtype=np.complex128)
    complement_projection = identity - projection
    residual = (hamiltonian - case.reference_energy * identity) @ candidate
    relevant_projection = np.outer(relevant, relevant.conj())

    s_lambda = float(abs(np.vdot(relevant, residual)))
    p_lambda = float(np.linalg.norm((identity - relevant_projection) @ residual))
    residual_norm = float(np.linalg.norm(residual))
    residual_sum = s_lambda + p_lambda
    residual_split_bound_holds = residual_norm <= residual_sum + 100.0 * tolerance

    complement_operator = complement_basis.conj().T @ hamiltonian @ complement_basis
    complement_eigenvalues = np.linalg.eigvalsh(complement_operator)
    g_star = float(np.min(np.abs(complement_eigenvalues - case.reference_energy)))
    actual_leakage = float(np.linalg.norm(complement_projection @ candidate))

    if g_star <= tolerance:
        residual_over_gap = float("inf")
        leakage_bound_holds = False
    else:
        residual_over_gap = residual_sum / g_star
        leakage_bound_holds = actual_leakage <= residual_over_gap + 100.0 * tolerance

    if not residual_split_bound_holds:
        diagnostic_status = "reject_bound_invariant"
    elif g_star < outer_gap_floor:
        diagnostic_status = "reject_outer_gap_collapse"
    elif not leakage_bound_holds:
        diagnostic_status = "reject_bound_invariant"
    elif residual_over_gap > residual_ratio_ceiling:
        diagnostic_status = "reject_residual_over_gap"
    else:
        diagnostic_status = "finite_control_pass"
    transfer_decision = decide_transfer(case, diagnostic_status)

    return AuditResult(
        case_id=case.case_id,
        role=case.role,
        complement_kind=case.complement_kind,
        cluster_predefined=case.cluster_predefined,
        projection_predefined=case.projection_predefined,
        cluster_fixed_before_leakage=case.cluster_fixed_before_leakage,
        complement_fixed_before_leakage=case.complement_fixed_before_leakage,
        cluster_definition_source=case.cluster_definition_source,
        predefinition_scope=case.predefinition_scope,
        predefinition_certificate_id=case.predefinition_certificate_id,
        predefinition_certificate_status=case.predefinition_certificate_status,
        mean_log_tau_B=float(case.mean_log_tau_B),
        reference_energy=float(case.reference_energy),
        self_adjoint_error=self_adjoint_error,
        projection_error=projection_error,
        reducing_error=reducing_error,
        s_lambda=s_lambda,
        p_lambda=p_lambda,
        residual_norm=residual_norm,
        residual_sum=residual_sum,
        g_star=g_star,
        residual_over_gap=float(residual_over_gap),
        actual_complement_leakage=actual_leakage,
        residual_split_bound_holds=bool(residual_split_bound_holds),
        leakage_bound_holds=bool(leakage_bound_holds),
        outer_gap_floor=float(outer_gap_floor),
        residual_ratio_ceiling=float(residual_ratio_ceiling),
        diagnostic_status=diagnostic_status,
        transfer_decision=transfer_decision,
        claim_status="diagnostic_only_no_yang_mills_claim",
        notes=case.notes,
    )


def _three_level_case(
    *,
    case_id: str,
    role: str,
    complement_kind: str,
    complement_energy: float,
    leakage: float,
    mean_log_tau_B: float,
    notes: str,
    cluster_predefined: bool = True,
    projection_predefined: bool = True,
    cluster_fixed_before_leakage: bool = True,
    complement_fixed_before_leakage: bool = True,
    cluster_definition_source: str = "synthetic_coordinate_basis_e0_e1",
    predefinition_scope: str = "synthetic_fixture",
    predefinition_certificate_id: str = "SYNTHETIC_FIXTURE_FIXED_BEFORE_MEASUREMENT",
    predefinition_certificate_status: str = "synthetic_fixture_only",
) -> AuditCase:
    if not 0.0 <= leakage < 1.0:
        raise ValueError("leakage muss im Intervall [0, 1) liegen")
    hamiltonian = np.diag([0.0, 1.0, complement_energy])
    candidate = np.array([0.0, np.sqrt(1.0 - leakage**2), leakage])
    relevant = np.array([0.0, 1.0, 0.0])
    projection = np.diag([1.0, 1.0, 0.0])
    return AuditCase(
        case_id=case_id,
        role=role,
        complement_kind=complement_kind,
        hamiltonian=hamiltonian,
        candidate=candidate,
        reference_energy=1.0,
        relevant_direction=relevant,
        cluster_projection=projection,
        cluster_predefined=cluster_predefined,
        projection_predefined=projection_predefined,
        cluster_fixed_before_leakage=cluster_fixed_before_leakage,
        complement_fixed_before_leakage=complement_fixed_before_leakage,
        cluster_definition_source=cluster_definition_source,
        predefinition_scope=predefinition_scope,
        predefinition_certificate_id=predefinition_certificate_id,
        predefinition_certificate_status=predefinition_certificate_status,
        mean_log_tau_B=mean_log_tau_B,
        notes=notes,
    )


def build_matched_controls() -> list[AuditCase]:
    """Return one positive and two mechanism-matched negative controls."""

    return [
        _three_level_case(
            case_id="finite_predefined_positive_control",
            role="positive_control",
            complement_kind="regular_outer_complement",
            complement_energy=4.0,
            leakage=0.05,
            mean_log_tau_B=-0.40,
            notes="Synthetic reducing cluster with g_*=3 and a small, exactly bounded leakage.",
        ),
        _three_level_case(
            case_id="negative_bad_scale_residual",
            role="matched_negative_control",
            complement_kind="bad_scale_complement",
            complement_energy=1.25,
            leakage=0.60,
            mean_log_tau_B=-0.35,
            notes="Negative mean_log_tau_B does not rescue a residual/gap ratio above the registered ceiling.",
        ),
        _three_level_case(
            case_id="negative_gribov_gap_collapse",
            role="matched_negative_control",
            complement_kind="gribov_near_zero_complement",
            complement_energy=1.00000001,
            leakage=0.10,
            mean_log_tau_B=-0.30,
            notes="A near-zero complement separation is rejected even when the finite leakage ratio looks small.",
        ),
        _three_level_case(
            case_id="negative_posthoc_cluster_good_numbers",
            role="matched_negative_control",
            complement_kind="posthoc_selected_complement",
            complement_energy=4.0,
            leakage=0.05,
            mean_log_tau_B=-0.40,
            notes="Numerically identical to the positive control, but selected after leakage inspection.",
            cluster_predefined=False,
            projection_predefined=False,
            cluster_fixed_before_leakage=False,
            complement_fixed_before_leakage=False,
            cluster_definition_source="selected_after_spectrum_and_leakage_inspection",
            predefinition_certificate_id="",
            predefinition_certificate_status="missing",
        ),
    ]


def build_results(
    *,
    outer_gap_floor: float = 0.1,
    residual_ratio_ceiling: float = 0.25,
    tolerance: float = 1e-10,
) -> list[AuditResult]:
    return [
        audit_case(
            case,
            outer_gap_floor=outer_gap_floor,
            residual_ratio_ceiling=residual_ratio_ceiling,
            tolerance=tolerance,
        )
        for case in build_matched_controls()
    ]


def _json_safe_row(result: AuditResult) -> dict[str, object]:
    row = asdict(result)
    if not np.isfinite(row["residual_over_gap"]):
        row["residual_over_gap"] = None
    return row


def write_outputs(results: Sequence[AuditResult], output_dir: Path) -> dict[str, Path]:
    """Write deterministic CSV, JSON and Markdown diagnostics."""

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "csv": output_dir / f"{STEM}.csv",
        "json": output_dir / f"{STEM}.json",
        "md": output_dir / f"{STEM}.md",
    }

    with paths["csv"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(asdict(result) for result in results)

    counts: dict[str, int] = {}
    transfer_counts: dict[str, int] = {}
    for result in results:
        counts[result.diagnostic_status] = counts.get(result.diagnostic_status, 0) + 1
        transfer_counts[result.transfer_decision] = transfer_counts.get(result.transfer_decision, 0) + 1
    payload = {
        "schema_version": 2,
        "date": DATE,
        "scope_status": "finite_synthetic_controls_only",
        "claim_pass": 0,
        "status_counts": counts,
        "transfer_decision_counts": transfer_counts,
        "rows": [_json_safe_row(result) for result in results],
    }
    paths["json"].write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    lines = [
        f"# Yang--Mills Coercive-Complement Ledger {DATE}",
        "",
        "Status: finite synthetische Gate-Kalibrierung; kein Yang--Mills- oder Kontinuumsclaim.",
        "",
        "| case | complement | mean log tau | s_lambda | p_lambda | g_* | (s+p)/g_* | leakage | numeric status | transfer decision |",
        "|---|---|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for result in results:
        ratio = "inf" if not np.isfinite(result.residual_over_gap) else f"{result.residual_over_gap:.8g}"
        lines.append(
            f"| {result.case_id} | {result.complement_kind} | {result.mean_log_tau_B:.6g} | "
            f"{result.s_lambda:.8g} | {result.p_lambda:.8g} | {result.g_star:.8g} | "
            f"{ratio} | {result.actual_complement_leakage:.8g} | {result.diagnostic_status} | "
            f"{result.transfer_decision} |"
        )
    lines.extend(
        [
            "",
            "## Guardrails",
            "",
            "- `mean_log_tau_B` is reported but never used as `g_*` or as a pass criterion.",
            "- Non-self-adjoint operators, non-orthogonal projectors, and non-reducing clusters fail closed.",
            "- The bad-scale control fails on `(s_lambda+p_lambda)/g_*`; the Gribov control fails on `g_*`.",
            "- Numerical status and transfer decision are separate: the post-hoc control has the same",
            "  good numbers as the positive fixture but is rejected as circular.",
            "- A usable transfer schema requires target cluster and complement to be fixed before leakage",
            "  measurement and an independently verified physical Yang--Mills predefinition certificate.",
            "- No bundled synthetic row is transfer-eligible; declaration fields are not certificates.",
            "- `claim_pass` remains 0.",
            "",
        ]
    )
    paths["md"].write_text("\n".join(lines), encoding="utf-8")
    return paths


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Finite Yang--Mills coercive-complement gate calibration")
    parser.add_argument("--outer-gap-floor", type=float, default=0.1)
    parser.add_argument("--residual-ratio-ceiling", type=float, default=0.25)
    parser.add_argument("--tolerance", type=float, default=1e-10)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "_results")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    results = build_results(
        outer_gap_floor=args.outer_gap_floor,
        residual_ratio_ceiling=args.residual_ratio_ceiling,
        tolerance=args.tolerance,
    )
    paths = write_outputs(results, args.output_dir)
    for result in results:
        print(
            f"{result.case_id}: g_*={result.g_star:.8g} "
            f"(s+p)/g_*={result.residual_over_gap:.8g} numeric={result.diagnostic_status} "
            f"transfer={result.transfer_decision}"
        )
    print(f"claim_pass=0 json={paths['json']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
