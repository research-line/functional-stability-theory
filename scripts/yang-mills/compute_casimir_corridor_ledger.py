"""Finite Casimir-corridor audit for the Yang--Mills mass-gap programme.

The Selberg donor pattern is a spectral window of a central Casimir inside a
predeclared symmetry sector.  The Yang--Mills analogue tested here is

    P_I = 1_I(C) P_phys,

where ``P_phys`` is the Gauss-law/gauge-invariant projector and ``I`` is fixed
before inspecting the Hamiltonian spectrum or its corridor commutator.

Non-linear Fourier analysis and spin-network/prepotential constructions give
orthonormal bases of the physical lattice-gauge Hilbert space (Burgio et al.,
arXiv:hep-lat/9906036; Mathur, arXiv:hep-lat/0405008).  This makes a finite
symmetry-sector diagnostic meaningful, but it does not imply that the full
Hamiltonian reduces an arbitrarily chosen Casimir window.  The ledger therefore
checks three logically separate gates:

1. the window is fixed non-circularly;
2. the Hamiltonian and Casimir preserve the physical gauge sector;
3. the Hamiltonian reduces the derived Casimir corridor.

The bundled rows comprise synthetic matched controls and one finite compact
U(1) strong-coupling loop-flux control with electric Casimir ``n^2``.  They
neither identify a physical four-dimensional non-Abelian Yang--Mills corridor
nor establish a continuum mass gap.
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
STEM = f"YM_CASIMIR_CORRIDOR_LEDGER_{DATE}"


@dataclass(frozen=True)
class AuditCase:
    """One finite gauge/Casimir corridor control."""

    case_id: str
    role: str
    hamiltonian: np.ndarray
    gauge_projection: np.ndarray
    casimir: np.ndarray
    window_lower: float
    window_upper: float
    window_predefined: bool
    window_fixed_before_hamiltonian_spectrum: bool
    window_fixed_before_commutator_measurement: bool
    window_definition_source: str
    predefinition_scope: str
    predefinition_certificate_id: str
    predefinition_certificate_status: str
    notes: str


@dataclass(frozen=True)
class AuditResult:
    """Serializable, fail-closed Casimir-corridor measurements."""

    case_id: str
    role: str
    window_lower: float
    window_upper: float
    window_predefined: bool
    window_fixed_before_hamiltonian_spectrum: bool
    window_fixed_before_commutator_measurement: bool
    window_definition_source: str
    predefinition_scope: str
    predefinition_certificate_id: str
    predefinition_certificate_status: str
    dimension: int
    physical_rank: int
    corridor_rank: int
    physical_exterior_rank: int
    casimir_spectral_separation: float
    hamiltonian_spectral_separation: float
    hamiltonian_self_adjoint_error: float
    casimir_self_adjoint_error: float
    gauge_projection_error: float
    hamiltonian_gauge_commutator: float
    casimir_gauge_commutator: float
    corridor_projection_error: float
    corridor_gauge_leakage: float
    hamiltonian_corridor_commutator: float
    commutator_tolerance: float
    definition_status: str
    gauge_casimir_structure_used: bool
    rg_coercivity_used: bool
    numerical_residual_used: bool
    symmetry_provenance_status: str
    gap_fingerprint_status: str
    transfer_decision: str
    claim_status: str
    notes: str


FIELDS = tuple(AuditResult.__dataclass_fields__)


def _finite_array(value: np.ndarray, *, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=np.complex128)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} contains non-finite values")
    return array


def _operator_error(operator: np.ndarray) -> float:
    return float(np.linalg.norm(operator - operator.conj().T, ord=2))


def _projection_error(projection: np.ndarray) -> float:
    adjoint_error = np.linalg.norm(projection - projection.conj().T, ord=2)
    idempotence_error = np.linalg.norm(projection @ projection - projection, ord=2)
    return float(max(adjoint_error, idempotence_error))


def _commutator_norm(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.linalg.norm(left @ right - right @ left, ord=2))


def _spectral_separation(
    operator: np.ndarray,
    corridor_basis: np.ndarray,
    exterior_basis: np.ndarray,
) -> float:
    """Return the distance between two reducing finite spectral blocks.

    A zero value is reported when either block is empty; the structural rank
    gates then reject that row before the separation can support a decision.
    """

    if corridor_basis.shape[1] == 0 or exterior_basis.shape[1] == 0:
        return 0.0
    corridor_spectrum = np.linalg.eigvalsh(corridor_basis.conj().T @ operator @ corridor_basis)
    exterior_spectrum = np.linalg.eigvalsh(exterior_basis.conj().T @ operator @ exterior_basis)
    distances = np.abs(corridor_spectrum[:, None] - exterior_spectrum[None, :])
    return float(np.min(distances))


def decide_transfer(
    case: AuditCase,
    *,
    physical_rank: int,
    corridor_rank: int,
    physical_exterior_rank: int,
    hamiltonian_gauge_commutator: float,
    casimir_gauge_commutator: float,
    corridor_projection_error: float,
    corridor_gauge_leakage: float,
    hamiltonian_corridor_commutator: float,
    casimir_spectral_separation: float,
    hamiltonian_spectral_separation: float,
    tolerance: float,
) -> str:
    """Return a non-circular, gauge-aware and dynamics-aware decision."""

    fixed_before_measurement = (
        case.window_predefined
        and case.window_fixed_before_hamiltonian_spectrum
        and case.window_fixed_before_commutator_measurement
    )
    if not fixed_before_measurement:
        return "reject_posthoc_or_adaptive_casimir_window"
    if physical_rank == 0:
        return "reject_empty_physical_gauge_sector"
    if hamiltonian_gauge_commutator > tolerance:
        return "reject_hamiltonian_breaks_gauge_sector"
    if casimir_gauge_commutator > tolerance:
        return "reject_casimir_not_gauge_compatible"
    if corridor_rank == 0:
        return "reject_empty_casimir_corridor"
    if physical_exterior_rank == 0:
        return "reject_no_physical_exterior_for_corridor_audit"
    if corridor_projection_error > tolerance or corridor_gauge_leakage > tolerance:
        return "reject_invalid_gauge_corridor_projection"
    if hamiltonian_corridor_commutator > tolerance:
        return "reject_hamiltonian_mixes_casimir_corridor"
    if casimir_spectral_separation <= tolerance or hamiltonian_spectral_separation <= tolerance:
        return "reject_corridor_spectral_separation_collapse"
    if case.predefinition_scope == "u1_strong_coupling_control":
        if (
            case.predefinition_certificate_status != "analytic_control_only"
            or not case.predefinition_certificate_id
        ):
            return "blocked_u1_control_provenance_not_verified"
        return "control_pass_u1_strong_coupling_no_yang_mills_claim"
    if (
        case.predefinition_certificate_status != "independently_verified"
        or not case.predefinition_certificate_id
    ):
        return "blocked_physical_predefinition_not_independently_verified"
    if case.predefinition_scope != "physical_yang_mills":
        return "blocked_nonphysical_predefinition_scope"
    return "eligible_for_analytic_review_no_claim"


def audit_case(case: AuditCase, *, tolerance: float = 1e-10) -> AuditResult:
    """Derive the physical Casimir window and audit its symmetry properties."""

    if not np.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("tolerance must be finite and strictly positive")
    if not np.isfinite(case.window_lower) or not np.isfinite(case.window_upper):
        raise ValueError("Casimir window bounds must be finite")
    if case.window_lower > case.window_upper:
        raise ValueError("window_lower must not exceed window_upper")
    if not case.window_definition_source:
        raise ValueError("window_definition_source must not be empty")
    allowed_scopes = {"synthetic_fixture", "u1_strong_coupling_control", "physical_yang_mills"}
    if case.predefinition_scope not in allowed_scopes:
        raise ValueError(f"unknown predefinition_scope: {case.predefinition_scope}")
    allowed_certificate_statuses = {
        "missing",
        "synthetic_fixture_only",
        "analytic_control_only",
        "independently_verified",
    }
    if case.predefinition_certificate_status not in allowed_certificate_statuses:
        raise ValueError(
            f"unknown predefinition_certificate_status: {case.predefinition_certificate_status}"
        )

    hamiltonian = _finite_array(case.hamiltonian, name="hamiltonian")
    gauge_projection = _finite_array(case.gauge_projection, name="gauge_projection")
    casimir = _finite_array(case.casimir, name="casimir")
    if hamiltonian.ndim != 2 or hamiltonian.shape[0] != hamiltonian.shape[1]:
        raise ValueError("hamiltonian must be square")
    dimension = hamiltonian.shape[0]
    if gauge_projection.shape != (dimension, dimension) or casimir.shape != (dimension, dimension):
        raise ValueError("gauge_projection and casimir must match the Hamiltonian dimension")

    hamiltonian_self_adjoint_error = _operator_error(hamiltonian)
    casimir_self_adjoint_error = _operator_error(casimir)
    gauge_projection_error = _projection_error(gauge_projection)
    if hamiltonian_self_adjoint_error > tolerance:
        raise ValueError("hamiltonian is not self-adjoint")
    if casimir_self_adjoint_error > tolerance:
        raise ValueError("casimir is not self-adjoint")
    if gauge_projection_error > tolerance:
        raise ValueError("gauge_projection is not an orthogonal projector")

    gauge_eigenvalues, gauge_eigenvectors = np.linalg.eigh(gauge_projection)
    physical_basis = gauge_eigenvectors[:, gauge_eigenvalues > 0.5]
    physical_rank = int(physical_basis.shape[1])
    corridor_basis = np.empty((dimension, 0), dtype=np.complex128)
    exterior_basis = np.empty((dimension, 0), dtype=np.complex128)
    if physical_rank:
        physical_casimir = physical_basis.conj().T @ casimir @ physical_basis
        casimir_eigenvalues, casimir_eigenvectors = np.linalg.eigh(physical_casimir)
        selected = (casimir_eigenvalues >= case.window_lower - tolerance) & (
            casimir_eigenvalues <= case.window_upper + tolerance
        )
        selected_vectors = casimir_eigenvectors[:, selected]
        corridor_basis = physical_basis @ selected_vectors
        exterior_basis = physical_basis @ casimir_eigenvectors[:, ~selected]
        corridor_projection = corridor_basis @ corridor_basis.conj().T
    else:
        corridor_projection = np.zeros_like(gauge_projection)

    corridor_rank = int(round(float(np.trace(corridor_projection).real)))
    physical_exterior_rank = physical_rank - corridor_rank
    hamiltonian_gauge_commutator = _commutator_norm(hamiltonian, gauge_projection)
    casimir_gauge_commutator = _commutator_norm(casimir, gauge_projection)
    corridor_projection_error = _projection_error(corridor_projection)
    identity = np.eye(dimension, dtype=np.complex128)
    corridor_gauge_leakage = float(
        np.linalg.norm((identity - gauge_projection) @ corridor_projection, ord=2)
    )
    hamiltonian_corridor_commutator = _commutator_norm(hamiltonian, corridor_projection)
    casimir_spectral_separation = _spectral_separation(casimir, corridor_basis, exterior_basis)
    hamiltonian_spectral_separation = _spectral_separation(
        hamiltonian,
        corridor_basis,
        exterior_basis,
    )
    definition_status = (
        "predefined_casimir_window"
        if (
            case.window_predefined
            and case.window_fixed_before_hamiltonian_spectrum
            and case.window_fixed_before_commutator_measurement
        )
        else "posthoc_or_adaptive_window"
    )
    transfer_decision = decide_transfer(
        case,
        physical_rank=physical_rank,
        corridor_rank=corridor_rank,
        physical_exterior_rank=physical_exterior_rank,
        hamiltonian_gauge_commutator=hamiltonian_gauge_commutator,
        casimir_gauge_commutator=casimir_gauge_commutator,
        corridor_projection_error=corridor_projection_error,
        corridor_gauge_leakage=corridor_gauge_leakage,
        hamiltonian_corridor_commutator=hamiltonian_corridor_commutator,
        casimir_spectral_separation=casimir_spectral_separation,
        hamiltonian_spectral_separation=hamiltonian_spectral_separation,
        tolerance=tolerance,
    )
    if (
        physical_rank > 0
        and corridor_rank > 0
        and physical_exterior_rank > 0
        and hamiltonian_gauge_commutator <= tolerance
        and casimir_gauge_commutator <= tolerance
        and hamiltonian_corridor_commutator <= tolerance
        and casimir_spectral_separation > tolerance
        and hamiltonian_spectral_separation > tolerance
    ):
        gap_fingerprint_status = "finite_reducing_corridor_gap_fingerprint"
    elif hamiltonian_corridor_commutator > tolerance:
        gap_fingerprint_status = "reject_nonreducing_corridor"
    elif corridor_rank > 0 and physical_exterior_rank > 0:
        gap_fingerprint_status = "reject_spectral_or_structural_gate"
    else:
        gap_fingerprint_status = "reject_degenerate_corridor"

    return AuditResult(
        case_id=case.case_id,
        role=case.role,
        window_lower=float(case.window_lower),
        window_upper=float(case.window_upper),
        window_predefined=case.window_predefined,
        window_fixed_before_hamiltonian_spectrum=case.window_fixed_before_hamiltonian_spectrum,
        window_fixed_before_commutator_measurement=case.window_fixed_before_commutator_measurement,
        window_definition_source=case.window_definition_source,
        predefinition_scope=case.predefinition_scope,
        predefinition_certificate_id=case.predefinition_certificate_id,
        predefinition_certificate_status=case.predefinition_certificate_status,
        dimension=dimension,
        physical_rank=physical_rank,
        corridor_rank=corridor_rank,
        physical_exterior_rank=physical_exterior_rank,
        casimir_spectral_separation=casimir_spectral_separation,
        hamiltonian_spectral_separation=hamiltonian_spectral_separation,
        hamiltonian_self_adjoint_error=hamiltonian_self_adjoint_error,
        casimir_self_adjoint_error=casimir_self_adjoint_error,
        gauge_projection_error=gauge_projection_error,
        hamiltonian_gauge_commutator=hamiltonian_gauge_commutator,
        casimir_gauge_commutator=casimir_gauge_commutator,
        corridor_projection_error=corridor_projection_error,
        corridor_gauge_leakage=corridor_gauge_leakage,
        hamiltonian_corridor_commutator=hamiltonian_corridor_commutator,
        commutator_tolerance=float(tolerance),
        definition_status=definition_status,
        gauge_casimir_structure_used=True,
        rg_coercivity_used=False,
        numerical_residual_used=False,
        symmetry_provenance_status="gauge_casimir_only_no_rg_or_residual_input",
        gap_fingerprint_status=gap_fingerprint_status,
        transfer_decision=transfer_decision,
        claim_status="diagnostic_only_no_yang_mills_claim",
        notes=case.notes,
    )


def _base_case(**overrides: object) -> AuditCase:
    values: dict[str, object] = {
        "case_id": "predefined_gauge_casimir_positive_control",
        "role": "positive_control",
        "hamiltonian": np.diag([0.0, 1.0, 4.0, 7.0, 8.0]),
        "gauge_projection": np.diag([1.0, 1.0, 1.0, 0.0, 0.0]),
        "casimir": np.diag([0.0, 0.75, 2.0, 0.75, 2.0]),
        "window_lower": 0.5,
        "window_upper": 1.0,
        "window_predefined": True,
        "window_fixed_before_hamiltonian_spectrum": True,
        "window_fixed_before_commutator_measurement": True,
        "window_definition_source": "synthetic_gauss_sector_casimir_interval",
        "predefinition_scope": "synthetic_fixture",
        "predefinition_certificate_id": "SYNTHETIC_CASIMIR_CORRIDOR_FIXTURE",
        "predefinition_certificate_status": "synthetic_fixture_only",
        "notes": "Exact finite gauge sector and reducing Casimir corridor; calibration only.",
    }
    values.update(overrides)
    return AuditCase(**values)


def build_matched_controls() -> list[AuditCase]:
    """Return two positives and five mechanism-matched negatives."""

    non_gauge_casimir = np.diag([0.0, 0.75, 2.0, 0.75, 2.0])
    non_gauge_casimir[1, 3] = non_gauge_casimir[3, 1] = 0.2
    gauge_breaking_hamiltonian = np.diag([0.0, 1.0, 4.0, 7.0, 8.0])
    gauge_breaking_hamiltonian[1, 3] = gauge_breaking_hamiltonian[3, 1] = 0.2
    corridor_mixing_hamiltonian = np.diag([0.0, 1.0, 4.0, 7.0, 8.0])
    corridor_mixing_hamiltonian[1, 2] = corridor_mixing_hamiltonian[2, 1] = 0.2
    collapsed_gap_hamiltonian = np.diag([0.0, 1.0, 1.0, 7.0, 8.0])
    u1_electric_casimir = np.diag([4.0, 1.0, 0.0, 1.0, 4.0])
    return [
        _base_case(),
        _base_case(
            case_id="u1_strong_coupling_loop_flux_positive_control",
            role="u1_strong_coupling_positive_control",
            hamiltonian=u1_electric_casimir,
            gauge_projection=np.eye(5),
            casimir=u1_electric_casimir,
            window_lower=0.5,
            window_upper=1.5,
            window_definition_source="analytic_compact_u1_loop_flux_electric_casimir_n_squared",
            predefinition_scope="u1_strong_coupling_control",
            predefinition_certificate_id="U1_LOOP_FLUX_ELECTRIC_CASIMIR_ANALYTIC",
            predefinition_certificate_status="analytic_control_only",
            notes=(
                "Gauge-reduced compact U(1) loop-flux basis n=-2,...,2 with electric energy "
                "n^2 in fixed units; the n=±1 corridor has a unit spectral separation."
            ),
        ),
        _base_case(
            case_id="negative_posthoc_window_same_operators",
            role="matched_negative_control",
            window_predefined=False,
            window_fixed_before_hamiltonian_spectrum=False,
            window_fixed_before_commutator_measurement=False,
            window_definition_source="selected_after_hamiltonian_and_commutator_inspection",
            predefinition_certificate_id="",
            predefinition_certificate_status="missing",
            notes="Numerically identical to the positive fixture, but the interval is post-hoc.",
        ),
        _base_case(
            case_id="negative_casimir_crosses_gauge_boundary",
            role="matched_negative_control",
            casimir=non_gauge_casimir,
            notes="The Casimir mixes physical and nonphysical vectors and cannot define a gauge symmetry sector.",
        ),
        _base_case(
            case_id="negative_hamiltonian_breaks_gauge_sector",
            role="matched_negative_control",
            hamiltonian=gauge_breaking_hamiltonian,
            notes="The Hamiltonian itself mixes the declared physical and nonphysical sectors.",
        ),
        _base_case(
            case_id="negative_gauge_invariant_but_corridor_mixed",
            role="matched_negative_control",
            hamiltonian=corridor_mixing_hamiltonian,
            notes="Gauge invariance holds, but the Hamiltonian couples the corridor to its physical exterior.",
        ),
        _base_case(
            case_id="negative_reducing_corridor_gap_collapsed",
            role="matched_negative_control",
            hamiltonian=collapsed_gap_hamiltonian,
            notes="The corridor reduces the Hamiltonian, but its energy is degenerate with the physical exterior.",
        ),
    ]


def build_results(*, tolerance: float = 1e-10) -> list[AuditResult]:
    return [audit_case(case, tolerance=tolerance) for case in build_matched_controls()]


def write_outputs(results: Sequence[AuditResult], output_dir: Path) -> dict[str, Path]:
    """Write deterministic CSV, JSON and Markdown audit artefacts."""

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

    decision_counts: dict[str, int] = {}
    for result in results:
        decision_counts[result.transfer_decision] = decision_counts.get(result.transfer_decision, 0) + 1
    payload = {
        "schema_version": 1,
        "date": DATE,
        "scope_status": "finite_synthetic_and_u1_strong_coupling_controls_only",
        "claim_pass": 0,
        "transfer_decision_counts": decision_counts,
        "rows": [asdict(result) for result in results],
    }
    paths["json"].write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    lines = [
        f"# Yang--Mills Casimir-Corridor Ledger {DATE}",
        "",
        "Status: finite synthetic gauge/Casimir calibration; no Yang--Mills or continuum claim.",
        "",
        "| case | phys rank | corridor rank | [H,P_g] | [C,P_g] | [H,P_I] | H separation | definition | decision |",
        "|---|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for result in results:
        lines.append(
            f"| {result.case_id} | {result.physical_rank} | {result.corridor_rank} | "
            f"{result.hamiltonian_gauge_commutator:.8g} | {result.casimir_gauge_commutator:.8g} | "
            f"{result.hamiltonian_corridor_commutator:.8g} | "
            f"{result.hamiltonian_spectral_separation:.8g} | {result.definition_status} | "
            f"{result.transfer_decision} |"
        )
    lines.extend(
        [
            "",
            "## Guardrails",
            "",
            "- `P_I` is derived from the declared Casimir interval inside `P_phys`; it is not supplied by the candidate.",
            "- The post-hoc control has the same matrices and interval as the positive fixture and is rejected only on provenance.",
            "- Gauge preservation (`[H,P_phys]=0`) does not imply Casimir-corridor preservation (`[H,P_I]=0`).",
            "- The compact U(1) loop-flux control has electric Casimir `n^2` and a unit first-flux gap fingerprint.",
            "- Corridor provenance uses gauge/Casimir structure only; RG coercivity and residual smallness are not inputs.",
            "- Even an exactly reducing corridor is rejected when its Hamiltonian separation from the physical exterior collapses.",
            "- A physical transfer still requires an independently verified Yang--Mills predefinition certificate.",
            "- The U(1) control passes only as a nonphysical analytic calibration; no bundled row is transfer-eligible.",
            "- No bundled row is a physical non-Abelian Yang--Mills certificate; `claim_pass` remains 0.",
            "",
        ]
    )
    paths["md"].write_text("\n".join(lines), encoding="utf-8")
    return paths


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Finite Yang--Mills gauge/Casimir corridor audit")
    parser.add_argument("--tolerance", type=float, default=1e-10)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "_results")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    results = build_results(tolerance=args.tolerance)
    paths = write_outputs(results, args.output_dir)
    for result in results:
        print(
            f"{result.case_id}: [H,P_g]={result.hamiltonian_gauge_commutator:.8g} "
            f"[C,P_g]={result.casimir_gauge_commutator:.8g} "
            f"[H,P_I]={result.hamiltonian_corridor_commutator:.8g} "
            f"decision={result.transfer_decision}"
        )
    print(f"claim_pass=0 json={paths['json']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
