"""Finite Hodge-inspired orbit-rigidity audit for Yang--Mills RG signals.

The audit asks whether good lattice/strong-coupling/RG signals travel through
one predeclared family of bounded-complexity gauge-orbit components.  It keeps
the four selector fields explicit:

``bounded_scale_component``
    A derived statement about component coherence and a registered complexity
    bound across all observed scales.
``gauge_orbit_certificate``
    A provenance certificate that must not be inferred from the successful gap
    signal itself.
``continuum_residue_gap``
    The worst registered separation from residual modes along the finite route.
``bad_component_escape``
    The largest measured escape into components outside the registered family.

This is a finite guardrail, not a continuum theorem.  Mosco convergence of
Dirichlet forms on varying Hilbert spaces and norm-resolvent spectral
convergence motivate separate limit-closure gates (Grothaus--Wittmann,
arXiv:2105.05140; Roesler, arXiv:1812.02525).  The script never infers those
gates from good-scale density, a finite gap floor, or a successful form test.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parent
DATE = "2026-08-26"
STEM = f"YM_ORBIT_RIGIDITY_LEDGER_{DATE}"


@dataclass(frozen=True)
class OrbitRigidityCase:
    """One finite scale-family and its provenance/closure declarations."""

    case_id: str
    role: str
    scale_component_ids: tuple[str, ...]
    component_complexities: tuple[int, ...]
    good_scale_mask: tuple[bool, ...]
    finite_gap_signals: tuple[float, ...]
    continuum_residue_gaps: tuple[float, ...]
    bad_component_escape_values: tuple[float, ...]
    component_complexity_bound: int
    good_scale_density_floor: float
    finite_gap_floor_required: float
    continuum_residue_gap_floor: float
    bad_component_escape_ceiling: float
    new_remainder_modes: int
    component_family_predefined: bool
    component_lineage_predefined: bool
    component_family_source: str
    gauge_orbit_certificate: str
    gauge_orbit_certificate_status: str
    gauge_orbit_certificate_independent_of_gap_signal: bool
    dirichlet_form_status: str
    mosco_liminf_status: str
    mosco_recovery_status: str
    compactness_status: str
    limit_closure_certificate_id: str
    limit_closure_certificate_status: str
    predefinition_scope: str
    notes: str


@dataclass(frozen=True)
class OrbitRigidityResult:
    """Serializable finite measurements and fail-closed transfer decision."""

    case_id: str
    role: str
    n_scales: int
    component_family_predefined: bool
    component_lineage_predefined: bool
    component_family_source: str
    component_family_coherent: bool
    component_complexity_bound: int
    observed_max_component_complexity: int
    complexity_overrun: int
    bounded_scale_component: bool
    good_scale_count: int
    good_scale_density: float
    good_scale_density_floor: float
    finite_gap_floor: float
    finite_gap_floor_required: float
    gauge_orbit_certificate: str
    gauge_orbit_certificate_status: str
    gauge_orbit_certificate_independent_of_gap_signal: bool
    dirichlet_form_status: str
    continuum_residue_gap: float
    continuum_residue_gap_floor: float
    bad_component_escape: float
    bad_component_escape_total: float
    bad_component_escape_ceiling: float
    new_remainder_modes: int
    mosco_liminf_status: str
    mosco_recovery_status: str
    compactness_status: str
    limit_closure_certificate_id: str
    limit_closure_certificate_status: str
    predefinition_scope: str
    form_test_status: str
    transfer_decision: str
    claim_status: str
    notes: str


FIELDS = tuple(OrbitRigidityResult.__dataclass_fields__)


def _finite_nonnegative(values: Sequence[float], *, name: str) -> tuple[float, ...]:
    converted = tuple(float(value) for value in values)
    if not converted:
        raise ValueError(f"{name} must not be empty")
    if not all(math.isfinite(value) and value >= 0.0 for value in converted):
        raise ValueError(f"{name} must contain finite nonnegative values")
    return converted


def _validate_case(case: OrbitRigidityCase) -> None:
    lengths = {
        len(case.scale_component_ids),
        len(case.component_complexities),
        len(case.good_scale_mask),
        len(case.finite_gap_signals),
        len(case.continuum_residue_gaps),
        len(case.bad_component_escape_values),
    }
    if lengths == {0} or len(lengths) != 1:
        raise ValueError("all scale sequences must have the same nonzero length")
    if not all(identifier for identifier in case.scale_component_ids):
        raise ValueError("scale_component_ids must not contain empty identifiers")
    if not all(isinstance(value, int) and not isinstance(value, bool) and value >= 1 for value in case.component_complexities):
        raise ValueError("component_complexities must be positive integers")
    if not all(isinstance(value, bool) for value in case.good_scale_mask):
        raise ValueError("good_scale_mask must contain booleans")
    _finite_nonnegative(case.finite_gap_signals, name="finite_gap_signals")
    _finite_nonnegative(case.continuum_residue_gaps, name="continuum_residue_gaps")
    _finite_nonnegative(case.bad_component_escape_values, name="bad_component_escape_values")
    if not isinstance(case.component_complexity_bound, int) or case.component_complexity_bound < 1:
        raise ValueError("component_complexity_bound must be a positive integer")
    if not 0.0 <= case.good_scale_density_floor <= 1.0:
        raise ValueError("good_scale_density_floor must lie in [0, 1]")
    if not math.isfinite(case.finite_gap_floor_required) or case.finite_gap_floor_required <= 0.0:
        raise ValueError("finite_gap_floor_required must be finite and strictly positive")
    if not math.isfinite(case.continuum_residue_gap_floor) or case.continuum_residue_gap_floor <= 0.0:
        raise ValueError("continuum_residue_gap_floor must be finite and strictly positive")
    if not math.isfinite(case.bad_component_escape_ceiling) or case.bad_component_escape_ceiling < 0.0:
        raise ValueError("bad_component_escape_ceiling must be finite and nonnegative")
    if not isinstance(case.new_remainder_modes, int) or case.new_remainder_modes < 0:
        raise ValueError("new_remainder_modes must be a nonnegative integer")
    if not case.component_family_source:
        raise ValueError("component_family_source must not be empty")

    allowed_gauge_statuses = {
        "missing",
        "derived_from_gap_signal",
        "analytic_control_only",
        "independently_verified",
    }
    if case.gauge_orbit_certificate_status not in allowed_gauge_statuses:
        raise ValueError(f"unknown gauge_orbit_certificate_status: {case.gauge_orbit_certificate_status}")
    allowed_form_statuses = {"missing", "finite_fixture_pass", "independently_verified"}
    if case.dirichlet_form_status not in allowed_form_statuses:
        raise ValueError(f"unknown dirichlet_form_status: {case.dirichlet_form_status}")
    allowed_limit_statuses = {"missing", "finite_fixture_pass", "independently_verified"}
    for field_name in ("mosco_liminf_status", "mosco_recovery_status", "compactness_status"):
        value = getattr(case, field_name)
        if value not in allowed_limit_statuses:
            raise ValueError(f"unknown {field_name}: {value}")
    allowed_closure_statuses = {"missing", "finite_fixture_only", "independently_verified"}
    if case.limit_closure_certificate_status not in allowed_closure_statuses:
        raise ValueError(
            f"unknown limit_closure_certificate_status: {case.limit_closure_certificate_status}"
        )
    if case.predefinition_scope not in {"finite_control", "physical_yang_mills"}:
        raise ValueError(f"unknown predefinition_scope: {case.predefinition_scope}")


def _form_test_status(
    case: OrbitRigidityCase,
    *,
    bounded_scale_component: bool,
    good_scale_density: float,
    finite_gap_floor: float,
    continuum_residue_gap: float,
    bad_component_escape: float,
) -> str:
    if not case.component_family_predefined or not case.component_lineage_predefined:
        return "finite_form_test_reject_posthoc_family"
    if not bounded_scale_component:
        return "finite_form_test_reject_unbounded_component"
    if (
        good_scale_density < case.good_scale_density_floor
        or finite_gap_floor < case.finite_gap_floor_required
    ):
        return "finite_form_test_reject_signal_floor"
    if case.dirichlet_form_status == "missing":
        return "finite_form_test_reject_missing_dirichlet_structure"
    if continuum_residue_gap < case.continuum_residue_gap_floor:
        return "finite_form_test_reject_residue_gap"
    if bad_component_escape > case.bad_component_escape_ceiling or case.new_remainder_modes > 0:
        return "finite_form_test_reject_component_escape"
    return "finite_form_test_pass"


def decide_transfer(
    case: OrbitRigidityCase,
    *,
    bounded_scale_component: bool,
    good_scale_density: float,
    finite_gap_floor: float,
    continuum_residue_gap: float,
    bad_component_escape: float,
) -> str:
    """Separate a correct finite form test from origin and limit closure."""

    if not case.component_family_predefined or not case.component_lineage_predefined:
        return "reject_posthoc_or_adaptive_component_family"
    if not bounded_scale_component:
        return "reject_unbounded_or_incoherent_scale_component"
    if (
        good_scale_density < case.good_scale_density_floor
        or finite_gap_floor < case.finite_gap_floor_required
    ):
        return "reject_insufficient_finite_gap_signal"
    if case.gauge_orbit_certificate_status == "missing" or not case.gauge_orbit_certificate:
        return "diagnostic_only_missing_independent_gauge_orbit_certificate"
    if (
        case.gauge_orbit_certificate_status == "derived_from_gap_signal"
        or not case.gauge_orbit_certificate_independent_of_gap_signal
    ):
        return "reject_gauge_orbit_certificate_derived_from_gap_signal"
    if case.dirichlet_form_status == "missing":
        return "diagnostic_only_missing_dirichlet_form_structure"
    if continuum_residue_gap < case.continuum_residue_gap_floor:
        return "reject_continuum_residue_gap_collapse"
    if bad_component_escape > case.bad_component_escape_ceiling or case.new_remainder_modes > 0:
        return "reject_bad_component_escape_or_new_remainder_modes"
    if case.predefinition_scope == "finite_control":
        return "form_test_pass_without_physical_origin_or_limit_closure"
    if case.gauge_orbit_certificate_status != "independently_verified":
        return "diagnostic_only_missing_independent_gauge_orbit_certificate"
    limit_statuses = (
        case.mosco_liminf_status,
        case.mosco_recovery_status,
        case.compactness_status,
    )
    if (
        any(status != "independently_verified" for status in limit_statuses)
        or case.limit_closure_certificate_status != "independently_verified"
        or not case.limit_closure_certificate_id
    ):
        return "diagnostic_only_missing_limit_closure_certificate"
    return "eligible_for_analytic_review_no_claim"


def audit_case(case: OrbitRigidityCase) -> OrbitRigidityResult:
    """Compute the four selector fields and a claim-neutral decision."""

    _validate_case(case)
    gaps = _finite_nonnegative(case.finite_gap_signals, name="finite_gap_signals")
    residue_gaps = _finite_nonnegative(
        case.continuum_residue_gaps,
        name="continuum_residue_gaps",
    )
    escapes = _finite_nonnegative(
        case.bad_component_escape_values,
        name="bad_component_escape_values",
    )
    n_scales = len(case.scale_component_ids)
    component_family_coherent = len(set(case.scale_component_ids)) == 1
    observed_max_component_complexity = max(case.component_complexities)
    complexity_overrun = max(0, observed_max_component_complexity - case.component_complexity_bound)
    bounded_scale_component = component_family_coherent and complexity_overrun == 0
    good_scale_count = sum(case.good_scale_mask)
    good_scale_density = good_scale_count / n_scales
    finite_gap_floor = min(gaps)
    continuum_residue_gap = min(residue_gaps)
    bad_component_escape = max(escapes)
    bad_component_escape_total = sum(escapes)
    form_test_status = _form_test_status(
        case,
        bounded_scale_component=bounded_scale_component,
        good_scale_density=good_scale_density,
        finite_gap_floor=finite_gap_floor,
        continuum_residue_gap=continuum_residue_gap,
        bad_component_escape=bad_component_escape,
    )
    transfer_decision = decide_transfer(
        case,
        bounded_scale_component=bounded_scale_component,
        good_scale_density=good_scale_density,
        finite_gap_floor=finite_gap_floor,
        continuum_residue_gap=continuum_residue_gap,
        bad_component_escape=bad_component_escape,
    )
    return OrbitRigidityResult(
        case_id=case.case_id,
        role=case.role,
        n_scales=n_scales,
        component_family_predefined=case.component_family_predefined,
        component_lineage_predefined=case.component_lineage_predefined,
        component_family_source=case.component_family_source,
        component_family_coherent=component_family_coherent,
        component_complexity_bound=case.component_complexity_bound,
        observed_max_component_complexity=observed_max_component_complexity,
        complexity_overrun=complexity_overrun,
        bounded_scale_component=bounded_scale_component,
        good_scale_count=good_scale_count,
        good_scale_density=good_scale_density,
        good_scale_density_floor=float(case.good_scale_density_floor),
        finite_gap_floor=finite_gap_floor,
        finite_gap_floor_required=float(case.finite_gap_floor_required),
        gauge_orbit_certificate=case.gauge_orbit_certificate,
        gauge_orbit_certificate_status=case.gauge_orbit_certificate_status,
        gauge_orbit_certificate_independent_of_gap_signal=(
            case.gauge_orbit_certificate_independent_of_gap_signal
        ),
        dirichlet_form_status=case.dirichlet_form_status,
        continuum_residue_gap=continuum_residue_gap,
        continuum_residue_gap_floor=float(case.continuum_residue_gap_floor),
        bad_component_escape=bad_component_escape,
        bad_component_escape_total=bad_component_escape_total,
        bad_component_escape_ceiling=float(case.bad_component_escape_ceiling),
        new_remainder_modes=case.new_remainder_modes,
        mosco_liminf_status=case.mosco_liminf_status,
        mosco_recovery_status=case.mosco_recovery_status,
        compactness_status=case.compactness_status,
        limit_closure_certificate_id=case.limit_closure_certificate_id,
        limit_closure_certificate_status=case.limit_closure_certificate_status,
        predefinition_scope=case.predefinition_scope,
        form_test_status=form_test_status,
        transfer_decision=transfer_decision,
        claim_status="diagnostic_only_no_yang_mills_claim",
        notes=case.notes,
    )


def _base_case(**overrides: object) -> OrbitRigidityCase:
    values: dict[str, object] = {
        "case_id": "finite_bounded_orbit_family_control",
        "role": "positive_control",
        "scale_component_ids": ("orbit_A",) * 5,
        "component_complexities": (2, 2, 3, 3, 3),
        "good_scale_mask": (True, True, False, True, True),
        "finite_gap_signals": (0.60, 0.56, 0.52, 0.48, 0.44),
        "continuum_residue_gaps": (0.35, 0.33, 0.31, 0.29, 0.27),
        "bad_component_escape_values": (0.01, 0.02, 0.02, 0.03, 0.04),
        "component_complexity_bound": 3,
        "good_scale_density_floor": 0.75,
        "finite_gap_floor_required": 0.40,
        "continuum_residue_gap_floor": 0.25,
        "bad_component_escape_ceiling": 0.05,
        "new_remainder_modes": 0,
        "component_family_predefined": True,
        "component_lineage_predefined": True,
        "component_family_source": "finite_fixture_predeclared_orbit_A",
        "gauge_orbit_certificate": "FINITE_FIXTURE_GAUSS_ORBIT_A",
        "gauge_orbit_certificate_status": "analytic_control_only",
        "gauge_orbit_certificate_independent_of_gap_signal": True,
        "dirichlet_form_status": "finite_fixture_pass",
        "mosco_liminf_status": "finite_fixture_pass",
        "mosco_recovery_status": "finite_fixture_pass",
        "compactness_status": "finite_fixture_pass",
        "limit_closure_certificate_id": "FINITE_FIXTURE_NO_CONTINUUM_CERTIFICATE",
        "limit_closure_certificate_status": "finite_fixture_only",
        "predefinition_scope": "finite_control",
        "notes": "Bounded finite family calibrating the four audit fields; no physical limit closure.",
    }
    values.update(overrides)
    return OrbitRigidityCase(**values)


def build_matched_controls() -> list[OrbitRigidityCase]:
    """Return one finite positive and seven mechanism-matched negatives."""

    return [
        _base_case(),
        _base_case(
            case_id="negative_high_good_density_unbounded_complexity",
            role="matched_negative_control",
            good_scale_mask=(True,) * 5,
            component_complexities=(2, 3, 5, 8, 13),
            notes="Perfect good-scale density cannot pay for unbounded component complexity.",
        ),
        _base_case(
            case_id="negative_stable_gap_missing_gauge_origin",
            role="matched_negative_control",
            gauge_orbit_certificate="",
            gauge_orbit_certificate_status="missing",
            notes="The finite form test passes, but no independent gauge-orbit origin is supplied.",
        ),
        _base_case(
            case_id="negative_gauge_certificate_derived_from_gap",
            role="matched_negative_control",
            gauge_orbit_certificate="SELECTED_FROM_SUCCESSFUL_GAP_SIGNAL",
            gauge_orbit_certificate_status="derived_from_gap_signal",
            gauge_orbit_certificate_independent_of_gap_signal=False,
            notes="The orbit label is reconstructed from the successful signal and is circular.",
        ),
        _base_case(
            case_id="negative_continuum_residue_gap_collapse",
            role="matched_negative_control",
            continuum_residue_gaps=(0.35, 0.24, 0.12, 0.03, 0.0),
            notes="Finite gaps remain good while the registered residue separation collapses.",
        ),
        _base_case(
            case_id="negative_bad_component_escape_new_modes",
            role="matched_negative_control",
            bad_component_escape_values=(0.01, 0.02, 0.08, 0.21, 0.46),
            new_remainder_modes=2,
            notes="New remainder modes and an escaping bad component block limit closure.",
        ),
        _base_case(
            case_id="negative_posthoc_component_family_same_numbers",
            role="matched_negative_control",
            component_family_predefined=False,
            component_lineage_predefined=False,
            component_family_source="selected_after_gap_and_residue_inspection",
            notes="Numerically identical to the positive fixture, but the family is selected post-hoc.",
        ),
        _base_case(
            case_id="negative_physical_label_missing_limit_closure",
            role="matched_negative_control",
            predefinition_scope="physical_yang_mills",
            gauge_orbit_certificate="DECLARED_PHYSICAL_GAUGE_ORBIT",
            gauge_orbit_certificate_status="independently_verified",
            dirichlet_form_status="independently_verified",
            mosco_liminf_status="missing",
            mosco_recovery_status="missing",
            compactness_status="missing",
            limit_closure_certificate_id="",
            limit_closure_certificate_status="missing",
            notes="A physical label and finite form pass do not replace Mosco/compact limit closure.",
        ),
    ]


def build_results() -> list[OrbitRigidityResult]:
    return [audit_case(case) for case in build_matched_controls()]


def write_outputs(results: Sequence[OrbitRigidityResult], output_dir: Path) -> dict[str, Path]:
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
        "scope_status": "finite_synthetic_orbit_rigidity_controls_only",
        "claim_pass": 0,
        "transfer_decision_counts": decision_counts,
        "rows": [asdict(result) for result in results],
    }
    paths["json"].write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    lines = [
        f"# Yang--Mills Orbit-Rigidity Ledger {DATE}",
        "",
        "Status: finite Hodge-inspired guardrail; no physical continuum or mass-gap claim.",
        "",
        "| case | good density | bounded_scale_component | gauge_orbit_certificate | continuum_residue_gap | bad_component_escape | form status | decision |",
        "|---|---:|---|---|---:|---:|---|---|",
    ]
    for result in results:
        certificate = result.gauge_orbit_certificate or "missing"
        lines.append(
            f"| {result.case_id} | {result.good_scale_density:.6g} | "
            f"{result.bounded_scale_component} | {certificate} | "
            f"{result.continuum_residue_gap:.6g} | {result.bad_component_escape:.6g} | "
            f"{result.form_test_status} | {result.transfer_decision} |"
        )
    lines.extend(
        [
            "",
            "## Guardrails",
            "",
            "- High good-scale density or a stable finite gap is insufficient without a bounded coherent component family.",
            "- `gauge_orbit_certificate` must be independent of the successful gap signal.",
            "- `continuum_residue_gap` is a finite route diagnostic, not a proof of continuum separation.",
            "- `bad_component_escape` and `new_remainder_modes` fail closed instead of being absorbed into an average.",
            "- Mosco liminf/recovery, compactness, and norm-resolvent-style closure remain separate certificate fields.",
            "- A correct finite form test without independent origin/limit closure remains diagnostic only.",
            "- No bundled row is physically transfer-eligible; `claim_pass` remains 0.",
            "",
        ]
    )
    paths["md"].write_text("\n".join(lines), encoding="utf-8")
    return paths


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Finite Yang--Mills orbit-rigidity/limit-closure audit")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "_results")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    results = build_results()
    paths = write_outputs(results, args.output_dir)
    for result in results:
        print(
            f"{result.case_id}: bounded={result.bounded_scale_component} "
            f"residue_gap={result.continuum_residue_gap:.6g} "
            f"escape={result.bad_component_escape:.6g} decision={result.transfer_decision}"
        )
    print(f"claim_pass=0 json={paths['json']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
