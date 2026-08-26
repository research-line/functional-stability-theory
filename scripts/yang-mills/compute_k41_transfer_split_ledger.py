"""K41-style local-positive-control versus continuum-transfer audit.

The K41 donor contributes a method, not a turbulence claim: a clean local or
static variational result must stay separate from the harder dynamical bridge.
Here the same wording discipline is applied to Yang--Mills.  A local
strong-coupling/lattice gap can be a positive control while the OS/continuum
step remains an explicitly open transfer hypothesis.

All bundled rows are finite or project-scope audit controls.  The ledger does
not prove a continuum Yang--Mills mass gap and transfers no turbulence claim.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Sequence


DATE = "2026-08-26"
STEM = f"YANG_MILLS_K41_TRANSFER_SPLIT_{DATE}"
ROOT = Path(__file__).resolve().parent


@dataclass(frozen=True)
class K41SplitCase:
    case_id: str
    role: str
    local_statement_id: str
    local_statement_scope: str
    local_proof_status: str
    local_gap_lower_bound: float
    budget_source_id: str
    evidence_scope: str
    claim_mode: str
    envelope_predefined: bool
    state_correlation_ratios: tuple[float, ...]
    envelope_min: float
    envelope_max: float
    defect_budget_terms: tuple[float, ...]
    defect_budget_ceiling: float
    warm_corridor_flags: tuple[bool, ...]
    rp_status: str
    os_reconstruction_status: str
    continuum_limit_status: str
    physical_normalization_status: str
    bridge_source_status: str
    notes: str


@dataclass(frozen=True)
class K41SplitResult:
    case_id: str
    role: str
    local_statement_id: str
    local_statement_scope: str
    local_proof_status: str
    local_gap_lower_bound: float
    budget_source_id: str
    evidence_scope: str
    claim_mode: str
    unconditional_positive_control: bool
    positive_control_status: str
    envelope_predefined: bool
    ratio_min: float
    ratio_max: float
    envelope_min: float
    envelope_max: float
    state_correlation_envelope: str
    defect_budget_total: float
    defect_budget_ceiling: float
    defect_budget_status: str
    warm_corridor_count: int
    warm_corridor_status: str
    rp_status: str
    os_reconstruction_status: str
    continuum_limit_status: str
    physical_normalization_status: str
    bridge_source_status: str
    transfer_hypothesis: str
    k41_style_split: bool
    split_status: str
    transfer_decision: str
    claim_status: str
    notes: str


FIELDS = tuple(K41SplitResult.__dataclass_fields__)

LOCAL_SCOPES = {
    "finite_fixture",
    "strong_coupling_lattice_local",
    "continuum_r4",
}
LOCAL_PROOF_STATUSES = {
    "missing",
    "assumed",
    "finite_control",
    "project_theorem",
    "independently_verified",
}
EVIDENCE_SCOPES = {
    "finite_fixture",
    "project_local_theorem",
    "physical_yang_mills",
}
CLAIM_MODES = {
    "local_positive_control",
    "transfer_hypothesis",
    "claimed_continuum_transfer",
}
RP_STATUSES = {
    "missing",
    "finite_control",
    "lattice_verified",
    "independently_verified",
}
BRIDGE_STATUSES = {
    "missing",
    "conditional",
    "finite_control",
    "independently_verified",
}
BRIDGE_SOURCE_STATUSES = {
    "missing",
    "independent",
    "derived_from_local_gap_signal",
}


def _finite_values(values: Sequence[float], *, name: str, positive: bool) -> tuple[float, ...]:
    if not values:
        raise ValueError(f"{name} must not be empty")
    converted = tuple(float(value) for value in values)
    if not all(math.isfinite(value) for value in converted):
        raise ValueError(f"{name} must contain only finite values")
    if positive and any(value <= 0.0 for value in converted):
        raise ValueError(f"{name} must contain strictly positive values")
    if not positive and any(value < 0.0 for value in converted):
        raise ValueError(f"{name} must contain nonnegative values")
    return converted


def _validate_case(case: K41SplitCase) -> None:
    if not case.case_id or not case.local_statement_id:
        raise ValueError("case_id and local_statement_id must not be empty")
    if not case.budget_source_id:
        raise ValueError("budget_source_id must not be empty")
    if case.local_statement_scope not in LOCAL_SCOPES:
        raise ValueError(f"unknown local_statement_scope: {case.local_statement_scope}")
    if case.local_proof_status not in LOCAL_PROOF_STATUSES:
        raise ValueError(f"unknown local_proof_status: {case.local_proof_status}")
    if case.evidence_scope not in EVIDENCE_SCOPES:
        raise ValueError(f"unknown evidence_scope: {case.evidence_scope}")
    if case.claim_mode not in CLAIM_MODES:
        raise ValueError(f"unknown claim_mode: {case.claim_mode}")
    if case.rp_status not in RP_STATUSES:
        raise ValueError(f"unknown rp_status: {case.rp_status}")
    for field_name in (
        "os_reconstruction_status",
        "continuum_limit_status",
        "physical_normalization_status",
    ):
        value = getattr(case, field_name)
        if value not in BRIDGE_STATUSES:
            raise ValueError(f"unknown {field_name}: {value}")
    if case.bridge_source_status not in BRIDGE_SOURCE_STATUSES:
        raise ValueError(f"unknown bridge_source_status: {case.bridge_source_status}")
    if not math.isfinite(case.local_gap_lower_bound) or case.local_gap_lower_bound < 0.0:
        raise ValueError("local_gap_lower_bound must be finite and nonnegative")
    if not math.isfinite(case.envelope_min) or case.envelope_min <= 0.0:
        raise ValueError("envelope_min must be finite and strictly positive")
    if not math.isfinite(case.envelope_max) or case.envelope_max < case.envelope_min:
        raise ValueError("envelope_max must be finite and at least envelope_min")
    if not math.isfinite(case.defect_budget_ceiling) or case.defect_budget_ceiling < 0.0:
        raise ValueError("defect_budget_ceiling must be finite and nonnegative")
    ratios = _finite_values(
        case.state_correlation_ratios,
        name="state_correlation_ratios",
        positive=True,
    )
    _finite_values(case.defect_budget_terms, name="defect_budget_terms", positive=False)
    if len(case.warm_corridor_flags) != len(ratios):
        raise ValueError("warm_corridor_flags must match state_correlation_ratios")
    if not all(isinstance(value, bool) for value in case.warm_corridor_flags):
        raise ValueError("warm_corridor_flags must contain booleans")


def _positive_control(case: K41SplitCase) -> tuple[bool, str]:
    local_scope = case.local_statement_scope in {
        "finite_fixture",
        "strong_coupling_lattice_local",
    }
    proof_ok = case.local_proof_status in {
        "finite_control",
        "project_theorem",
        "independently_verified",
    }
    if not local_scope:
        return False, "not_a_local_or_finite_statement"
    if not proof_ok:
        return False, "missing_or_assumed_local_result"
    if case.local_gap_lower_bound <= 0.0:
        return False, "nonpositive_local_gap_lower_bound"
    return True, "local_positive_control_pass"


def _transfer_hypothesis(case: K41SplitCase) -> str:
    if case.claim_mode == "local_positive_control":
        return "not_invoked_local_result_only"
    if case.claim_mode == "transfer_hypothesis":
        return "open_os_continuum_transfer_hypothesis"
    return "asserted_continuum_transfer"


def _all_physical_bridge_fields_verified(case: K41SplitCase) -> bool:
    return (
        case.evidence_scope == "physical_yang_mills"
        and case.rp_status == "independently_verified"
        and case.os_reconstruction_status == "independently_verified"
        and case.continuum_limit_status == "independently_verified"
        and case.physical_normalization_status == "independently_verified"
        and case.bridge_source_status == "independent"
    )


def _split_status(case: K41SplitCase, *, physical_bridge_verified: bool) -> tuple[bool, str]:
    if case.claim_mode == "local_positive_control":
        return True, "split_respected_local_positive_control"
    if case.claim_mode == "transfer_hypothesis":
        if physical_bridge_verified:
            return True, "split_respected_bridge_ready_for_analytic_review"
        return True, "split_respected_open_transfer_hypothesis"
    if physical_bridge_verified:
        return True, "split_respected_but_claim_requires_external_review"
    return False, "split_violated_local_result_relabelled_as_continuum"


def decide_transfer(
    case: K41SplitCase,
    *,
    unconditional_positive_control: bool,
    state_correlation_envelope: str,
    defect_budget_status: str,
    warm_corridor_count: int,
    physical_bridge_verified: bool,
) -> str:
    """Fail closed while keeping a valid local theorem visible."""

    if not unconditional_positive_control:
        return "reject_missing_local_positive_control"
    if not case.envelope_predefined:
        return "reject_posthoc_state_correlation_envelope"
    if warm_corridor_count:
        return "reject_warm_corridor_outside_local_control"
    if state_correlation_envelope != "within_predefined_envelope":
        return "reject_state_correlation_envelope_escape"
    if defect_budget_status != "within_predefined_defect_budget":
        return "reject_cross_scale_defect_budget_overrun"
    if case.bridge_source_status == "derived_from_local_gap_signal":
        return "reject_circular_bridge_derived_from_local_gap"
    if case.claim_mode == "local_positive_control":
        return "positive_control_pass_local_lattice_no_continuum_claim"
    if case.claim_mode == "claimed_continuum_transfer" and not physical_bridge_verified:
        return "reject_scope_inflation_local_result_relabelled_as_continuum"
    if physical_bridge_verified:
        return "eligible_for_analytic_review_no_claim"
    if case.evidence_scope != "physical_yang_mills":
        return "diagnostic_only_no_physical_yang_mills_bridge_evidence"
    return "diagnostic_only_open_os_continuum_transfer_hypothesis"


def audit_case(case: K41SplitCase) -> K41SplitResult:
    """Audit one proposed local-result/continuum-transfer split."""

    _validate_case(case)
    ratios = tuple(float(value) for value in case.state_correlation_ratios)
    defects = tuple(float(value) for value in case.defect_budget_terms)
    ratio_min = min(ratios)
    ratio_max = max(ratios)
    within_envelope = ratio_min >= case.envelope_min and ratio_max <= case.envelope_max
    if not case.envelope_predefined:
        envelope_status = "posthoc_envelope"
    elif within_envelope:
        envelope_status = "within_predefined_envelope"
    else:
        envelope_status = "outside_predefined_envelope"

    defect_budget_total = math.fsum(defects)
    defect_budget_status = (
        "within_predefined_defect_budget"
        if defect_budget_total <= case.defect_budget_ceiling + 1e-12
        else "predefined_defect_budget_overrun"
    )
    warm_corridor_count = sum(case.warm_corridor_flags)
    warm_corridor_status = (
        "warm_corridor_clear" if warm_corridor_count == 0 else "warm_corridor_detected"
    )
    unconditional_positive_control, positive_control_status = _positive_control(case)
    physical_bridge_verified = _all_physical_bridge_fields_verified(case)
    k41_style_split, split_status = _split_status(
        case,
        physical_bridge_verified=physical_bridge_verified,
    )
    transfer_decision = decide_transfer(
        case,
        unconditional_positive_control=unconditional_positive_control,
        state_correlation_envelope=envelope_status,
        defect_budget_status=defect_budget_status,
        warm_corridor_count=warm_corridor_count,
        physical_bridge_verified=physical_bridge_verified,
    )
    return K41SplitResult(
        case_id=case.case_id,
        role=case.role,
        local_statement_id=case.local_statement_id,
        local_statement_scope=case.local_statement_scope,
        local_proof_status=case.local_proof_status,
        local_gap_lower_bound=float(case.local_gap_lower_bound),
        budget_source_id=case.budget_source_id,
        evidence_scope=case.evidence_scope,
        claim_mode=case.claim_mode,
        unconditional_positive_control=unconditional_positive_control,
        positive_control_status=positive_control_status,
        envelope_predefined=case.envelope_predefined,
        ratio_min=ratio_min,
        ratio_max=ratio_max,
        envelope_min=float(case.envelope_min),
        envelope_max=float(case.envelope_max),
        state_correlation_envelope=envelope_status,
        defect_budget_total=defect_budget_total,
        defect_budget_ceiling=float(case.defect_budget_ceiling),
        defect_budget_status=defect_budget_status,
        warm_corridor_count=warm_corridor_count,
        warm_corridor_status=warm_corridor_status,
        rp_status=case.rp_status,
        os_reconstruction_status=case.os_reconstruction_status,
        continuum_limit_status=case.continuum_limit_status,
        physical_normalization_status=case.physical_normalization_status,
        bridge_source_status=case.bridge_source_status,
        transfer_hypothesis=_transfer_hypothesis(case),
        k41_style_split=k41_style_split,
        split_status=split_status,
        transfer_decision=transfer_decision,
        claim_status="diagnostic_only_no_yang_mills_or_turbulence_claim",
        notes=case.notes,
    )


@lru_cache(maxsize=1)
def _rg_positive_fixture() -> dict[str, object]:
    """Read the existing finite RG control instead of inventing a second baseline."""

    source_path = ROOT / "compute_os_capacity_ledger.py"
    spec = importlib.util.spec_from_file_location("ym_os_capacity_k41_source", source_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load RG budget source: {source_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    rows = [
        row
        for row in module.generated_v2_control_rows(k_max=10)
        if row["scenario"] == "strong_coupling_positive_control"
    ]
    if len(rows) != 10:
        raise RuntimeError("expected ten strong-coupling RG control rows")
    summary = module.summarize_v2_group("strong_coupling_positive_control", rows)
    if summary["decision"] != "control_pass_summable_no_claim":
        raise RuntimeError("existing strong-coupling RG control no longer passes its local gate")

    contractions = tuple(float(row["tau_B"]) for row in rows)
    geometric_mean = math.exp(math.fsum(math.log(value) for value in contractions) / len(rows))
    return {
        "local_gap_lower_bound": 1.0 - max(contractions),
        "state_correlation_ratios": tuple(value / geometric_mean for value in contractions),
        "defect_budget_terms": tuple(
            float(row["epsilon_safe"]) + float(row["eta_os_danger"])
            for row in rows
        ),
        "warm_corridor_flags": tuple(bool(module.is_bad_v2_row(row)) for row in rows),
    }


def _base_case(**overrides: object) -> K41SplitCase:
    fixture = _rg_positive_fixture()
    values: dict[str, object] = {
        "case_id": "strong_coupling_local_gap_positive_control",
        "role": "positive_control",
        "local_statement_id": "YM_LOCAL_STRONG_COUPLING_GAP",
        "local_statement_scope": "strong_coupling_lattice_local",
        "local_proof_status": "project_theorem",
        "local_gap_lower_bound": fixture["local_gap_lower_bound"],
        "budget_source_id": "compute_os_capacity_ledger.py:strong_coupling_positive_control:k_max=10",
        "evidence_scope": "project_local_theorem",
        "claim_mode": "local_positive_control",
        "envelope_predefined": True,
        "state_correlation_ratios": fixture["state_correlation_ratios"],
        "envelope_min": 0.98,
        "envelope_max": 1.02,
        "defect_budget_terms": fixture["defect_budget_terms"],
        "defect_budget_ceiling": 0.02,
        "warm_corridor_flags": fixture["warm_corridor_flags"],
        "rp_status": "lattice_verified",
        "os_reconstruction_status": "conditional",
        "continuum_limit_status": "missing",
        "physical_normalization_status": "missing",
        "bridge_source_status": "independent",
        "notes": "Local strong-coupling/lattice control kept separate from the open continuum step.",
    }
    values.update(overrides)
    return K41SplitCase(**values)


def build_matched_controls() -> list[K41SplitCase]:
    """Return one local positive control and eight mechanism-matched audits."""

    return [
        _base_case(),
        _base_case(
            case_id="open_os_continuum_transfer_hypothesis",
            role="scope_control",
            claim_mode="transfer_hypothesis",
            notes="The local theorem is retained while OS reconstruction and physical normalization stay open.",
        ),
        _base_case(
            case_id="negative_posthoc_envelope_same_values",
            role="matched_negative_control",
            envelope_predefined=False,
            notes="The same benign ratios fail when their envelope is chosen after inspection.",
        ),
        _base_case(
            case_id="negative_warm_corridor_escape",
            role="matched_negative_control",
            state_correlation_ratios=(0.994, 0.997, 1.82, 1.001, 1.002),
            warm_corridor_flags=(False, False, True, False, False),
            budget_source_id="derived_stress:single_warm_corridor",
            notes="One warm corridor cannot be averaged into a global mass-gap hint.",
        ),
        _base_case(
            case_id="negative_nonwarm_envelope_escape",
            role="matched_negative_control",
            state_correlation_ratios=(0.62, 0.997, 0.999, 1.001, 1.002),
            warm_corridor_flags=(False, False, False, False, False),
            budget_source_id="derived_stress:nonwarm_envelope_escape",
            notes="A state-ratio escape fails even without a warm-corridor label.",
        ),
        _base_case(
            case_id="negative_cross_scale_defect_budget_overrun",
            role="matched_negative_control",
            defect_budget_terms=(0.03, 0.04, 0.05, 0.06, 0.07),
            budget_source_id="derived_stress:defect_budget_overrun",
            notes="A positive local gap does not pay an unbounded cross-scale defect bill.",
        ),
        _base_case(
            case_id="negative_bridge_derived_from_local_gap",
            role="matched_negative_control",
            claim_mode="transfer_hypothesis",
            bridge_source_status="derived_from_local_gap_signal",
            notes="Reusing the successful local gap as its own continuum certificate is circular.",
        ),
        _base_case(
            case_id="negative_local_result_relabelled_continuum",
            role="matched_negative_control",
            claim_mode="claimed_continuum_transfer",
            notes="An open transfer hypothesis is asserted as a continuum theorem without the missing bridge.",
        ),
        _base_case(
            case_id="negative_assumed_local_gap",
            role="matched_negative_control",
            local_proof_status="assumed",
            notes="Even the positive-control side of the split must be evidenced rather than assumed.",
        ),
    ]


def build_results() -> list[K41SplitResult]:
    return [audit_case(case) for case in build_matched_controls()]


def _json_safe_row(result: K41SplitResult) -> dict[str, object]:
    return asdict(result)


def _markdown(results: Sequence[K41SplitResult]) -> str:
    lines = [
        f"# Yang--Mills K41-Style Transfer Split {DATE}",
        "",
        "Status: wording and defect-budget audit only; no Yang--Mills or turbulence claim.",
        "",
        "K41 is used only as a method example: its static, reference-normalized",
        "variational result is kept separate from dynamical cascade claims.  The",
        "Yang--Mills analogue is a local strong-coupling/lattice positive control",
        "kept separate from OS reconstruction, physical normalization and the",
        "continuum limit.",
        "",
        "| case | local positive | envelope | defect budget | warm corridors | split | decision |",
        "|---|---|---|---|---:|---|---|",
    ]
    for result in results:
        lines.append(
            f"| {result.case_id} | {result.positive_control_status} | "
            f"{result.state_correlation_envelope} | {result.defect_budget_status} | "
            f"{result.warm_corridor_count} | {result.split_status} | "
            f"{result.transfer_decision} |"
        )
    lines.extend(
        [
            "",
            "## Review waterline",
            "",
            "- Unconditional positive control: only the stated finite/local lattice result.",
            "- Transfer hypothesis: RP/OS reconstruction, physical normalization and the continuum limit are separate fields.",
            "- A predeclared state/correlation envelope and a finite cross-scale defect budget are necessary audit gates, not a proof.",
            "- Warm corridors, post-hoc envelopes and bridge certificates derived from the local gap fail closed.",
            "- K41 contributes wording discipline only; no turbulence dynamics are mapped onto gauge theory.",
            "- `claim_pass` is fixed to zero for all bundled controls.",
            "",
        ]
    )
    return "\n".join(lines)


def write_outputs(
    results: Sequence[K41SplitResult],
    output_dir: Path,
) -> dict[str, Path]:
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
        writer.writerows(_json_safe_row(result) for result in results)

    decision_counts: dict[str, int] = {}
    for result in results:
        decision_counts[result.transfer_decision] = decision_counts.get(result.transfer_decision, 0) + 1
    payload = {
        "schema_version": 1,
        "date": DATE,
        "scope_status": "k41_method_transfer_only_no_turbulence_or_yang_mills_claim",
        "claim_pass": 0,
        "k41_method_source": "https://doi.org/10.5281/zenodo.20131305",
        "rg_budget_source": "compute_os_capacity_ledger.py:strong_coupling_positive_control:k_max=10",
        "yang_mills_requirement_source": "https://www.claymath.org/library/monographs/MPPc.pdf",
        "transfer_decision_counts": decision_counts,
        "rows": [_json_safe_row(result) for result in results],
    }
    paths["json"].write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    paths["md"].write_text(_markdown(results), encoding="utf-8")
    return paths


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit a K41-style Yang--Mills transfer split")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "_results")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    results = build_results()
    paths = write_outputs(results, args.output_dir)
    for result in results:
        print(
            f"{result.case_id}: local={result.unconditional_positive_control} "
            f"envelope={result.state_correlation_envelope} "
            f"defects={result.defect_budget_total:.6g} decision={result.transfer_decision}"
        )
    print(f"claim_pass=0 json={paths['json']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
