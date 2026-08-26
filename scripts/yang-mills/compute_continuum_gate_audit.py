"""Audit the open Yang--Mills continuum-transfer gates without upgrading a claim.

The ledger turns four proof-note questions into machine-readable gates:

* the current T0--T3 status;
* the quantitative gap between the Fathi--Mikulincer--Shenfeld (FMS)
  existence bound and the contraction required by the proposed SU(N)
  typical-set RG route;
* the extra, independently summable GT2 defect budget that Kingman's theorem
  does not supply; and
* the precise point at which Gribov copies matter for a
  Dobrushin--Zegarlinski argument.

All bundled rows are analytical controls or project-status records.  They do
not establish a continuum Yang--Mills theory or a mass gap.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import sys
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Sequence


DATE = "2026-08-26"
STEM = f"YANG_MILLS_CONTINUUM_GATE_AUDIT_{DATE}"
ROOT = Path(__file__).resolve().parent

FMS_SOURCE = "https://arxiv.org/abs/2305.03786v3"
KINGMAN_SOURCE = "https://doi.org/10.1111/j.2517-6161.1968.tb00749.x"
ZEGA_SOURCE = "https://doi.org/10.1016/0022-1236(92)90073-R"
SINGER_SOURCE = "https://doi.org/10.1007/BF01609471"
CLAY_SOURCE = "https://www.claymath.org/wp-content/uploads/2022/02/MPPc.pdf"


@dataclass(frozen=True)
class ContinuumGateStatus:
    gate_id: str
    requirement: str
    current_status: str
    gate_pass: bool
    established_evidence: str
    missing_certificate: str
    failure_consequence: str


@dataclass(frozen=True)
class FMSTransportAudit:
    case_id: str
    beta: float
    full_log_lipschitz_scale: float
    typical_log_lipschitz_scale: float
    full_log_log_bound_scale: float
    typical_log_log_bound_scale: float
    target_contraction: float
    available_bound_class: str
    target_bound_class: str
    source_manifold_hypotheses_status: str
    typical_set_restriction_status: str
    su_n_conditional_instantiation_status: str
    rg_kernel_identification_status: str
    contraction_certificate_pass: bool
    decision: str
    claim_pass: bool


@dataclass(frozen=True)
class KingmanGT2Audit:
    case_id: str
    role: str
    mean_log_tau: float | None
    mean_contraction_status: str
    kingman_hypotheses_status: str
    gap_recursion_status: str
    defect_tail_model: str
    defect_budget_summable: bool
    defect_budget_upper_bound: float | None
    conditional_product_lower_bound: float
    physical_gt2_certificate: str
    decision: str
    claim_pass: bool


@dataclass(frozen=True)
class GribovDZAudit:
    case_id: str
    gauge_representation: str
    global_gauge_section_status: str
    product_local_specification_status: str
    uniform_single_site_lsi_status: str
    influence_matrix_status: str
    seam_or_horizon_capacity_status: str
    gribov_effect: str
    dobrushin_zegarlinski_status: str
    decision: str
    claim_pass: bool


def build_t0_t3_status() -> list[ContinuumGateStatus]:
    """Return the fail-closed proof-note status of the four continuum gates."""

    return [
        ContinuumGateStatus(
            gate_id="T0",
            requirement=(
                "OS-positive continuum Schwinger functions for a dense gauge-invariant class "
                "with uniform physical normalization"
            ),
            current_status="open_no_continuum_reconstruction_certificate",
            gate_pass=False,
            established_evidence="finite-lattice reflection positivity only",
            missing_certificate=(
                "continuum existence, OS axioms, nontriviality, and normalized observable convergence"
            ),
            failure_consequence="no physical continuum Hilbert space or Hamiltonian is reconstructed",
        ),
        ContinuumGateStatus(
            gate_id="T1",
            requirement=(
                "uniform SU(N) typical-set tail control and a quantitative transport or functional-inequality bound"
            ),
            current_status="partial_single_link_restricted_control_only",
            gate_pass=False,
            established_evidence=(
                "restricted finite/single-link typical-set estimates and an abstract FMS transport theorem"
            ),
            missing_certificate=(
                "FMS manifold-hypothesis instantiation for gauge conditionals, restriction/leakage control, "
                "and a contractive 1-O(beta^-1/2) RG map"
            ),
            failure_consequence="the weak-coupling Dobrushin/LSI input is not uniform across scales",
        ),
        ContinuumGateStatus(
            gate_id="T2",
            requirement=(
                "independently derived recursion kappa_(k+1) >= kappa_k(1-epsilon_k) with "
                "0 <= epsilon_k <= 1/2 and sum epsilon_k < infinity"
            ),
            current_status="open_gt2_defect_sequence_not_derived",
            gate_pass=False,
            established_evidence=(
                "conditional infinite-product lemma and finite negative-mean Kingman diagnostics"
            ),
            missing_certificate=(
                "stationary/ergodic RG cocycle, physical gap recursion, and independently summable defects"
            ),
            failure_consequence="average RG contraction cannot preserve a positive limiting PI/LSI constant",
        ),
        ContinuumGateStatus(
            gate_id="T3",
            requirement=(
                "reflection-compatible gauge-covariant blocking that preserves one-sided support and the OS cone "
                "at every scale"
            ),
            current_status="partial_finite_lattice_rp_only",
            gate_pass=False,
            established_evidence="finite Wilson/heat-kernel RP and a formal RP-compatible kernel criterion",
            missing_certificate=(
                "independent multiscale kernel, OS-crossing repair, and limit-stable cone preservation"
            ),
            failure_consequence="blocking may destroy reflection positivity before the continuum limit",
        ),
    ]


def audit_fms_transport(
    beta: float,
    *,
    target_coefficient: float = 0.5,
) -> FMSTransportAudit:
    """Expose the asymptotic existence-versus-contraction gap.

    The FMS manifold estimate has a double-exponential envelope in the
    log-Lipschitz scale L: ``log(log(K_bound)) = O(L^2)``.  The normalized
    proxies below record only that asymptotic class.  They are deliberately
    not presented as the theorem's numerical constants.
    """

    if not math.isfinite(beta) or beta <= 0.0:
        raise ValueError("beta must be finite and strictly positive")
    if not math.isfinite(target_coefficient) or target_coefficient <= 0.0:
        raise ValueError("target_coefficient must be finite and strictly positive")
    target = 1.0 - target_coefficient / math.sqrt(beta)
    if not 0.0 < target < 1.0:
        raise ValueError("target contraction must lie strictly between zero and one")

    full_scale = beta
    typical_scale = math.sqrt(beta)
    return FMSTransportAudit(
        case_id=f"fms_su_n_typical_beta_{beta:g}",
        beta=float(beta),
        full_log_lipschitz_scale=full_scale,
        typical_log_lipschitz_scale=typical_scale,
        full_log_log_bound_scale=full_scale**2,
        typical_log_log_bound_scale=typical_scale**2,
        target_contraction=target,
        available_bound_class="K <= exp(exp(O(L^2))); existence/nonexpansion bound only",
        target_bound_class="Lip(R_beta|T_beta) <= 1-c/sqrt(beta) < 1",
        source_manifold_hypotheses_status="abstract_theorem_not_yet_checked_for_rg_conditionals",
        typical_set_restriction_status="missing_boundary_and_bad_set_extension_certificate",
        su_n_conditional_instantiation_status="missing_uniform_in_boundary_volume_and_scale",
        rg_kernel_identification_status="missing_fms_map_not_identified_with_block_rg_kernel",
        contraction_certificate_pass=False,
        decision="open_fms_bound_does_not_imply_required_contraction",
        claim_pass=False,
    )


def build_fms_audits() -> list[FMSTransportAudit]:
    return [audit_fms_transport(beta) for beta in (4.0, 16.0, 64.0, 256.0)]


@lru_cache(maxsize=1)
def _existing_harmonic_control() -> object:
    """Reuse the existing matched Kingman false positive."""

    source_path = ROOT / "compute_k41_transfer_split_ledger.py"
    spec = importlib.util.spec_from_file_location("ym_k41_continuum_gate_source", source_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load Kingman control source: {source_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    matches = [
        result
        for result in module.build_results()
        if result.case_id == "negative_harmonic_worst_scale_corridor"
    ]
    if len(matches) != 1:
        raise RuntimeError("existing harmonic Kingman control is missing or duplicated")
    result = matches[0]
    if result.transfer_decision != "reject_non_summable_worst_scale_corridor":
        raise RuntimeError("existing harmonic Kingman control no longer fails closed")
    return result


def build_kingman_gt2_audits() -> list[KingmanGT2Audit]:
    """Separate Kingman drift from the independently summable GT2 input."""

    harmonic = _existing_harmonic_control()
    p2_budget = 0.1 * math.pi**2 / 6.0
    p2_product_lower_bound = math.exp(-2.0 * p2_budget)
    return [
        KingmanGT2Audit(
            case_id="existing_negative_mean_harmonic_corridor",
            role="matched_project_negative",
            mean_log_tau=float(harmonic.mean_log_rg_contraction),
            mean_contraction_status="finite_diagnostic_negative",
            kingman_hypotheses_status="not_verified_for_physical_rg_cocycle",
            gap_recursion_status="not_derived_from_rg_kernel",
            defect_tail_model="harmonic_non_summable",
            defect_budget_summable=False,
            defect_budget_upper_bound=None,
            conditional_product_lower_bound=0.0,
            physical_gt2_certificate="missing",
            decision="reject_negative_mean_without_summable_gt2_budget",
            claim_pass=False,
        ),
        KingmanGT2Audit(
            case_id="analytic_p2_gt2_control",
            role="conditional_positive_control",
            mean_log_tau=None,
            mean_contraction_status="not_needed_for_deterministic_product_lemma",
            kingman_hypotheses_status="not_invoked",
            gap_recursion_status="assumed_control_recursion",
            defect_tail_model="epsilon_k=0.1/k^2",
            defect_budget_summable=True,
            defect_budget_upper_bound=p2_budget,
            conditional_product_lower_bound=p2_product_lower_bound,
            physical_gt2_certificate="missing_control_only",
            decision="conditional_product_lemma_pass_no_physical_rg_input",
            claim_pass=False,
        ),
        KingmanGT2Audit(
            case_id="analytic_harmonic_gt2_negative",
            role="matched_analytic_negative",
            mean_log_tau=-0.1,
            mean_contraction_status="negative_by_construction",
            kingman_hypotheses_status="synthetic_stationary_diagnostic_only",
            gap_recursion_status="assumed_control_recursion",
            defect_tail_model="epsilon_k=0.1/k",
            defect_budget_summable=False,
            defect_budget_upper_bound=None,
            conditional_product_lower_bound=0.0,
            physical_gt2_certificate="missing_control_only",
            decision="reject_harmonic_defect_even_with_negative_mean",
            claim_pass=False,
        ),
    ]


def build_gribov_dz_audits() -> list[GribovDZAudit]:
    """Clarify when the Gribov ambiguity enters a DZ proof.

    Singer's obstruction concerns a global continuous gauge section.  A
    gauge-unfixed lattice specification on the compact product SU(N)^E does
    not choose such a section, so Gribov copies are not a direct input to its
    conditional laws.  This does not prove the DZ bounds.  Gauge-fixed or
    quotient implementations need additional seam/capacity and uniformity
    certificates and fail closed here.
    """

    return [
        GribovDZAudit(
            case_id="unfixed_compact_link_product",
            gauge_representation="gauge-unfixed SU(N)^E link variables",
            global_gauge_section_status="not_used",
            product_local_specification_status="available_in_principle_finite_lattice",
            uniform_single_site_lsi_status="open_quantitative_bound",
            influence_matrix_status="open_outside_strong_coupling",
            seam_or_horizon_capacity_status="not_applicable_without_gauge_slice",
            gribov_effect="no_direct_obstruction_because_no_global_gauge_choice_is_made",
            dobrushin_zegarlinski_status="not_blocked_by_gribov_but_analytic_bounds_open",
            decision="diagnostic_only_keep_unfixed_route",
            claim_pass=False,
        ),
        GribovDZAudit(
            case_id="global_gauge_fixed_slice",
            gauge_representation="single global nonabelian gauge slice",
            global_gauge_section_status="obstructed_by_singer_gribov_ambiguity",
            product_local_specification_status="not_certified_global_or_local_after_fixing",
            uniform_single_site_lsi_status="missing",
            influence_matrix_status="missing_nonlocal_gauge_fixing_can_change_conditionals",
            seam_or_horizon_capacity_status="missing",
            gribov_effect="global_section_failure_and_possible_nonlocal_or_singular_conditionals",
            dobrushin_zegarlinski_status="blocked_no_uniform_local_specification_certificate",
            decision="reject_global_gauge_slice_shortcut",
            claim_pass=False,
        ),
        GribovDZAudit(
            case_id="patchwise_orbit_space",
            gauge_representation="orbit-space atlas or fundamental modular patches",
            global_gauge_section_status="patchwise_only",
            product_local_specification_status="requires_overlap_and_measurable_disintegration_proof",
            uniform_single_site_lsi_status="missing_uniform_patch_bound",
            influence_matrix_status="missing_uniform_cross_patch_bound",
            seam_or_horizon_capacity_status="missing_summable_capacity_certificate",
            gribov_effect="patch seams may concentrate transport_or_gradient_defects",
            dobrushin_zegarlinski_status="blocked_until_patch_and_capacity_controls_close",
            decision="diagnostic_only_patchwise_route_open",
            claim_pass=False,
        ),
        GribovDZAudit(
            case_id="measure_zero_horizon_shortcut",
            gauge_representation="gauge-fixed typical set with horizon deleted",
            global_gauge_section_status="still_not_global",
            product_local_specification_status="not_restored_by_measure_zero_statement",
            uniform_single_site_lsi_status="missing_capacity_or_hardy_control",
            influence_matrix_status="missing",
            seam_or_horizon_capacity_status="measure_zero_only_not_a_capacity_bound",
            gribov_effect="null_mass_does_not_control_gradients_conditionals_or_transport_cost",
            dobrushin_zegarlinski_status="blocked_measure_zero_is_insufficient",
            decision="reject_measure_zero_gribov_shortcut",
            claim_pass=False,
        ),
    ]


def _summary_rows(payload: dict[str, object]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for family, key, status_key in (
        ("t0_t3", "t0_t3", "current_status"),
        ("fms_transport", "fms_transport", "decision"),
        ("kingman_gt2", "kingman_gt2", "decision"),
        ("gribov_dz", "gribov_dz", "decision"),
    ):
        for item in payload[key]:
            rows.append(
                {
                    "audit_family": family,
                    "case_id": item.get("gate_id", item.get("case_id")),
                    "status": item[status_key],
                    "claim_pass": item.get("claim_pass", item.get("gate_pass", False)),
                }
            )
    return rows


def build_payload() -> dict[str, object]:
    t0_t3 = [asdict(item) for item in build_t0_t3_status()]
    fms = [asdict(item) for item in build_fms_audits()]
    kingman = [asdict(item) for item in build_kingman_gt2_audits()]
    gribov = [asdict(item) for item in build_gribov_dz_audits()]
    return {
        "schema_version": 1,
        "date": DATE,
        "scope_status": "proof_status_audit_only_no_continuum_yang_mills_claim",
        "claim_pass": 0,
        "proof_note_writeback_status": "blocked_onedrive_cloud_lock_and_host_artifacts",
        "sources": {
            "fms": FMS_SOURCE,
            "kingman": KINGMAN_SOURCE,
            "dobrushin_zegarlinski": ZEGA_SOURCE,
            "gribov_singer": SINGER_SOURCE,
            "clay_problem": CLAY_SOURCE,
        },
        "t0_t3": t0_t3,
        "fms_transport": fms,
        "kingman_gt2": kingman,
        "gribov_dz": gribov,
    }


def _markdown(payload: dict[str, object]) -> str:
    lines = [
        f"# Yang--Mills Continuum-Gate Audit {DATE}",
        "",
        "Status: proof-note audit only; no continuum Yang--Mills or mass-gap claim.",
        "",
        "The OneDrive proof-note writeback is blocked by the current cloud-lock/host-artifact gate.",
        "This generated review is the deterministic source for a later safe writeback.",
        "",
        "## T0--T3 status",
        "",
        "| gate | status | established | missing |",
        "|---|---|---|---|",
    ]
    for item in payload["t0_t3"]:
        lines.append(
            f"| {item['gate_id']} | {item['current_status']} | {item['established_evidence']} | "
            f"{item['missing_certificate']} |"
        )
    lines.extend(
        [
            "",
            "## FMS quantitative gap",
            "",
            "FMS proves existence of globally Lipschitz transports under its stated geometric",
            "hypotheses. Its available general-manifold envelope is non-contractive and has",
            "`log log K = O(L^2)`. Substituting the project scalings `L=O(beta)` and",
            "`L_typ=O(sqrt(beta))` changes the envelope from `exp(exp(O(beta^2)))` to",
            "`exp(exp(O(beta)))`; it does not yield the required",
            "`1-c/sqrt(beta) < 1`. Restriction to a typical set, uniform SU(N) conditional",
            "instantiation, and identification with the RG kernel remain separate certificates.",
            "",
            "| beta | full log-log scale | typical log-log scale | target | decision |",
            "|---:|---:|---:|---:|---|",
        ]
    )
    for item in payload["fms_transport"]:
        lines.append(
            f"| {item['beta']:g} | {item['full_log_log_bound_scale']:g} | "
            f"{item['typical_log_log_bound_scale']:g} | {item['target_contraction']:.6f} | "
            f"{item['decision']} |"
        )
    lines.extend(
        [
            "",
            "## Kingman versus GT2",
            "",
            "Kingman's theorem can identify an almost-sure asymptotic rate for a stationary",
            "ergodic subadditive process. A negative rate does not derive the physical recursion",
            "`kappa_(k+1) >= kappa_k(1-epsilon_k)` and does not imply",
            "`sum epsilon_k < infinity`. The p-series control shows only the conditional",
            "infinite-product lemma; the harmonic controls show why mean contraction is not enough.",
            "",
            "| case | mean status | defect tail | summable | product lower bound | decision |",
            "|---|---|---|---|---:|---|",
        ]
    )
    for item in payload["kingman_gt2"]:
        lines.append(
            f"| {item['case_id']} | {item['mean_contraction_status']} | "
            f"{item['defect_tail_model']} | {item['defect_budget_summable']} | "
            f"{item['conditional_product_lower_bound']:.6f} | {item['decision']} |"
        )
    lines.extend(
        [
            "",
            "## Gribov copies and Dobrushin--Zegarlinski",
            "",
            "Singer's obstruction targets a global continuous non-Abelian gauge section. It does",
            "not directly enter a gauge-unfixed finite-lattice specification on `SU(N)^E`, because",
            "that route chooses no global slice. The single-site LSI and influence bounds are still",
            "open outside the verified regime. A global gauge-fixed or patchwise quotient route",
            "must additionally control nonlocal/singular conditionals and seam or horizon capacity;",
            "deleting a measure-zero horizon is not such a control.",
            "",
            "| case | Gribov effect | DZ status | decision |",
            "|---|---|---|---|",
        ]
    )
    for item in payload["gribov_dz"]:
        lines.append(
            f"| {item['case_id']} | {item['gribov_effect']} | "
            f"{item['dobrushin_zegarlinski_status']} | {item['decision']} |"
        )
    lines.extend(
        [
            "",
            "## Claim waterline",
            "",
            "- All T0--T3 gates remain open or partial; none passes.",
            "- FMS supplies an existence theorem, not the requested RG contraction certificate.",
            "- Kingman mean contraction and GT2 summability are independent proof obligations.",
            "- Gribov ambiguity is avoided, not solved, by the gauge-unfixed product route.",
            "- `claim_pass` is fixed to zero for every bundled record.",
            "",
            "## Primary sources",
            "",
            f"- FMS transport: {FMS_SOURCE}",
            f"- Kingman subadditive theorem: {KINGMAN_SOURCE}",
            f"- Dobrushin--Zegarlinski/LSI: {ZEGA_SOURCE}",
            f"- Singer Gribov obstruction: {SINGER_SOURCE}",
            f"- Clay Yang--Mills problem statement: {CLAY_SOURCE}",
            "",
        ]
    )
    return "\n".join(lines)


def write_outputs(payload: dict[str, object], output_dir: Path) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "csv": output_dir / f"{STEM}.csv",
        "json": output_dir / f"{STEM}.json",
        "md": output_dir / f"{STEM}.md",
    }
    paths["json"].write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    paths["md"].write_text(_markdown(payload), encoding="utf-8")
    with paths["csv"].open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("audit_family", "case_id", "status", "claim_pass"),
        )
        writer.writeheader()
        writer.writerows(_summary_rows(payload))
    return paths


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit the Yang--Mills continuum transfer gates")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "_results")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    payload = build_payload()
    paths = write_outputs(payload, args.output_dir)
    statuses = ", ".join(
        f"{item['gate_id']}={item['current_status']}" for item in payload["t0_t3"]
    )
    print(statuses)
    print(
        f"fms_rows={len(payload['fms_transport'])} "
        f"kingman_rows={len(payload['kingman_gt2'])} "
        f"gribov_rows={len(payload['gribov_dz'])} claim_pass=0 json={paths['json']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
