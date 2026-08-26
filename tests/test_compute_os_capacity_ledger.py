import csv
import importlib.util
import json
import sys
from collections import defaultdict
from pathlib import Path

import pytest


SCRIPT_PATH = (
    Path(__file__).parents[1] / "scripts" / "yang-mills" / "compute_os_capacity_ledger.py"
)
SPEC = importlib.util.spec_from_file_location("compute_os_capacity_ledger_under_test", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_generated_window_ledger_is_predefined_and_normalized():
    rows = MODULE.generated_v2_control_rows(k_max=20)
    payload = MODULE.v2_payload(rows, source_path=None, date_tag="test")

    expected_fields = {
        "rg_window_id",
        "window_predefined",
        "scale_occupancy",
        "safe_signal_share",
        "os_capacity_share",
        "defect_share_over_occupancy",
        "bad_run_switch_share",
        "alternate_blocking_control_ratio",
        "nonlocal_tail_cost",
        "transfer_status",
    }
    by_scenario = defaultdict(list)
    for window in payload["window_ledger"]:
        assert expected_fields <= window.keys()
        assert window["window_predefined"] is True
        assert 0.0 < window["scale_occupancy"] <= 1.0
        assert 0.0 <= window["bad_run_switch_share"] <= 1.0
        by_scenario[window["scenario"]].append(window)

    assert len(by_scenario) == 5
    for windows in by_scenario.values():
        assert sum(window["scale_occupancy"] for window in windows) == pytest.approx(1.0)
        assert sum(window["safe_signal_share"] for window in windows) == pytest.approx(1.0)
        assert sum(window["os_capacity_share"] for window in windows) == pytest.approx(1.0)

    early_summable = next(
        window
        for window in by_scenario["summable_os_capacity_control"]
        if window["rg_window_id"] == "rg_early_10pct"
    )
    assert early_summable["os_capacity_share"] > early_summable["scale_occupancy"]
    assert early_summable["defect_share_over_occupancy"] > 1.0
    assert early_summable["transfer_status"] == "flagged_bad_channel_control_only"


def test_post_hoc_window_is_blocked_even_when_base_control_passes():
    rows = [
        dict(row)
        for row in MODULE.generated_v2_control_rows(k_max=10)
        if row["scenario"] == "summable_os_capacity_control"
    ]
    rows[0]["window_predefined"] = False

    payload = MODULE.v2_payload(rows, source_path=None, date_tag="test")
    early = next(
        window
        for window in payload["window_ledger"]
        if window["rg_window_id"] == "rg_early_10pct"
    )

    assert early["window_predefined"] is False
    assert early["transfer_status"] == "blocked_post_hoc_window"


def test_noncontiguous_window_ids_are_rejected():
    rows = [
        dict(row)
        for row in MODULE.generated_v2_control_rows(k_max=6)
        if row["scenario"] == "summable_os_capacity_control"
    ]
    rows[0]["rg_window_id"] = "scattered"
    rows[2]["rg_window_id"] = "scattered"

    with pytest.raises(ValueError, match="not contiguous"):
        MODULE.v2_payload(rows, source_path=None, date_tag="test")


def test_isolated_cli_run_writes_window_schema(tmp_path):
    data_dir = tmp_path / "data"
    output_dir = tmp_path / "results"

    exit_code = MODULE.main(
        [
            "--rows",
            "10",
            "--date-tag",
            "test",
            "--data-dir",
            str(data_dir),
            "--output-dir",
            str(output_dir),
        ]
    )

    assert exit_code == 0
    input_path = data_dir / "OS_CAPACITY_LEDGER_V2_CONTROL_INPUT_test.csv"
    json_path = output_dir / "OS_CAPACITY_LEDGER_V2_test.json"
    window_path = output_dir / "OS_CAPACITY_LEDGER_V2_WINDOWS_test.csv"
    assert input_path.is_file()
    assert json_path.is_file()
    assert window_path.is_file()

    roundtrip_rows = MODULE.read_input_csv(input_path)
    assert len(roundtrip_rows) == 50
    assert all(row["window_predefined"] is True for row in roundtrip_rows)
    assert all(row["alternate_blocking_control"] >= 0.0 for row in roundtrip_rows)

    payload = json.loads(json_path.read_text(encoding="utf-8"))
    assert payload["status"] == "ledger/control run only; no Yang-Mills mass-gap claim"
    assert payload["window_ledger"]

    with window_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == MODULE.WINDOW_LEDGER_FIELDS
        assert len(list(reader)) == len(payload["window_ledger"])


@pytest.mark.parametrize("value", ["yes", "1", "", "pre_registered"])
def test_window_predefined_parser_is_strict(value):
    with pytest.raises(ValueError, match="Expected boolean"):
        MODULE.parse_bool(value)
