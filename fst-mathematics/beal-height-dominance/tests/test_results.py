"""Offline checks of the curated Beal data; no proof or claim promotion."""
import csv
import hashlib
import json
import pathlib
import re
import unittest
from collections import Counter

ROOT = pathlib.Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "evidence"

def rows(name):
    with (EVIDENCE / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))

class CuratedResults(unittest.TestCase):
    def test_manifest_identities_and_csv_schemas(self):
        manifest = json.loads((ROOT / "manifest.json").read_text(encoding="utf-8"))
        actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob("*")
                  if p.is_file() and p.name != "manifest.json"
                  and "__pycache__" not in p.parts}
        self.assertEqual(actual, {item["path"] for item in manifest["files"]})
        for item in manifest["files"]:
            with self.subTest(path=item["path"]):
                path = ROOT / item["path"]
                self.assertEqual(path.stat().st_size, item["bytes"])
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), item["sha256"])
                if path.suffix == ".csv":
                    with path.open(encoding="utf-8", newline="") as stream:
                        reader = csv.DictReader(stream)
                        data = list(reader)
                    self.assertEqual(reader.fieldnames, item["columns"])
                    self.assertEqual(len(data), item["rows"])
                    self.assertTrue(all(None not in row for row in data))

    def test_near_miss_arithmetic_and_rank_identity(self):
        data = rows("near_power_b250_e3-7_z3-9_k5000.csv")
        self.assertEqual(len(data), 5000)
        self.assertEqual([int(row["rank"]) for row in data], list(range(1, 5001)))
        for row in data:
            A, x, B, y, C, z = (int(row[k]) for k in ("A", "x", "B", "y", "C", "z"))
            N, H, gap = (int(row[k]) for k in ("N", "H_Cz", "gap_signed"))
            self.assertEqual(N, A**x + B**y)
            self.assertEqual(H, C**z)
            self.assertEqual(gap, N-H)
            self.assertEqual(int(row["gap_abs"]), abs(gap))

    def test_variable_q_input_alignment_and_accounting(self):
        source = {row["rank"]: row for row in rows("near_power_b250_e3-7_z3-9_k5000.csv")}
        expected = [str(n) for n in range(1, 1001)]
        for Q in (251, 503, 1009, 2003):
            data = rows(f"variable_Q_sieve_Q{Q}_n1000.csv")
            self.assertEqual([row["rank"] for row in data], expected)
            for row in data:
                for key in ("A", "x", "B", "y", "C", "z"):
                    self.assertEqual(row[key], source[row["rank"]][key])
                self.assertEqual(int(row["passed"])+int(row["failed"]), int(row["nontrivial_tests"]))

    def test_root_ledger_control_classes_are_separate(self):
        names = [f"root_branch_drift_ledger_Q{Q}_n250.csv" for Q in (251, 503, 1009)]
        names += [f"ROOT_BRANCH_SURVIVOR_SCAN_2026-06-01_Q{Q}_n250.csv" for Q in (7, 13, 31, 61)]
        for name in names:
            data = rows(name)
            self.assertEqual(Counter(row["control_class"] for row in data),
                             {"actual_near_miss": 250, "random_root_sets": 250, "exact_power_anchor": 250})
            for kind in ("actual_near_miss", "random_root_sets", "exact_power_anchor"):
                self.assertEqual({int(row["rank"]) for row in data if row["control_class"] == kind},
                                 set(range(1, 251)))

    def test_rank8_prime_channel_matches_the_arithmetic_input(self):
        source = rows("near_power_b250_e3-7_z3-9_k5000.csv")[7]
        A, x, B, y, z = (int(source[k]) for k in ("A", "x", "B", "y", "z"))
        data = rows("ROOT_BRANCH_RANK8_MICRO_AUDIT_2026-06-02.csv")
        empty = []
        for row in data:
            q = int(row["q"])
            target = (pow(A, x, q)+pow(B, y, q)) % q
            roots = [c for c in range(q) if pow(c, z, q) == target]
            self.assertEqual(target, int(row["value_mod_q"]))
            self.assertEqual(roots, [int(c) for c in row["actual_roots"].split()])
            self.assertEqual(len(roots), int(row["actual_root_count"]))
            if not roots:
                empty.append(q)
        self.assertEqual(min(empty), 37)

    def test_natural_survivors_and_synthetic_controls_keep_their_scope(self):
        natural = rows("ABOVE_THRESHOLD_SURVIVOR_SEARCH_2026-06-02.csv")
        self.assertEqual(len(natural), 5)
        self.assertTrue(all(int(row["best_pass_Q"]) == 61 for row in natural))
        self.assertTrue(all(61 < int(row["first_empty_after_best"]) <= 97 for row in natural))
        self.assertTrue(all(int(row["low_sections_at_best"]) == 0 for row in natural))
        data = rows("Q97_PLUS_SYNTHETIC_CONTROL_2026-06-23.csv")
        self.assertEqual(Counter(row["control_type"] for row in data),
                         {"actual_natural_survivor": 5, "exact_registered_section": 5,
                          "shifted_registered_section": 5, "shuffled_branch_map": 5})
        self.assertTrue(all(row["claim_level"] == "no_claim_upgrade" for row in data))
        actual = [row for row in data if row["control_type"] == "actual_natural_survivor"]
        self.assertTrue(all(row["q97_local_pass"] == "False" for row in actual))

    def test_package_contains_no_internal_working_documents_or_local_paths(self):
        forbidden = {"BEWEISNOTIZ.md", "ZENODO_CREDENTIALS.md", "TODO.md", "PROJEKTPLAN.md",
                     "KONZEPT.md", "AGENTS.md", "CLAUDE.md"}
        for path in ROOT.rglob("*"):
            if not path.is_file() or "__pycache__" in path.parts:
                continue
            self.assertNotIn(path.name, forbidden)
            self.assertFalse(set(path.parts) & {"proof_notes", "_sources", "_results", "_reviews"})
            if path.suffix in {".tex", ".md", ".py", ".csv", ".json"}:
                content = path.read_text(encoding="utf-8")
                self.assertIsNone(re.search(r"[A-Za-z]:[\\/](?:Users|_Local_DEV)", content))

if __name__ == "__main__":
    unittest.main()
