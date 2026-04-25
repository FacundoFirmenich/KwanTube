import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "src" / "scripts" / "bayesian_heom_hierarchy_v2.py"
INPUT_CSV = REPO_ROOT / "src" / "heom_bayes_input_current.csv"


class TestBayesianHeomHierarchyV2Smoke(unittest.TestCase):
    def test_v2_runs_and_produces_expected_core_outputs(self):
        self.assertTrue(SCRIPT.exists())
        self.assertTrue(INPUT_CSV.exists())

        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp) / "heom_bayes_out_v2_test"
            cmd = [
                sys.executable,
                str(SCRIPT),
                str(INPUT_CSV),
                "--output-dir",
                str(out_dir),
                "--draws",
                "4000",
                "--seed",
                "20260424",
            ]
            subprocess.run(cmd, check=True, cwd=REPO_ROOT)

            required = [
                "diagnostics_v2.txt",
                "group_loglinear_summary.csv",
                "hierarchy_global_contraction.csv",
                "hierarchy_group_shrinkage.csv",
                "extrapolated_jumps.csv",
                "level_reference_checks.csv",
            ]
            for name in required:
                self.assertTrue((out_dir / name).exists(), f"Missing v2 output: {name}")

            # Minimal numerical sanity on global contraction ratio.
            with (out_dir / "hierarchy_global_contraction.csv").open(newline="", encoding="utf-8") as f:
                rows = list(csv.DictReader(f))
            global_r = next(r for r in rows if r["parameter"] == "global_r")
            mean_r = float(global_r["mean"])
            self.assertGreater(mean_r, 0.0)
            self.assertLess(mean_r, 1.0)


if __name__ == "__main__":
    unittest.main()
