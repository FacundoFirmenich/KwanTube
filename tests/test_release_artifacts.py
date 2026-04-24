import json
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class TestReleaseArtifacts(unittest.TestCase):
    def test_core_release_files_exist(self):
        required = [
            REPO_ROOT / "README.md",
            REPO_ROOT / "paper.md",
            REPO_ROOT / "paper.bib",
            REPO_ROOT / "CITATION.cff",
            REPO_ROOT / "reproduce_paper_results.py",
            REPO_ROOT / "generate_paper_figures.py",
            REPO_ROOT / "heom_acceptance_criteria.md",
        ]
        for path in required:
            self.assertTrue(path.exists(), f"Missing release artifact: {path.name}")

    def test_validation_report_has_required_sections(self):
        report_path = REPO_ROOT / "validation_report.json"
        self.assertTrue(report_path.exists(), "validation_report.json not found")

        data = json.loads(report_path.read_text(encoding="utf-8"))
        for key in ["_metadata", "_validation", "model_selection", "sbc", "lattice"]:
            self.assertIn(key, data)


if __name__ == "__main__":
    unittest.main()
