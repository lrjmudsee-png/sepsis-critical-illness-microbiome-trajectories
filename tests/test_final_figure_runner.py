"""Synthetic-input CLI safety checks; no patient data or R execution required."""
import csv
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


class FigureRunnerSafetyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        (self.root / "tools").mkdir()
        (self.root / "docs").mkdir()
        self.inputs = self.root / "inputs"
        self.inputs.mkdir()
        self.output = self.root / "new_output"
        self.runner = self.root / "tools" / "REPRODUCE_FINAL_FIGURES.py"
        shutil.copyfile(Path(__file__).resolve().parents[1] / "tools" / self.runner.name, self.runner)
        self.fixture = self.inputs / "synthetic.csv"
        self.fixture.write_text("group,value\nsynthetic,1\n", encoding="utf-8")
        with (self.root / "docs" / "FINAL_FIGURE_INPUT_MANIFEST.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["file", "bytes", "sha256"])
            writer.writeheader()
            writer.writerow({"file": self.fixture.name, "bytes": self.fixture.stat().st_size,
                             "sha256": hashlib.sha256(self.fixture.read_bytes()).hexdigest()})

    def tearDown(self):
        self.temp.cleanup()

    def invoke(self, *extra):
        return subprocess.run([sys.executable, str(self.runner), "--inputs", str(self.inputs),
                               "--output", str(self.output), *extra], capture_output=True, text=True)

    def test_check_only_writes_nothing(self):
        result = self.invoke("--check-only")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(self.output.exists())

    def test_missing_input_refused_before_build(self):
        self.fixture.unlink()
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Missing input", result.stderr)
        self.assertFalse(self.output.exists())

    def test_changed_input_refused_before_build(self):
        self.fixture.write_text("group,value\nsynthetic,2\n", encoding="utf-8")
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("checksum differs", result.stderr)
        self.assertFalse(self.output.exists())

    def test_existing_output_preserved(self):
        self.output.mkdir()
        sentinel = self.output / "preserve.txt"
        sentinel.write_text("unchanged", encoding="utf-8")
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Output already exists", result.stderr)
        self.assertEqual(sentinel.read_text(), "unchanged")

    def test_missing_rscript_refused(self):
        result = self.invoke("--rscript", "nonexistent_Rscript_for_safety_test")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Rscript not found", result.stderr)
        self.assertFalse(self.output.exists())


if __name__ == "__main__":
    unittest.main()
