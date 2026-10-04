"""Software filename/link inventory; no clinical inputs or processing required."""
import csv
from pathlib import Path
import re
import unittest


class NamingInventoryTests(unittest.TestCase):
    def setUp(self):
        self.repo = Path(__file__).resolve().parents[1]
        with (self.repo / 'docs/CODE_FILENAME_MAP.csv').open(encoding='utf-8-sig', newline='') as handle:
            self.rows = list(csv.DictReader(handle))

    def test_all_original_programs_have_unique_existing_targets(self):
        self.assertEqual(len(self.rows), 364)
        self.assertEqual(len({r['old_path'] for r in self.rows}), len(self.rows))
        self.assertEqual(len({r['new_path'].lower() for r in self.rows}), len(self.rows))
        for row in self.rows:
            self.assertTrue((self.repo / row['new_path']).is_file(), row['new_path'])
            self.assertFalse((self.repo / row['old_path']).exists(), row['old_path'])

    def test_every_distributed_program_uses_the_numbered_convention(self):
        for path in self.repo.rglob('*'):
            if path.is_file() and path.suffix.lower() in ('.r', '.py', '.ps1'):
                self.assertRegex(path.name, r'^(?:test_)?\d{2}_')

    def test_current_figure_runner_dependencies_exist(self):
        code = (self.repo / 'tools/04_reproduce_final_figures.py').read_text(encoding='utf-8')
        for name in ('10_001_rebuild_submission_figures.R', '10_002_generate_figure_legends.R',
                     '10_003_assemble_supplementary_pdf.py'):
            self.assertIn(name, code)
            self.assertTrue((self.repo / 'analysis/final_figures' / name).is_file())
        renderer = (self.repo / 'analysis/final_figures/10_001_rebuild_submission_figures.R').read_text(encoding='utf-8')
        self.assertIn('10_004_figure_display_settings.R', renderer)

    def test_current_code_does_not_reference_an_original_script_basename(self):
        # Historical evidence and the explicit old-to-new map intentionally retain old names.
        old_names = {Path(r['old_path']).name for r in self.rows}
        old_regex = re.compile('|'.join(re.escape(n) for n in sorted(old_names, key=len, reverse=True)))
        for row in self.rows:
            code = (self.repo / row['new_path']).read_text(encoding='utf-8-sig')
            self.assertIsNone(old_regex.search(code), row['new_path'])


if __name__ == '__main__':
    unittest.main()
