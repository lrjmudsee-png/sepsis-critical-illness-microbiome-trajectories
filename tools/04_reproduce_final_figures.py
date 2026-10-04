"""Rebuild the submission figures from externally supplied frozen CSVs.

This rendering route is not raw-read preprocessing or statistical inference.
Only code, documentation and input checksums are distributed in the repository.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--font-dir", type=Path, help="Arial font directory; enables supplementary PDF assembly")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    inputs = args.inputs.resolve(strict=True)
    output = args.output.resolve()
    manifest = repo / "docs" / "FINAL_FIGURE_INPUT_MANIFEST.csv"
    failures = []
    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        records = list(csv.DictReader(handle))
    for record in records:
        path = inputs / record["file"]
        if not path.is_file():
            failures.append(f"Missing input: {record['file']}")
        elif hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
            failures.append(f"Input checksum differs: {record['file']}")
    if failures:
        parser.error("\n".join(failures))
    print(f"Verified {len(records)} frozen input CSVs.")
    if args.check_only:
        return
    if output.exists():
        parser.error("Output already exists; choose a new directory. Nothing will be overwritten.")
    rscript = shutil.which(args.rscript)
    if not rscript:
        parser.error("Rscript not found; install the documented environment or provide --rscript.")
    if args.font_dir:
        for name in ("arial.ttf", "arialbd.ttf"):
            if not (args.font_dir / name).is_file():
                parser.error(f"Missing font: {name}. Fonts are not redistributed with this repository.")
        try:
            import PIL  # noqa: F401
            import reportlab  # noqa: F401
        except ImportError as exc:
            parser.error(f"Install requirements-figures.txt before PDF assembly: {exc}")
    scripts = repo / "analysis" / "final_figures"
    commands = [
        [rscript, str(scripts / "10_001_rebuild_submission_figures.R"), str(inputs), str(output)],
        [rscript, str(scripts / "10_002_generate_figure_legends.R"), str(inputs), str(output)],
    ]
    if args.font_dir:
        commands.append([sys.executable, str(scripts / "10_003_assemble_supplementary_pdf.py"),
                         str(output), str(output / "Supplementary_Information.pdf"),
                         "--font-dir", str(args.font_dir.resolve())])
    # The first R stage owns creation of the output folder and never overwrites it.
    for index, command in enumerate(commands, 1):
        result = subprocess.run(command, text=True, encoding="utf-8", errors="replace", capture_output=True)
        if output.is_dir():
            (output / f"reproduction_stage_{index}.log").write_text(result.stdout + result.stderr, encoding="utf-8")
        if result.returncode:
            raise RuntimeError(f"Reproduction stage {index} failed: {result.stderr or result.stdout}")
    rendered = list((output / "02_MAIN_FIGURES").glob("*")) + list((output / "03_SUPPLEMENTARY_FIGURES").glob("*"))
    if len(rendered) != 33:
        raise RuntimeError(f"Expected 33 figure exports, got {len(rendered)}.")
    report = {"input_files_verified": len(records), "figure_exports": len(rendered),
              "inference_rerun": False, "raw_read_processing_rerun": False,
              "supplementary_pdf_assembled": bool(args.font_dir),
              "output_sha256": {p.relative_to(output).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sorted(rendered)}}
    (output / "REPRODUCTION_REPORT.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Complete: {len(rendered)} figure exports. See {output / 'REPRODUCTION_REPORT.json'}")


if __name__ == "__main__":
    main()
