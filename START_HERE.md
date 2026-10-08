# Start here

This repository distributes software and aggregate audit evidence for **Longitudinal gut-microbiome trajectories in sepsis and critical illness**. It does not redistribute raw reads or reviewer-only patient inputs.

## Current reproduction entry points

| Task | Entry point | Required inputs |
| --- | --- | --- |
| Selected statistics | `tools/01_reproduce_selected_statistics.R` | External flat-file statistical inputs described in `docs/SELECTED_STATISTICS_REPRODUCTION.md` |
| Specimen sensitivities | `tools/02_reproduce_specimen_sensitivity.R` | External specimen and pair inputs, plus the documented metafor environment |
| Aggregate forest example | `tools/03_render_direction_forest.R` | The public aggregate forest source table |
| Final figures | `tools/04_reproduce_final_figures.py` | The 24 external frozen CSV files in `docs/FINAL_FIGURE_INPUT_MANIFEST.csv` |
| Prepared terminal chain | `tools/05_run_prepared_terminal_chain.ps1` | A separately prepared run root and its ordered 11-script chain; do not supply the entire historical analysis folder |

Example figure rendering:

```text
python tools/04_reproduce_final_figures.py --inputs /path/to/frozen_csvs --output /path/to/new_output --rscript /path/to/Rscript --font-dir /path/to/Arial_fonts
```

## Naming convention and catalogue

Analysis programs follow `PP_NNN_descriptive_function.ext`: `PP` is a workflow phase, and `NNN` is an inventory order within that phase. A PowerShell launcher normally shares its program's prefix and adds `_run`. Portable tools have a two-digit tool number. Superseded versions are retained with descriptive correction/revision names rather than letter-coded stage names.

**These numbers are a catalogue, not a command to execute every script.** The historical source collection includes alternatives, diagnostic trials, superseded versions and paths to frozen external inputs. Use the documented current entry point for your task. Running every numbered historical script sequentially is neither required nor validated.

| Phase | Scope |
| --- | --- |
| 01 | Metadata collection, patient/sample/time reconciliation and freezes |
| 02 | Read processing, taxonomy assignment and sequence QC |
| 03 | Longitudinal ecology and paired contrasts |
| 04 | Clinical metadata and context |
| 05 | Taxonomic trajectories and exploratory index development |
| 06 | Infection-source analysis |
| 07 | External validation, reconciliation and governance |
| 08 | Secondary direction-analysis workflow |
| 09 | Historical manuscript and table assembly |
| 10 | Current final-figure display layer |
| 11 | Workflow utilities |

`docs/CODE_FILENAME_MAP.csv` maps all 364 old script paths to their renamed paths and records original SHA256 values. `docs/RUN_ORDER.csv` is the current complete script catalogue. `docs/LEGACY_RUN_ORDER.csv` preserves the prior historical workflow overview.

## Preservation and limits

This migration changes script names and references to software locations. Statistical formulas, input/output filenames, frozen cohort selection, historical results and audit evidence are not renamed. Old stage labels inside result directories, scientific report fields and historical evidence still identify the original analyses; they are not stale executable filenames.

Historical scripts remain partly dependent on the original Windows project, package environment and frozen inputs. Missing historical source/assets are not reconstructed by renaming. Historical evidence records and `MANIFEST_SHA256.csv` intentionally preserve their original names and hashes. Use the current `CODE_FILE_MANIFEST.csv` for this distribution. The new `evidence/NAMING_MIGRATION_20261004/` reports tests of the renamed software, not a new full FASTQ/DADA2 rerun.

Software authorship remains Lu, Rongji and the original software license remains MIT. Software version `1.0.0`, from commit `1e998aaa9b19e48643996c95d78b7cf7b74b1807`, is published at [Zenodo DOI 10.5281/zenodo.23134779](https://doi.org/10.5281/zenodo.23134779). The archive is immutable; later documentation updates on the default branch do not change this archived software snapshot. Do not create a duplicate archive for the same snapshot. Additional reporting-recovery helpers are supplied with the manuscript's Supplementary Data 5 rather than the archived software.

