# Final submission figure reproduction

The authoritative final display layer is `analysis/final_figures/`, derived from the corrected 2 October 2026 submission figure build. Use `tools/04_reproduce_final_figures.py`, not the historical cosmetic scripts, to reproduce the five main figures and six supplementary figures. The new interface requires explicit input and output directories and does not depend on a local project root or an E: drive.

This route renders frozen results; it does not fit models, recompute confidence intervals, run statistical tests, change patient mappings, or process FASTQ reads. Selected-statistics reproduction remains a separate route documented in `SELECTED_STATISTICS_REPRODUCTION.md` and uses Supplementary Data 5.

## Inputs and access

Provide the 24 original CSV files listed in `FINAL_FIGURE_INPUT_MANIFEST.csv` in one external folder. SHA256 verification is mandatory. This manifest contains filenames, byte counts and checksums only, not patient records. Some figure inputs contain repeated patient observations or sample-level distances. They are not published in this software repository, and the MIT licence does not apply to third-party data. The corresponding frozen figure-input bundle is held by the authors for reviewer access pending redistribution review. Public cohort accessions alone do not supply these derived inputs. A repository clone alone is therefore not a self-contained reproduction package.

## Environment

Tested with R 4.4.0 on Windows and readr 2.2.0, dplyr 1.2.1, tidyr 1.3.2, tibble 3.3.1, ragg 1.5.2, ggplot2 4.0.3 and patchwork 1.3.2. The scripts do not install packages automatically. Python 3.12 with `requirements-figures.txt` is used for optional supplementary PDF assembly. Arial must be installed for matching figure typography; licensed Arial regular and bold font files must be supplied separately for PDF assembly. No fonts are redistributed. Other operating systems and font environments have not been validated for identical rendering.

## Run

From any working directory, use absolute paths or paths relative to that directory:

```text
python /path/to/repository/tools/04_reproduce_final_figures.py --inputs /path/to/frozen_csvs --output /path/to/new_figure_build --rscript /path/to/Rscript --font-dir /path/to/Arial_fonts
```

Omit `--font-dir` for the 33 figure exports and generated legend files without supplementary PDF assembly. Use `--check-only` to verify inputs without writing outputs. Existing output directories are refused. Input checksum failures are reported before a build begins.

The build writes PDF, 300-dpi PNG and LZW-compressed 600-dpi TIFF for each of 11 figures, export metadata, R session information, generated full-length audit legends and a checksum report. The full-length main-figure audit legends are not the condensed journal-facing manuscript legends; scientific values and plot definitions are shared. Supplementary PDF assembly uses the generated supplementary legends.

## Scope and historical provenance

The 4 October naming migration renames programs and software-location references only; scientific calculations and frozen data/output names are preserved. Previously missing historical 96F2-F11 cosmetic scripts are not recovered by this update; the corrected final display layer provides the current figure-rendering route. This update does not establish complete raw-to-manuscript reproduction, clinical roster adjudication, or permission to redistribute input data. No GitHub release, version tag or published Zenodo DOI is created by this update.
