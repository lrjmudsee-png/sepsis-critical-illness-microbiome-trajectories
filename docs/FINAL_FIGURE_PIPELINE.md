# Final submission figure reproduction

The authoritative final display layer is `analysis/final_figures/`, derived from the corrected 2 October 2026 submission figure build. Use `tools/04_reproduce_final_figures.py`, not the historical cosmetic scripts, to reproduce the five main figures and six supplementary figures. The new interface requires explicit input and output directories and does not depend on a local project root or an E: drive.

This route renders frozen results; it does not fit models, recompute confidence intervals, run statistical tests, change patient mappings, or process FASTQ reads. Selected-statistics reproduction remains a separate route documented in `SELECTED_STATISTICS_REPRODUCTION.md` and uses Supplementary Data 5.

## Inputs and access

Provide the 24 original CSV files listed in `FINAL_FIGURE_INPUT_MANIFEST.csv` in one external folder. SHA256 verification is mandatory. This manifest contains filenames, byte counts and checksums only, not patient records. Some figure inputs contain repeated public-study patient observations or sample-level distances. They are not published in this software repository, and the MIT licence does not apply to third-party data. The manuscript's Data Set S6 supplies the 24 derived frozen CSVs, figure code and checksums for editorial and peer-review reproduction. Public cohort accessions alone do not supply these derived inputs. A repository clone alone therefore needs the separately supplied figure-input bundle.

## Environment

Tested with R 4.4.0 on Windows and readr 2.2.0, dplyr 1.2.1, tidyr 1.3.2, tibble 3.3.1, ragg 1.5.2, ggplot2 4.0.3 patchwork 1.3.2, svglite and jsonlite. Exact current versions are recorded in the generated R session information. The scripts do not install packages automatically. Python 3.12.14 with `requirements-figures.txt` is used for optional supplementary PDF assembly. Arial must be installed for matching figure typography; licensed Arial regular and bold font files must be supplied separately for PDF assembly. No fonts are redistributed. Other operating systems and font environments have not been validated for identical rendering.

## Run

From any working directory, use absolute paths or paths relative to that directory:

```text
python /path/to/repository/tools/04_reproduce_final_figures.py --inputs /path/to/frozen_csvs --output /path/to/new_figure_build --rscript /path/to/Rscript --font-dir /path/to/Arial_fonts
```

Omit `--font-dir` for the 44 figure exports and generated legend files without supplementary PDF assembly. Use `--check-only` to verify inputs without writing outputs. Existing output directories are refused. Input checksum failures are reported before a build begins.

The build writes PDF, editable SVG, 300-dpi PNG and LZW-compressed 600-dpi TIFF for each of 11 figures, export metadata, R session information, synchronized figure legends and a checksum report. The final main-figure legends are synchronized with the mSystems manuscript, including continuation paragraphs in Fig. 3. Supplemental legends describe the same frozen data. The scientific figure plates retain the 26 October dimensions and values; submission PDFs may add surrounding review legends and page furniture without changing the plates. Supplementary PDF assembly uses the generated supplementary legends.

## Scope and historical provenance

The 4 October naming migration renames programs and software-location references only; scientific calculations and frozen data/output names are preserved. Previously missing historical 96F2-F11 cosmetic scripts are not recovered by this update; the corrected final display layer provides the current figure-rendering route. This update does not establish complete raw-to-manuscript reproduction, clinical roster adjudication, or permission to redistribute input data. No GitHub release, version tag or published Zenodo DOI is created by this update.

## Display revision for mSystems preparation — 28 October 2026

The current R display layer reproduces the reviewed five main and six supplemental figures. It uses neutral study-design boxes, adjusted label clearance, an evidence-led mixed panel layout and editable PDF/SVG text. The plotted data and 24 frozen input checksums are unchanged. PDF text uses standard Helvetica/Symbol; TIFF/PNG/SVG use Arial. TIFF is 600-dpi RGB with LZW compression, PNG is 300 dpi, and the figure plates are 7.25 inches wide. ASM's [initial submission checklist](https://journals.asm.org/initial-submission-checklist) accepts format-neutral figures; its [figure size guidance](https://journals.asm.org/writing-your-paper) recommends 7 inches as a file-size tip rather than stating a hard maximum.

The renderer includes a portable R viewport measurement helper and a Python standard-library checker for shared physical plot-area boundaries at 1.5 pt. The launcher supplies its own Python executable through `FIGURE_PYTHON`. It does not depend on an author-specific home directory or an author-local tool installation. PDF/SVG exports retain text; no statistics are recomputed. Optional PDF assembly merges the native R figure PDFs, avoiding rasterization of scientific content.

The immutable Zenodo software 1.0.0 snapshot remains the original analysis archive. This display update is a separate GitHub commit; it is not a new Zenodo publication or a change to historical analysis results.
