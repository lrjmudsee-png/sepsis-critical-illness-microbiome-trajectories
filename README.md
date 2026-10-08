# Longitudinal gut-microbiome trajectories in sepsis and critical illness

Code-availability repository for the study of longitudinal gut-microbiome ecological displacement and heterogeneous taxonomic trajectories across sepsis and critical-illness cohorts.

**Start with [START_HERE.md](START_HERE.md).** All 364 distributed programs and launchers now have descriptive numbered names. [The filename map](docs/CODE_FILENAME_MAP.csv) connects the old names to the new names; historical data/output labels remain unchanged. Numbered inventory order does not mean every historical alternative should be executed.

## License and current availability

Original code and associated software documentation are licensed under the [MIT License](LICENSE), copyright (c) 2026 Lu Rongji. This does not relicense public source data, research outputs or third-party software; see [licensing scope](docs/LICENSING.md).

This is a code-and-audit repository, not yet a self-contained reproduction package for every manuscript result. Selected statistical runners require processed inputs distributed with the manuscript as Supplementary Data 5; those inputs are not duplicated in this code-only repository. The public aggregate forest example can be rendered using the command in [selected-statistics reproduction](docs/SELECTED_STATISTICS_REPRODUCTION.md). Software citation metadata are provided in `CITATION.cff`. Software version `1.0.0`, from commit `1e998aaa9b19e48643996c95d78b7cf7b74b1807`, is published at [Zenodo DOI 10.5281/zenodo.23134779](https://doi.org/10.5281/zenodo.23134779). The archive is immutable; later documentation updates on the default branch do not change this archived software snapshot.

## Repository contents

- `analysis/`: Retained historical source scripts plus the secondary direction-analysis snapshot and PRJNA1125274 validation/governance scripts; not a guarantee that every historical script was recovered.
- `docs/RUN_ORDER.csv`: ordered map from metadata recovery and sequence processing to statistics and manuscript assembly.
- `docs/DATA_ACCESS.md`: public cohort accessions and local checkpoint definitions.
- `docs/ENVIRONMENT.md`: verified software and package versions.
- `evidence/`: records from the independently isolated terminal reproduction run performed on 2026-09-04, plus the downstream PRJNA1125274 validation audit completed on 2026-09-11.
- `tools/05_run_prepared_terminal_chain.ps1`: runner for the 11-step frozen-input-to-manuscript terminal chain.

Raw sequencing reads, patient-level data, analysis workspaces, generated result trees, caches, and local directory links are not included in this repository.

## Public data

The study uses publicly available cohorts from NCBI SRA/BioProject, ENA, and NGDC-GSA:

`PRJNA691455`, `PRJNA516701`, `PRJNA851469`, `PRJNA578267`, `PRJNA430161`, `PRJNA1166732`, `PRJNA978257`, `PRJNA1010969`, `PRJNA1125274`, `PRJEB82425`, and `CRA002354`.

Use the frozen sample/run inclusion manifests rather than downloading every run under each project accession. See `docs/DATA_ACCESS.md`.

## Verified terminal reproduction

The frozen-input terminal chain was rerun in an isolated result root on Windows 11 with R 4.4.0.

- 11/11 analysis and manuscript-integration steps exited successfully.
- 132/132 generated CSV files matched the frozen results after run-root normalization; 112 were identical without normalization.
- 10/10 PDFs were identical after rendering and SHA256 comparison.
- 8/8 PNG files were SHA256-identical.
- 42/42 text and completion-marker files matched after normalizing paths and timestamps.

See `evidence/REPRODUCTION_REPORT_20260904.md` and the machine-readable comparison tables in `evidence/`.

## PRJNA1125274 downstream validation

The Step98 series adds the external longitudinal cohort, patient-metadata reconciliation, frozen-result governance, non-target-feature and ASV-threshold sensitivities, complete-case selection audit, statistical-assumption checks, and reproducibility verification. The final report and machine-readable summaries are in `evidence/PRJNA1125274_20260911/`.

These downstream checks did not rerun FASTQ preprocessing or DADA2. The author-adjudicated 132-person clinical roster and clinical covariates needed for a definitive 132-versus-134 resolution and clinical attrition/confounding analysis are not available in this repository.

## Reproduction levels

### 1. Frozen inputs to manuscript outputs

This route was tested with prepared frozen metadata and analysis objects; those full inputs are not supplied by cloning this repository alone. It runs the following terminal chain:

1. trajectory ecology robustness;
2. longitudinal evidence freeze;
3. external validation evidence freeze;
4. infection-source evidence freeze;
5. cross-cohort taxonomic reproducibility;
6. global evidence synthesis;
7. manuscript source-pack assembly;
8. manuscript architecture and main tables;
9. Table 1 population correction and Results assembly;
10. final Results/Table 1 sentence audit;
11. supplementary tables and Discussion assembly.

The PowerShell runner expects a prepared Windows run root containing `_scripts/`, `data/`, `metadata/`, `processed_data/`, and the required frozen upstream `results/` directories.

```powershell
powershell -ExecutionPolicy Bypass -File tools/05_run_prepared_terminal_chain.ps1 -RunRoot E:\path\to\prepared_run_root
```

### 2. Raw public reads to manuscript outputs

The retained historical workflow is mapped in `docs/RUN_ORDER.csv`. It includes metadata reconstruction, DADA2 processing, the VSEARCH 97% OTU route for CRA002354, SILVA 138.2 taxonomy assignment, analysis-object freezing, longitudinal models, robustness analyses, and manuscript integration.

This full route has not been rerun after packaging and remains partly Windows-path-bound. Several historical scripts also install packages automatically. A future software version should parameterize the remaining paths and lock the environment with `renv` or a container.

Historical 96F2-F11 cosmetic figure scripts have not been recovered. For current submission figure styling, use the corrected final display layer in `analysis/final_figures/`, documented in [final figure reproduction](docs/FINAL_FIGURE_PIPELINE.md). Its external frozen CSV inputs are not included in this code-only repository.

## Important provenance note

The repository preserves superseded and repair scripts because they document how the final frozen workflow was reached. Use `START_HERE.md`, the current `docs/RUN_ORDER.csv` catalogue and `docs/CODE_FILENAME_MAP.csv` when locating the renamed sources. Historical alternatives and superseded revisions are retained for provenance, not selected automatically as authoritative outputs.

## Code availability statement

The 4 October 2026 code update adds the actual corrected submission figure programs, explicit input/output arguments, mandatory input checksum verification and output-overwrite protection. Run `tools/04_reproduce_final_figures.py` to render all five main and six supplementary figures from the 24 frozen CSV inputs. This is a rendering route, not a new statistical analysis or raw-read reprocessing. See [final figure reproduction](docs/FINAL_FIGURE_PIPELINE.md) for input-access limits, environment versions and commands. The submission snapshot is archived as Zenodo version 1.0.0 at the DOI above; no new software archive is created by the 8 October documentation update.

Analysis code, run order, environment information, public data accessions, and machine-readable terminal reproduction checks for this study are available at <https://github.com/lrjmudsee-png/sepsis-critical-illness-microbiome-trajectories>. Raw sequence data remain available from their originating public archives under the accessions listed above.

See [CODE_AVAILABILITY.md](CODE_AVAILABILITY.md) for the full statement and its reproduction limits. `CODE_FILE_MANIFEST.csv` records the current distributed files, excluding itself and `.git/`; `MANIFEST_SHA256.csv` is a historical snapshot and is not the current manifest.


## Secondary direction analysis and selected-statistics reproduction

The `analysis/direction_analysis/` snapshot contains the retrospective fixed family-balance analysis and healthy-reference audit from 13 September 2026. Both Holm-adjusted balance tests were inconclusive (p=1); no common pathobiome direction or host mechanism was established. Aggregate result and QC records are in `evidence/DIRECTION_ANALYSIS_20260913/`. Historical runners in that directory retain their original Windows/frozen-input paths.

`tools/01_reproduce_selected_statistics.R` accepts relative flat-file input/output directories and provides a smaller tested statistical route. Required processed inputs are supplied with the manuscript as Supplementary Data 5 and are not duplicated in this public code-only repository. Thus the repository alone does not presently rerun every manuscript result. See `docs/PROCESSED_INPUT_AVAILABILITY.md`. The full raw-processing historical route remains separately documented and was not rerun in this update.


## Final specimen governance and portable tests

The 2026-09-15 post-hoc PRJNA516701 audit found mixed stool/rectal specimens and five of 15 pairs with a type change at an available visit. Aggregate same-specimen and cohort-exclusion sensitivity results, plus 98 selected-statistics checks and isolated ZIP-run checks, are in `evidence/FINALIZATION_20260915/`. The original primary sets/results were preserved; sensitivities are nominal, not replacement primary tests. `tools/02_reproduce_specimen_sensitivity.R` requires the documented metafor environment. See `docs/SELECTED_STATISTICS_REPRODUCTION.md`. The repository now includes an active `CITATION.cff` aligned with the manuscript CRediT statement, with Rongji Lu listed as the software creator. Software version `1.0.0`, from commit `1e998aaa9b19e48643996c95d78b7cf7b74b1807`, is published at [Zenodo DOI 10.5281/zenodo.23134779](https://doi.org/10.5281/zenodo.23134779). The archive is immutable; later documentation updates on the default branch do not change this archived software snapshot. MIT software licensing does not resolve the separate data-redistribution review.

## Reporting addenda supplied with the manuscript (8 October 2026)

Supplementary Data 5 also supplies `REPORT_POOLED_EFFECT_TESTS.R` and `REPORT_DISCOVERY_SIGNED_RANK.R`, with their aggregate reporting tables and task-specific R session record. These downstream reporting-recovery helpers are not included in the archived 1.0.0 snapshot. They do not rerun FASTQ/DADA2 or change the frozen primary results. Cite the fixed software DOI above and retain the manuscript supplement when reproducing these addenda.

