# Longitudinal gut-microbiome trajectories in sepsis and critical illness

Code-availability repository for the study of longitudinal gut-microbiome ecological displacement and heterogeneous taxonomic trajectories across sepsis and critical-illness cohorts.

## License and current availability

Original code and associated software documentation are licensed under the [MIT License](LICENSE), copyright (c) 2026 Lu Rongji. This does not relicense public source data, research outputs or third-party software; see [licensing scope](docs/LICENSING.md).

This is a code-and-audit repository, not yet a self-contained reproduction package for every manuscript result. Selected statistical runners require processed inputs currently held in a local author-review bundle. The public aggregate forest example can be rendered using the command in [selected-statistics reproduction](docs/SELECTED_STATISTICS_REPRODUCTION.md). Software citation metadata and Zenodo publication remain pending; no new archival release is made by this update.

## Repository contents

- `analysis/`: Retained historical source scripts plus the secondary direction-analysis snapshot and PRJNA1125274 validation/governance scripts; not a guarantee that every historical script was recovered.
- `docs/RUN_ORDER.csv`: ordered map from metadata recovery and sequence processing to statistics and manuscript assembly.
- `docs/DATA_ACCESS.md`: public cohort accessions and local checkpoint definitions.
- `docs/ENVIRONMENT.md`: verified software and package versions.
- `evidence/`: records from the independently isolated terminal reproduction run performed on 2026-09-04, plus the downstream PRJNA1125274 validation audit completed on 2026-09-11.
- `tools/V2_REPRODUCE_TERMINAL_CHAIN.ps1`: runner for the 11-step frozen-input-to-manuscript terminal chain.

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
powershell -ExecutionPolicy Bypass -File tools/V2_REPRODUCE_TERMINAL_CHAIN.ps1 -RunRoot E:\path\to\prepared_run_root
```

### 2. Raw public reads to manuscript outputs

The retained historical workflow is mapped in `docs/RUN_ORDER.csv`. It includes metadata reconstruction, DADA2 processing, the VSEARCH 97% OTU route for CRA002354, SILVA 138.2 taxonomy assignment, analysis-object freezing, longitudinal models, robustness analyses, and manuscript integration.

This full route has not been rerun after packaging and remains partly Windows-path-bound. Several historical scripts also install packages automatically. A future archival release should parameterize the remaining paths and lock the environment with `renv` or a container.

Historical 96F2-F11 cosmetic figure scripts have not been recovered. The smaller public forest-rendering entry does not replace these missing scripts or reproduce all final figure styling.

## Important provenance note

The repository preserves superseded and repair scripts because they document how the final frozen workflow was reached. Use `docs/RUN_ORDER.csv`, the `FINAL`, `FIXED`, and freeze-labelled scripts, and the 11-step terminal chain above when identifying authoritative outputs.

## Code availability statement

Analysis code, run order, environment information, public data accessions, and machine-readable terminal reproduction checks for this study are available at <https://github.com/lrjmudsee-png/sepsis-critical-illness-microbiome-trajectories>. Raw sequence data remain available from their originating public archives under the accessions listed above.

See [CODE_AVAILABILITY.md](CODE_AVAILABILITY.md) for the full statement and its reproduction limits. `CODE_FILE_MANIFEST.csv` records the current distributed files, excluding itself and `.git/`; `MANIFEST_SHA256.csv` is a historical snapshot and is not the current manifest.


## Secondary direction analysis and selected-statistics reproduction

The `analysis/99_DIRECTION_UPGRADE/` snapshot contains the retrospective fixed family-balance analysis and healthy-reference audit from 13 September 2026. Both Holm-adjusted balance tests were inconclusive (p=1); no common pathobiome direction or host mechanism was established. Aggregate result and QC records are in `evidence/DIRECTION_ANALYSIS_20260913/`. Historical runners in that directory retain their original Windows/frozen-input paths.

`tools/REPRODUCE_SELECTED_STATISTICS.R` accepts relative flat-file input/output directories and provides a smaller tested statistical route. Required processed inputs are staged locally for author review and are not redistributed in this public code-only repository. Thus the repository alone does not presently rerun every manuscript result. See `docs/PROCESSED_INPUT_AVAILABILITY.md`. The full raw-processing historical route remains separately documented and was not rerun in this update.


## Final specimen governance and portable tests

The 2026-09-15 post-hoc PRJNA516701 audit found mixed stool/rectal specimens and five of 15 pairs with a type change at an available visit. Aggregate same-specimen and cohort-exclusion sensitivity results, plus 98 selected-statistics checks and isolated ZIP-run checks, are in `evidence/FINALIZATION_20260915/`. The original primary sets/results were preserved; sensitivities are nominal, not replacement primary tests. `tools/REPRODUCE_SPECIMEN_SENSITIVITY.R` requires the documented metafor environment. See `docs/SELECTED_STATISTICS_REPRODUCTION.md`. No new GitHub release or Zenodo DOI was made; `CITATION.cff.template` is inactive pending confirmed software creators and citation metadata. MIT software licensing does not resolve the separate data-redistribution review.
