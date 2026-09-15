# Study data access and analysis checkpoints

## Public cohorts used in the manuscript

### NCBI SRA / BioProject

- PRJNA691455
- PRJNA516701
- PRJNA851469
- PRJNA578267
- PRJNA430161
- PRJNA1166732
- PRJNA978257
- PRJNA1010969
- PRJNA1125274

### ENA

- PRJEB82425

### NGDC-GSA

- CRA002354

The sample/run-level inclusion set must be taken from the frozen metadata and manifests in the final results, rather than downloading every run associated with each project accession without filtering.

## Local raw and analysis-ready data

- Raw/publicly downloaded data root: `data/`
- Canonical analysis-ready hard-link workspace: `data/_V2_ANALYSIS_READY/`
- Final repaired frozen metadata: `results/V2_21D_FINAL_FREEZE_CLINICAL_REPAIRED/`
- Final analysis objects: `results/V2_27B_ANALYSIS_OBJECTS_AND_REPLICATE_FREEZE/`
- Final depth-filtered analysis-ready state: `results/V2_26B_DEPTH2000_FINAL_ANALYSIS_READY/`

`data/_V2_ANALYSIS_READY` contains hard links. Do not delete or move it as though it were an independent disposable copy.

## Reviewer-facing recommendation

Deposit the following in a repository with a permanent DOI before submission or revision:

1. Code snapshot and run order.
2. Final frozen metadata and data dictionary, after checking that no restricted/private fields are present.
3. Small processed analysis objects or flat tables needed to rerun manuscript statistics.
4. Final figures/tables and SHA256 manifest.
5. Accession-to-sample/run inclusion manifest.

The repository name, release tag, URL and DOI remain author-supplied fields and must be inserted into the manuscript, cover letter and Data/Code Availability statement.
