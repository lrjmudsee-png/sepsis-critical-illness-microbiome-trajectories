# Release notes

## Repository-only licensing update 2026-09-18

- Added the standard MIT License for original software and associated documentation, with the maintainer-confirmed copyright holder Lu Rongji.
- Documented the separate status of third-party software, public research data, processed inputs and research outputs; no data license was assigned.
- Aligned README, code/data availability and citation-template notes. Removed overbroad claims of complete script recovery and clarified the prepared-input scope of prior reproduction checks.
- Refreshed `CODE_FILE_MANIFEST.csv`; the original `MANIFEST_SHA256.csv` remains a historical snapshot.
- No scientific scripts, input data or frozen statistical results were changed or rerun. The creator list remains unconfirmed; no active CFF citation, GitHub release or Zenodo publication was created.

## Repository-only finishing update 2026-09-15

Added the 99 direction snapshot, two parameterized statistical entries, aggregate negative results, post-hoc specimen-consistency/cohort-exclusion sensitivity and isolated-directory test evidence. Original scientific inputs/primary results remain frozen. Processed patient inputs are local author-review only, software citation is an inactive template, and no GitHub release or Zenodo publication was created.


## Post-release validation update — 2026-09-11

- Added 16 Step98 source scripts covering PRJNA1125274 external validation, patient-level metadata reconciliation, final governance, and remaining downstream checks.
- Added a curated evidence bundle for ASV-threshold and non-target-feature sensitivity, complete-case selection, statistical assumptions, frozen-result numeric reproduction, Supplementary Figure S1 source integrity, and an 11-item statistical fallacy scan.
- Confirmed that the frozen 20/3 primary results are numerically reproducible and remain directionally significant after excluding Eukaryota or restricting to Bacteria/Archaea.
- Documented the 5/2 all-feature CZM row-loss limitation and the unresolved need for an author-adjudicated 132-person roster and clinical covariates.
- This is a repository update only; no GitHub release or Zenodo version was created.

## Code-availability snapshot — 2026-09-04

- Curated 326 executable source files from the complete V2 analysis history.
- Excluded raw/processed data, directory links, caches, R session files, duplicate archives, and generated figure assets.
- Included public-accession, environment, run-order, result-status, and static-QA documentation.
- Included an isolated 11-step terminal reproduction report and machine-readable comparisons.
- Fixed Windows R source-encoding handling by clearing incompatible `C.UTF-8` environment variables and passing `--encoding=UTF-8`.
- Fixed the Step95A2 source locator so later manuscript/reproducibility directories cannot make its candidate audit self-referential.
