# Release notes

## V2 code-availability release — 2026-09-04

- Curated 326 executable source files from the complete V2 analysis history.
- Excluded raw/processed data, directory links, caches, R session files, duplicate archives, and generated figure assets.
- Included public-accession, environment, run-order, result-status, and static-QA documentation.
- Included an isolated 11-step terminal reproduction report and machine-readable comparisons.
- Fixed Windows R source-encoding handling by clearing incompatible `C.UTF-8` environment variables and passing `--encoding=UTF-8`.
- Fixed the Step95A2 source locator so later manuscript/reproducibility directories cannot make its candidate audit self-referential.

