# Verified analysis environment

## Final figure/table package session

Source: `results/V2_97A_FINAL_FIGURE_AND_TABLE_PACKAGE/04_PROVENANCE_AND_QC/sessionInfo_STEP97A.txt`

- Operating system: Windows 11 x64, build 26200
- Time zone: Asia/Hong_Kong
- R: 4.4.0 (x86_64-w64-mingw32/x64)
- DADA2: 1.34.0
- ggplot2: 4.0.3
- dplyr: 1.2.1
- readr: 2.2.0
- tibble: 3.3.1
- stringr: 1.6.0
- patchwork: 1.3.2
- openxlsx: 4.2.8.1

## Statistical packages verified in the local R 4.4.0 library

- vegan: 2.7.3
- lme4: 2.0.1
- lmerTest: 3.2.1
- permute: 0.9.10

## External software and references

- VSEARCH: 2.31.0
- SILVA: release 138.2 DADA2 training/reference files
- Main paired-read route: DADA2
- CRA002354 route: VSEARCH 97% OTU workflow because the public FASTA files lack per-base quality scores

## Environment limitations still to resolve before public code release

- No project-level `renv.lock`, Conda lock file or container image is currently present.
- Python metadata-processing dependencies were not captured at the time of every historical step.
- Some package logs report packages built under a later R patch version than the R 4.4.0 runtime. Final output comparisons did not show numerical differences, but a locked release environment is still recommended.
- Exact executable paths are Windows-specific in the historical scripts.
