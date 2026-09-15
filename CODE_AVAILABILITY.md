# Code and data availability

## Manuscript-ready statement

All analysis code supporting this study is available at <https://github.com/lrjmudsee-png/sepsis-critical-illness-microbiome-trajectories>, including scripts for public-metadata reconstruction, sequence processing, taxonomy assignment, longitudinal ecological and taxonomic analyses, robustness analyses, figure and table generation, and manuscript-result assembly. The repository also contains the ordered workflow, verified software environment, public data accessions, and machine-readable evidence from an isolated frozen-input terminal rerun in which all 11 final steps completed successfully and all 192 generated files passed content or rendered-output comparison against the frozen results.

The repository additionally contains the Step98 PRJNA1125274 external-validation and governance workflow. This downstream audit reproduced the frozen primary statistics, tested ASV-threshold and non-target-feature sensitivity, checked complete-case selection and model assumptions, and documented the unresolved public-metadata discrepancy between 134 reconstructed identifiers and the paper-reported 132 patients. It did not rerun FASTQ preprocessing or DADA2.

Raw sequencing data were obtained from public repositories under accessions PRJNA691455, PRJNA516701, PRJNA851469, PRJNA578267, PRJNA430161, PRJNA1166732, PRJNA978257, PRJNA1010969, PRJNA1125274, PRJEB82425, and CRA002354. The exact sample/run inclusion set should be taken from the study's frozen inclusion manifests rather than inferred from all runs associated with an accession.

Large raw-read files and patient-level data are not redistributed in the code repository. The tested repository-level reproduction route begins from frozen metadata and processed analysis objects. The full historical raw-read workflow is included for provenance but remains partly dependent on the original Windows directory layout.
