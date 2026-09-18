# Selected statistical reproduction

The local author-review input bundle has 10 flat CSVs, a field dictionary, input SHA256 manifest and two command-line entries. Inputs are not publicly redistributed pending author approval. Obtain the author-approved bundle rather than infer or fabricate absent patient values. Public archive accessions alone are not the missing processed inputs.

    Rscript tools/REPRODUCE_SELECTED_STATISTICS.R INPUT_DIRECTORY NEW_OUTPUT_DIRECTORY
    Rscript tools/REPRODUCE_SPECIMEN_SENSITIVITY.R INPUT_DIRECTORY NEW_SPECIMEN_OUTPUT_DIRECTORY

The first requires base R only. The second requires metafor (tested R4.4.0, metafor5.0-1; dependencies recorded in the sensitivity session text). Neither installs packages, accesses a network, overwrites existing output, runs FASTQ/DADA2 or reconstructs distances from ASVs. The selected route reruns external 24/30 paired contrasts, B from raw family sums, its REML/conservative-HK/Holm family, and the ordinary healthy-reference contrast. Two-level healthy bootstrap and full raw-to-manuscript workflow have separate historical scopes.

On 2026-09-15 the ZIP was extracted under a different directory and both entries actually executed. All 98 numeric comparisons and 21 package/run checks passed. Three specimen-sensitivity CSVs matched their initial rerun exactly; only aggregate tables are public. The primary specimen audit is post-hoc: original frozen results remain unchanged, and all subset tests are nominal. These are isolated-directory tests on the same machine, not an independent-machine full raw-data validation.

The original software is now MIT-licensed; see [licensing scope](LICENSING.md). Remaining archival preparation includes confirmed software creators and version metadata, plus a separately approved redistribution/licensing plan for any processed data to be deposited. An unpublished article does not yet need a final article DOI for software archiving. `CITATION.cff.template` is not a validated active CFF, and Zenodo remains deferred at the maintainer's request.

## Direction forest layout repair

    Rscript tools/RENDER_DIRECTION_FOREST.R evidence/FINALIZATION_20260915/B14_forest_plot_source.csv NEW_FOREST.pdf

This base-R entry only draws the five frozen aggregate estimates and pointwise intervals; it does not refit statistics. The submission S5 margin/label repair was reproduced with identical rendered pixels (PDF metadata timestamps may differ). The historical 99B snapshot includes the corresponding layout change; original frozen result tables were not changed. Historical 96F2-F11 cosmetic figure scripts are not recovered by this entry.
