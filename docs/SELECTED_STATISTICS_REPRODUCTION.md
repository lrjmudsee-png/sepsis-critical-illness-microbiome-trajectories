# Selected statistical reproduction

The manuscript's Supplementary Data 5 supplies ten flat CSV inputs, a field dictionary, an input SHA256 manifest and two selected-statistics command-line entries. Obtain that manuscript supplement (including its updated manifest) rather than infer or fabricate absent patient values. These processed inputs are not duplicated in the code-only repository. Public raw-data accessions alone do not replace the required processed inputs.

    Rscript tools/01_reproduce_selected_statistics.R INPUT_DIRECTORY NEW_OUTPUT_DIRECTORY
    Rscript tools/02_reproduce_specimen_sensitivity.R INPUT_DIRECTORY NEW_SPECIMEN_OUTPUT_DIRECTORY

The first requires base R only. The second requires metafor (tested R4.4.0, metafor5.0-1; dependencies recorded in the sensitivity session text). Neither installs packages, accesses a network, overwrites existing output, runs FASTQ/DADA2 or reconstructs distances from ASVs. The selected route reruns external 24/30 paired contrasts, B from raw family sums, its REML/conservative-HK/Holm family, and the ordinary healthy-reference contrast. Two-level healthy bootstrap and full raw-to-manuscript workflow have separate historical scopes.

On 2026-09-15 the ZIP was extracted under a different directory and both entries actually executed. All 98 numeric comparisons and 21 package/run checks passed. Three specimen-sensitivity CSVs matched their initial rerun exactly; only aggregate tables are public. The primary specimen audit is post-hoc: original frozen results remain unchanged, and all subset tests are nominal. These are isolated-directory tests on the same machine, not an independent-machine full raw-data validation.

The original software is now MIT-licensed; see [licensing scope](LICENSING.md). Software version `1.0.0`, from commit `1e998aaa9b19e48643996c95d78b7cf7b74b1807`, is published at [Zenodo DOI 10.5281/zenodo.23134779](https://doi.org/10.5281/zenodo.23134779). The archive is immutable; later documentation updates on the default branch do not change this archived software snapshot. Active `CITATION.cff` lists Rongji Lu as software creator. Any further processed-data deposit requires its own redistribution/licensing review.

## Direction forest layout repair

    Rscript tools/03_render_direction_forest.R evidence/FINALIZATION_20260915/B14_forest_plot_source.csv NEW_FOREST.pdf

This base-R entry only draws the five frozen aggregate estimates and pointwise intervals; it does not refit statistics. The submission S5 margin/label repair was reproduced with identical rendered pixels (PDF metadata timestamps may differ). The historical 99B snapshot includes the corresponding layout change; original frozen result tables were not changed. Historical 96F2-F11 cosmetic figure scripts are not recovered by this entry.

## Statistical-reporting addenda (8 October 2026)

Supplementary Data 5 adds `REPORT_POOLED_EFFECT_TESTS.R` and `REPORT_DISCOVERY_SIGNED_RANK.R`, their aggregate CSV reports and an R session record. The first recovers the existing REML/Hartung-Knapp pooled-effect tests, not the heterogeneity Q tests; the second reports retained paired signed-rank statistics. These helpers do not overwrite primary results or rerun read processing. They are distributed in the manuscript supplement and are not part of archived software 1.0.0. The documented 98-check selected-statistics route remains a limited reproducibility test, not complete coverage of every exploratory genus screen.

