## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: validate
- Origin Date: 2026-09-11T13:23:59+0800
- Verification Status: VERIFIED
- Version Label: validation_v1

# Remaining validation report

- Source: frozen PRJNA1125274 Step98D object, Step98E displacement output, Step98H results, and frozen cross-cohort taxonomic concordance tables
- Overall Confidence: CAUTION
- Scope: downstream validation only; no FASTQ preprocessing, DADA2, patient mapping, or frozen 98H file was changed

## Statistical findings

| Analysis | Result | Confidence |
|---|---|---|
| Executed 20/3 Bray primary | n=24, mean delta=0.192, 95% CI 0.076 to 0.308, Hedges g_z=0.674 (95% CI 0.219 to 1.129), paired-t p=0.00237, Wilcoxon p=0.0177 | CAUTION: medium effect; paired differences are non-normal, but Wilcoxon corroborates |
| Executed 20/3 Aitchison primary | n=24, mean delta=16.518, 95% CI 4.388 to 28.649, Hedges g_z=0.556 (95% CI 0.115 to 0.997), paired-t p=0.00978, Wilcoxon p=0.0258 | SOLID within analysed 24-person population; medium effect |
| Exclude Eukaryota Bray | n=24, mean delta=0.194, 95% CI 0.076 to 0.312, Hedges g_z=0.672 (95% CI 0.217 to 1.127), paired-t p=0.00243, Wilcoxon p=0.0191 | CAUTION sensitivity with paired-t/Wilcoxon concordance |
| Exclude Eukaryota Aitchison | n=24, mean delta=15.228, 95% CI 4.850 to 25.605, Hedges g_z=0.599 (95% CI 0.153 to 1.045), paired-t p=0.00588, Wilcoxon p=0.014 | SOLID sensitivity |
| Bacteria/Archaea-only Bray | n=24, mean delta=0.194, 95% CI 0.076 to 0.312, Hedges g_z=0.672 (95% CI 0.217 to 1.127), paired-t p=0.00243, Wilcoxon p=0.0191 | CAUTION sensitivity with paired-t/Wilcoxon concordance |
| Bacteria/Archaea-only Aitchison | n=24, mean delta=14.131, 95% CI 4.502 to 23.760, Hedges g_z=0.599 (95% CI 0.153 to 1.045), paired-t p=0.00587, Wilcoxon p=0.0164 | SOLID sensitivity |
| Historical-document 5/2 Aitchison | n=23, mean delta=17.223, 95% CI 3.784 to 30.663, Hedges g_z=0.535 (95% CI 0.087 to 0.983), paired-t p=0.0144, Wilcoxon p=0.0333 | CAUTION: one date-clean patient lost numerically |

## Technical findings

- The executed 20/3 table contained 9373 ASVs; 3977 were labelled Eukaryota and accounted for 8.81% of retained reads.
- ASV IDs, count columns, full-seqtab sequences, and taxonomy rows aligned exactly.
- Removing Eukaryota or retaining only Bacteria/Archaea preserved the positive direction, 95% CI exclusion of zero, paired-t inference, and Wilcoxon inference for both primary metrics.
- Bray paired differences departed from normality (Shapiro-Wilk W=0.8884, p=0.0123); Aitchison paired differences did not show evidence of departure (W=0.9594, p=0.4262). The Bray conclusion is retained because the paired Wilcoxon result is concordant.
- The old 5/2 rule retained 22,031 ASVs. CZM-Aitchison retained 23/24 date-clean primary patients and 28/30 complete cases. The affected complete patients were SAL_4 and SAL_51.
- SAL_51 T1 and SAL_4 T0 produced non-finite CLR rows after CZM in the very sparse 5/2 table. The Bacteria/Archaea-only 5/2 sensitivity restored 24/24 and 30/30 while retaining positive inference.
- Therefore, 20/3 remains the executed primary rule; 5/2 is a documentation-concordant sensitivity with a disclosed numerical limitation, not a replacement primary analysis.

## Complete-case selection

- Depth-qualified T0 patients: 116; complete T0/T1/T2: 30; date-clean primary: 24.
- Complete versus baseline-incomplete hospital distribution: Fisher p=0.186; T0 sequencing depth: Wilcoxon p=0.967.
- Complete-case patients necessarily had more qualifying longitudinal runs; that comparison is structural rather than evidence against selection bias.
- Age, sex, severity, antibiotics, nutrition, comorbidity and outcome are absent from the analysis object, so clinical attrition bias cannot be tested.

## Supplementary Figure S1

- Integrity verdict: PASS. The six displayed pairwise primary values and six displayed common-anchor values are source-identical; the figure code reads distinct rows and does not reuse a panel.
- PRJNA516701 has small differences in eligible within-cohort taxa/effects, but these do not change the pairwise summary values used by the figure. The identical panels are therefore real but add limited visual information.

## Reproducibility

- Method: deterministic numeric rerun from the frozen count/taxonomy/metadata object using the same seed and functions.
- Verdict: REPRODUCIBLE. Four frozen primary/complete metric rows matched in patient count and all compared numeric fields within 1e-10; maximum absolute difference=1.110223e-16.
- File-level byte identity is not expected because Step98I uses a new output schema; the frozen Step98H files were not rewritten.

## Warnings

| Type | Detail |
|---|---|
| Survivorship/attrition | Only 30/116 depth-qualified-baseline patients had complete trajectories; 24/116 entered the primary analysis. |
| Centre heterogeneity | Overall validation must not be restated as successful replication in each hospital. |
| Metadata reconciliation | The public 134-group versus paper 132-person discrepancy still requires an author-adjudicated roster. |
| Governance | Local timestamps identify 20/3 as the executed rule, but do not constitute immutable public preregistration. |
| Multiplicity | Sensitivity p-values are robustness diagnostics, not independent discovery tests. |
| Distributional assumption | Bray paired differences were non-normal; interpret the paired t-test together with the concordant Wilcoxon sensitivity. |

## Fallacy scan

- Coverage: 11/11 statistical fallacy types checked.
- Detailed findings: STATISTICAL_VALIDATION_FALLACY_SCAN.csv.

## Final verdict

The main PRJNA1125274 ecological-trajectory conclusion survives every same-population non-target-feature sensitivity and the prokaryote-only 5/2 sensitivity. Confidence remains CAUTION rather than unrestricted SOLID because complete-case attrition, centre heterogeneity, unavailable clinical covariates, non-public threshold governance, and the unresolved 132-versus-134 roster issue limit generalisation.

The remaining unresolved items cannot be answered by further computation on the current files: an author-adjudicated 132-person clinical roster and clinical covariates for attrition/confounding assessment.
