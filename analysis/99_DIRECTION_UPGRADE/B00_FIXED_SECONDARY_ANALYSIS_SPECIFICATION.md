# Fixed secondary-analysis specification for stage B

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: plan-to-run handoff
- Freeze date: 2026-09-13 Asia/Hong_Kong
- Verification status: FIXED_BEFORE_STAGE_B_EFFECT_ESTIMATION
- Version: direction_stage_b_sap_v1
- Scope statement: This is a prospectively fixed specification for a new secondary analysis after prior study results and the stage-A feasibility audit were known. It is not a historical preregistration and must not be described as one.

## Primary question

Does a fixed opportunist-associated versus commensal-associated Family balance increase from each study's prespecified early to late follow-up in the three natural-history cohorts, with directionally consistent evidence in the frozen 24-patient PRJNA1125274 primary population?

## Fixed taxonomy dictionary

- P: exact Family label `Enterobacteriaceae` or `Enterococcaceae`.
- C: exact Family label `Lachnospiraceae`, `Ruminococcaceae`, or `Oscillospiraceae`.
- OTHER: every other label, including unresolved Family labels.
- Automated normalization is limited to trimming whitespace and removing a leading `f__` or `F__`.
- No synonym expansion, genus rescue, abundance filter, prevalence filter, or result-dependent revision is permitted.

## Fixed populations and contrasts

- PRJNA691455: frozen ALL_PAIRED, Day 3 versus Day 7; expected maximum 9 pairs.
- PRJNA851469: frozen ALL_PAIRED, Day-3 versus Day-7; expected maximum 14 pairs.
- PRJNA516701: frozen ALL_PAIRED, DAY_3 versus DAY_7; expected maximum 15 pairs. COMMON_ANCHOR is separate.
- PRJNA1125274 primary: `in_prespecified_primary`, T1 versus T2; expected maximum 24 patients. The 30 complete T0/T1/T2 patients form a sensitivity population only.
- A pair is not replaced when the new metric is unavailable.

## Fixed metric

For each sample, aggregate integer counts before adding a pseudocount:

`B = ln[(P_count + 0.5)/(C_count + 0.5)]`

For each patient:

`delta_B = B(late) - B(early)`

A positive value means that P increased relative to C. It is not evidence of absolute bacterial expansion, metabolite production, translocation, or causation.

- If P+C=0 at either contrast sample, B is uninterpretable and the pair is missing. `ln(0.5/0.5)=0` is not used as an ecological value.
- One-sided zero is calculable under the pseudocount rule and is disclosed.
- Fixed sensitivities: pseudocount 0.1, pseudocount 1, and no pseudocount among samples with both P and C positive.

## Cohort estimates

- Primary estimate: mean paired `delta_B` and two-sided 95% t interval.
- Also report n, SD, SE, median, Cohen dz, two-sided one-sample t p value, two-sided paired Wilcoxon p value, Shapiro-Wilk p value, and a 5,000-replicate patient bootstrap percentile interval.
- Random seed: 20260913 plus a deterministic analysis-label offset.
- With n<3 or zero variance, do not force standard inference.

## Natural-history meta-analysis

- Include only PRJNA691455, PRJNA851469, and PRJNA516701.
- Pool raw mean `delta_B` with within-study variance `SD(delta_B)^2/n`.
- Random-effects REML.
- Hartung-Knapp inference using the larger of the conventional random-effects SE and the raw HK SE, with t critical value on k-1 degrees of freedom.
- Report pooled mean, 95% CI, p value, tau-squared, I-squared, Cochran Q, Q p value, and weights.
- Fixed sensitivities: common-anchor population, pseudocount 0.1, pseudocount 1, no-pseudocount nonzero subset, and leave-one-cohort-out. k=2 leave-one-out results are instability diagnostics only.

## External test and multiplicity

- The external primary estimate uses the frozen 24-patient population.
- Sensitivities: 30 complete cases, SAL_8 exclusion membership, hospital-specific primary subsets, pseudocount 0.1, pseudocount 1, and nonzero-only calculation.
- The primary family contains exactly two p values: the three-cohort natural-history meta-analysis and the 24-patient external test.
- Apply Holm adjustment to these two p values. Confidence intervals are pointwise and are not simultaneous-coverage intervals.
- Strong common-direction support requires both point estimates to be positive and both Holm-adjusted p values below 0.05. Otherwise report partial support, directional heterogeneity, or no support without revising the dictionary.

## Healthy-reference analysis

- Use the 13 healthy samples already co-processed in the frozen PRJNA851469 object.
- Compute Bray-Curtis on within-object relative abundance.
- For each patient sample, H is the mean Bray-Curtis distance to all 13 healthy references.
- Evaluate Day-7 minus Day-3 H in the frozen paired population.
- Use a 5,000-replicate two-level bootstrap resampling healthy individuals and patient pairs to propagate reference and patient uncertainty.
- This is a secondary same-cohort reference, not an independent validation.

## Interpretive analyses

- Join the existing frozen Bray displacement difference `delta_D` to `delta_B` by exact project, patient, and run identifiers.
- Report Spearman associations by cohort. These are same-data-source consistency analyses, not independent validation.
- Define descriptive quadrants only by the signs of delta_D and delta_B. Do not optimize thresholds or train a patient classifier.
- Retain all negative and opposing-direction results.

## Prohibited actions

- Do not rerun DADA2, FASTQ preprocessing, Step98D, taxonomy assignment, or patient reconciliation.
- Do not alter the frozen patient map, original objects, A-stage outputs, `V2-SUBMISSION`, GitHub, or Zenodo.
- Do not substitute EII as an independent outcome.
- Do not infer patient-level host measurements from published correlation matrices.
- Do not add tests or alter populations because of observed p values.

## Human ratification status

The user authorized automatic continuation. These defaults therefore govern computation. Final manuscript promotion, biological naming, causal language, author contact, and acceptance of the analysis as a paper claim remain human decisions and will be listed separately.

