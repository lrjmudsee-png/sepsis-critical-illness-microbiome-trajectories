# Stage-B statistical validation and governance report

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: validation and interpretation
- Verification Status: VALIDATED_WITH_CAUTION
- Version: direction-stage-B-governance-v1

## Confirmed numerical result

- Natural-history REML/Hartung-Knapp: estimate -0.464, 95% CI -3.463 to 2.536, raw p=0.5743, Holm p=1.0000.
- Frozen PRJNA1125274 primary population (n=24): estimate 0.287, 95% CI -0.956 to 1.529, raw p=0.6378, Holm p=1.0000.
- Rule-based conclusion: `DIRECTIONALLY_HETEROGENEOUS`.
- PRJNA851469 healthy-reference secondary result: mean change 0.006, ordinary 95% CI -0.026 to 0.037, two-level bootstrap 95% CI -0.028 to 0.039; patients moved away from the same-cohort healthy reference on average.

## Statistical fallacy audit

1. **Simpson's paradox — checked.** External center-specific means were compared with the overall mean. All-center reversal detected: FALSE. This does not eliminate other center heterogeneity.
2. **Ecological fallacy — caution.** Published immune-mediator tables are aggregate family-correlation matrices and cannot substitute for patient-level host measurements.
3. **Berkson's bias — caution.** ICU/sepsis cohort selection may induce associations that do not generalize outside enrolled populations.
4. **Collider bias — limited exposure.** No new covariate-adjusted host model was fitted, avoiding unsupported adjustment; selection into paired complete cases can still act as a collider.
5. **Base-rate neglect — not a diagnostic classifier analysis.** No claims about prediction or post-test probability are made.
6. **Regression to the mean — monitored.** Change was defined between fixed visits rather than selecting individuals by extreme baseline balance, but two-time-point change remains noisy.
7. **Survivorship/complete-case bias — caution.** Paired follow-up requires observed later samples. PRJNA516701 additionally excludes three double-zero balance pairs under the fixed estimand.
8. **Look-elsewhere effect — partly controlled.** The P/C family dictionary and two-test Holm family were fixed before examining new effect estimates; secondary correlations and quadrant summaries remain exploratory.
9. **Garden of forking paths — caution.** This is a time-stamped post-hoc secondary SAP, not preregistration. Pseudocount, common-anchor, leave-one-out, center, and complete-case results must remain labelled sensitivity analyses.
10. **Correlation-causation fallacy — caution.** The balance is a compositional direction marker, not evidence of absolute bacterial expansion or mechanism.
11. **Reverse causality — caution.** Treatment, illness severity, feeding, bowel function, and recovery may drive the microbial pattern.

## Governance decision

The analysis is acceptable as a reproducible secondary direction-of-change layer if all FAIL counts remain zero. It cannot be promoted to patient-level host-mechanism validation without joinable host data. The healthy-reference analysis is same-cohort supportive context, not an independent cohort.
