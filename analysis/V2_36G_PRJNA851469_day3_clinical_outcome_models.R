
# ============================================================
# V2_36G PRJNA851469 DAY-3 LANDMARK CLINICAL OUTCOME MODELS
#
# Scientific question:
# Does Day-1 -> Day-3 within-patient ecological displacement predict
# subsequent nosocomial infection or death through Day 30?
#
# Design:
# - Day-3 landmark
# - Exclude infection/death on or before Day 3
# - Outcome: first subsequent nosocomial infection or death
# - Primary exposure: standardized Day-3 Bray displacement
# - Primary model: Cox proportional hazards
# - Sensitivity: admission-antibiotic adjustment
# - Additional sensitivity: sepsis-admission adjustment
# - Logistic model as horizon-based sensitivity
# - Bootstrap + permutation robustness
#
# IMPORTANT:
# This analysis is exploratory/supportive.
# No previously frozen inference is changed.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(survival)
  library(ggplot2)
})

set.seed(20260825)

ROOT <- "E:/sepsis_project"
RESULTS <- file.path(ROOT, "results")
IN_DIR <- file.path(
  RESULTS,
  "V2_36F_PRJNA851469_CLINICAL_RECOVERY_LINKAGE"
)
OUT <- file.path(
  RESULTS,
  "V2_36G_PRJNA851469_DAY3_CLINICAL_OUTCOME_MODELS"
)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

IN_FILE <- file.path(IN_DIR, "03_day3_landmark_candidate.csv")
if (!file.exists(IN_FILE)) stop("Step36F landmark file not found: ", IN_FILE)

d <- read_csv(IN_FILE, show_col_types=FALSE) %>%
  filter(landmark_eligible) %>%
  mutate(
    event = as.integer(subsequent_composite_event),
    abx = as.integer(antibiotics_at_icu_admission),
    sepsis = as.integer(sepsis_admission),
    z_bray = as.numeric(z_day3_bray),
    bray = as.numeric(bray_from_patient_baseline),
    followup = as.numeric(followup_days_from_day3)
  )

req <- c("event","abx","sepsis","z_bray","bray","followup")
if (any(!complete.cases(d[,req]))) {
  stop("Unexpected missingness in required Step36G analysis variables.")
}

if (nrow(d) < 25 || sum(d$event) < 8) {
  stop("Landmark dataset no longer meets prespecified readiness threshold.")
}

write_csv(d, file.path(OUT, "01_day3_landmark_analysis_population.csv"))

# ------------------------------------------------------------
# 1. Descriptive comparison
# ------------------------------------------------------------
desc <- d %>%
  group_by(event) %>%
  summarise(
    n=n(),
    bray_mean=mean(bray),
    bray_sd=sd(bray),
    bray_median=median(bray),
    bray_q1=quantile(bray,.25),
    bray_q3=quantile(bray,.75),
    antibiotic_yes=sum(abx==1),
    sepsis_admission=sum(sepsis==1),
    .groups="drop"
  )
write_csv(desc, file.path(OUT, "02_descriptive_by_subsequent_event.csv"))

mw <- wilcox.test(bray ~ event, data=d, exact=FALSE)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
tidy_cox <- function(model, model_name) {
  s <- summary(model)
  cf <- as.data.frame(s$coefficients)
  ci <- as.data.frame(s$conf.int)
  tibble(
    model=model_name,
    term=rownames(cf),
    beta=cf$coef,
    HR=ci$`exp(coef)`,
    CI_low=ci$`lower .95`,
    CI_high=ci$`upper .95`,
    p=cf$`Pr(>|z|)`
  )
}

tidy_logit <- function(model, model_name) {
  sm <- summary(model)$coefficients
  tibble(
    model=model_name,
    term=rownames(sm),
    beta=sm[,1],
    SE=sm[,2],
    OR=exp(sm[,1]),
    CI_low=exp(sm[,1]-1.96*sm[,2]),
    CI_high=exp(sm[,1]+1.96*sm[,2]),
    p=sm[,4]
  )
}

# ------------------------------------------------------------
# 3. Cox landmark models
# ------------------------------------------------------------
surv_obj <- with(d, Surv(followup, event))

cox1 <- coxph(surv_obj ~ z_bray, data=d, ties="efron", x=TRUE)
cox2 <- coxph(surv_obj ~ z_bray + abx, data=d, ties="efron", x=TRUE)
cox3 <- coxph(surv_obj ~ z_bray + abx + sepsis, data=d, ties="efron", x=TRUE)

cox_results <- bind_rows(
  tidy_cox(cox1, "M1_PRIMARY_zBray"),
  tidy_cox(cox2, "M2_PLUS_ADMISSION_ANTIBIOTICS"),
  tidy_cox(cox3, "M3_PLUS_ANTIBIOTICS_AND_SEPSIS_ADMISSION")
)
write_csv(cox_results, file.path(OUT, "03_cox_landmark_models.csv"))

# proportional hazards diagnostics
zph1 <- cox.zph(cox1)
zph2 <- cox.zph(cox2)
zph3 <- cox.zph(cox3)

zph_to_df <- function(z, name) {
  x <- as.data.frame(z$table)
  tibble(
    model=name,
    term=rownames(x),
    chisq=x[,1],
    p=x[,3]
  )
}

ph <- bind_rows(
  zph_to_df(zph1,"M1_PRIMARY_zBray"),
  zph_to_df(zph2,"M2_PLUS_ADMISSION_ANTIBIOTICS"),
  zph_to_df(zph3,"M3_PLUS_ANTIBIOTICS_AND_SEPSIS_ADMISSION")
)
write_csv(ph, file.path(OUT, "04_cox_PH_diagnostics.csv"))

# ------------------------------------------------------------
# 4. Logistic horizon sensitivities
# ------------------------------------------------------------
log1 <- glm(event ~ z_bray, data=d, family=binomial())
log2 <- glm(event ~ z_bray + abx, data=d, family=binomial())
log3 <- glm(event ~ z_bray + abx + sepsis, data=d, family=binomial())

log_results <- bind_rows(
  tidy_logit(log1,"L1_PRIMARY_zBray"),
  tidy_logit(log2,"L2_PLUS_ADMISSION_ANTIBIOTICS"),
  tidy_logit(log3,"L3_PLUS_ANTIBIOTICS_AND_SEPSIS_ADMISSION")
)
write_csv(log_results, file.path(OUT, "05_logistic_horizon_models.csv"))

# ------------------------------------------------------------
# 5. Bootstrap robustness
# Patient-level resampling; primary Cox HR for z_bray
# ------------------------------------------------------------
B <- 3000
boot_beta <- rep(NA_real_, B)

for (b in seq_len(B)) {
  idx <- sample(seq_len(nrow(d)), replace=TRUE)
  db <- d[idx, , drop=FALSE]

  fit <- try(
    coxph(Surv(followup,event) ~ z_bray, data=db, ties="efron"),
    silent=TRUE
  )

  if (!inherits(fit,"try-error")) {
    cc <- coef(fit)
    if ("z_bray" %in% names(cc) && is.finite(cc["z_bray"])) {
      boot_beta[b] <- cc["z_bray"]
    }
  }
}

boot_beta <- boot_beta[is.finite(boot_beta)]
boot_summary <- tibble(
  B_requested=B,
  B_valid=length(boot_beta),
  beta_median=median(boot_beta),
  HR_median=exp(median(boot_beta)),
  HR_CI_low=exp(quantile(boot_beta,.025)),
  HR_CI_high=exp(quantile(boot_beta,.975)),
  proportion_HR_gt_1=mean(boot_beta > 0)
)
write_csv(boot_summary, file.path(OUT, "06_bootstrap_primary_cox.csv"))

# ------------------------------------------------------------
# 6. Permutation robustness
# Shuffle ecological displacement across patients while preserving
# event times/censoring. Statistic = absolute Cox z-statistic.
# ------------------------------------------------------------
obs_z <- summary(cox1)$coefficients["z_bray","z"]
P <- 5000
perm_z <- rep(NA_real_, P)

for (i in seq_len(P)) {
  dp <- d
  dp$z_bray <- sample(dp$z_bray, replace=FALSE)

  fit <- try(
    coxph(Surv(followup,event) ~ z_bray, data=dp, ties="efron"),
    silent=TRUE
  )

  if (!inherits(fit,"try-error")) {
    zz <- summary(fit)$coefficients["z_bray","z"]
    if (is.finite(zz)) perm_z[i] <- zz
  }
}

perm_z <- perm_z[is.finite(perm_z)]
perm_p <- (1 + sum(abs(perm_z) >= abs(obs_z))) / (1 + length(perm_z))

perm_summary <- tibble(
  permutations_requested=P,
  permutations_valid=length(perm_z),
  observed_z=obs_z,
  permutation_p_two_sided=perm_p
)
write_csv(perm_summary, file.path(OUT, "07_permutation_primary_cox.csv"))

# ------------------------------------------------------------
# 7. Prespecified evidence gate
# ------------------------------------------------------------
get_term <- function(tbl, model, term) {
  tbl %>% filter(.data$model == model, .data$term == term)
}

p1 <- get_term(cox_results,"M1_PRIMARY_zBray","z_bray")
p2 <- get_term(cox_results,"M2_PLUS_ADMISSION_ANTIBIOTICS","z_bray")
p3 <- get_term(cox_results,"M3_PLUS_ANTIBIOTICS_AND_SEPSIS_ADMISSION","z_bray")

if (nrow(p1)!=1 || nrow(p2)!=1 || nrow(p3)!=1) {
  stop("Could not retrieve z_bray model terms.")
}

same_direction <- sign(p1$beta) == sign(p2$beta) &&
                  sign(p1$beta) == sign(p3$beta)

primary_significant <- p1$p < 0.05
adjusted_supported <- p2$p < 0.05
bootstrap_excludes_null <- boot_summary$HR_CI_low > 1 ||
                           boot_summary$HR_CI_high < 1
perm_supported <- perm_p < 0.05

tier <- case_when(
  primary_significant &&
    adjusted_supported &&
    same_direction &&
    bootstrap_excludes_null &&
    perm_supported ~
      "SUPPORTIVE_CLINICAL_OUTCOME_ASSOCIATION_ROBUST",

  primary_significant &&
    same_direction ~
      "EXPLORATORY_CLINICAL_OUTCOME_ASSOCIATION",

  TRUE ~
      "NO_CLEAR_CLINICAL_OUTCOME_ASSOCIATION"
)

decision <- tibble(
  n_landmark=nrow(d),
  n_events=sum(d$event),
  primary_HR=p1$HR,
  primary_CI_low=p1$CI_low,
  primary_CI_high=p1$CI_high,
  primary_p=p1$p,
  antibiotic_adjusted_HR=p2$HR,
  antibiotic_adjusted_CI_low=p2$CI_low,
  antibiotic_adjusted_CI_high=p2$CI_high,
  antibiotic_adjusted_p=p2$p,
  fully_adjusted_HR=p3$HR,
  fully_adjusted_CI_low=p3$CI_low,
  fully_adjusted_CI_high=p3$CI_high,
  fully_adjusted_p=p3$p,
  bootstrap_HR_low=boot_summary$HR_CI_low,
  bootstrap_HR_high=boot_summary$HR_CI_high,
  permutation_p=perm_p,
  same_direction_across_models=same_direction,
  evidence_tier=tier
)

write_csv(decision, file.path(OUT, "08_EVIDENCE_DECISION.csv"))

# ------------------------------------------------------------
# 8. Figure: continuous displacement, event status
# Descriptive only; inference comes from Cox models.
# ------------------------------------------------------------
plot_dat <- d %>%
  mutate(
    event_label=ifelse(event==1,
      "Subsequent infection/death",
      "No event through day 30")
  )

p <- ggplot(plot_dat, aes(x=event_label, y=bray)) +
  geom_boxplot(outlier.shape=NA, width=.5) +
  geom_jitter(width=.10, height=0, size=2, alpha=.8) +
  labs(
    x=NULL,
    y="Day-3 Bray displacement from patient baseline",
    title="PRJNA851469 Day-3 landmark analysis"
  ) +
  theme_classic(base_size=11)

ggsave(
  file.path(OUT,"09_Figure_candidate_day3_displacement_event.pdf"),
  p, width=6.2, height=4.6
)
ggsave(
  file.path(OUT,"09_Figure_candidate_day3_displacement_event.png"),
  p, width=6.2, height=4.6, dpi=300
)

# ------------------------------------------------------------
# 9. Manuscript-safe summary
# ------------------------------------------------------------
summary_lines <- c(
  "V2 STEP36G PRJNA851469 DAY-3 LANDMARK CLINICAL OUTCOME ANALYSIS",
  "",
  paste0("Landmark population: n=", nrow(d),
         "; subsequent infection/death events=", sum(d$event), "."),
  paste0(
    "Primary Cox model, per 1-SD greater Day-3 Bray displacement: HR=",
    sprintf("%.3f",p1$HR),
    " (95% CI ",sprintf("%.3f",p1$CI_low),
    " to ",sprintf("%.3f",p1$CI_high),
    "), p=",sprintf("%.4g",p1$p), "."
  ),
  paste0(
    "Adjusted for antibiotics at ICU admission: HR=",
    sprintf("%.3f",p2$HR),
    " (95% CI ",sprintf("%.3f",p2$CI_low),
    " to ",sprintf("%.3f",p2$CI_high),
    "), p=",sprintf("%.4g",p2$p), "."
  ),
  paste0(
    "Further sensitivity including sepsis admission: HR=",
    sprintf("%.3f",p3$HR),
    " (95% CI ",sprintf("%.3f",p3$CI_low),
    " to ",sprintf("%.3f",p3$CI_high),
    "), p=",sprintf("%.4g",p3$p), "."
  ),
  paste0(
    "Bootstrap 95% HR interval: ",
    sprintf("%.3f",boot_summary$HR_CI_low),
    " to ",sprintf("%.3f",boot_summary$HR_CI_high), "."
  ),
  paste0("Permutation p=",sprintf("%.4g",perm_p), "."),
  "",
  paste0("EVIDENCE TIER: ",tier),
  "",
  "Interpretation guardrail:",
  "This is an exploratory external clinical-context analysis in one ICU-background cohort.",
  "Do not describe the ecological displacement measure as a validated prognostic biomarker.",
  "Admission-antibiotic status is a binary published-timeline variable and must not be interpreted causally.",
  "No patient-level SOFA adjustment was possible from public data."
)

writeLines(summary_lines, file.path(OUT, "10_MANUSCRIPT_SAFE_SUMMARY.txt"))

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP36G COMPLETE",
    paste0("Evidence tier: ",tier),
    "No previously frozen V2 result was modified."
  ),
  file.path(OUT, "_STEP36G_COMPLETE.txt")
)

cat("\nSTEP36G COMPLETE\n")
cat("Evidence tier:",tier,"\n")
