
# ============================================================
# V2_36G2 PRJNA851469 DAY-3 LANDMARK CLINICAL OUTCOME MODELS
# FIXED VERSION
#
# Fix:
#   - Corrects dplyr data-mask bug in evidence gate.
#   - Uses explicit .env / renamed arguments.
#   - Preserves all prespecified analyses from Step36G.
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
  "V2_36G2_PRJNA851469_DAY3_CLINICAL_OUTCOME_MODELS_FIXED"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

IN_FILE <- file.path(IN_DIR, "03_day3_landmark_candidate.csv")

if (!file.exists(IN_FILE)) {
  stop("Step36F landmark file not found: ", IN_FILE)
}

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
  stop("Unexpected missingness in required Step36G2 variables.")
}

if (nrow(d) < 25 || sum(d$event) < 8) {
  stop("Landmark dataset does not meet prespecified threshold.")
}

write_csv(
  d,
  file.path(OUT, "01_day3_landmark_analysis_population.csv")
)

# ------------------------------------------------------------
# Descriptive
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

write_csv(
  desc,
  file.path(OUT, "02_descriptive_by_subsequent_event.csv")
)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
tidy_cox <- function(fit, fit_name) {
  s <- summary(fit)
  cf <- as.data.frame(s$coefficients)
  ci <- as.data.frame(s$conf.int)

  tibble(
    model=fit_name,
    term=rownames(cf),
    beta=cf$coef,
    HR=ci$`exp(coef)`,
    CI_low=ci$`lower .95`,
    CI_high=ci$`upper .95`,
    p=cf$`Pr(>|z|)`
  )
}

tidy_logit <- function(fit, fit_name) {
  sm <- summary(fit)$coefficients

  tibble(
    model=fit_name,
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
# Cox models
# ------------------------------------------------------------
cox1 <- coxph(
  Surv(followup,event) ~ z_bray,
  data=d,
  ties="efron",
  x=TRUE
)

cox2 <- coxph(
  Surv(followup,event) ~ z_bray + abx,
  data=d,
  ties="efron",
  x=TRUE
)

cox3 <- coxph(
  Surv(followup,event) ~ z_bray + abx + sepsis,
  data=d,
  ties="efron",
  x=TRUE
)

cox_results <- bind_rows(
  tidy_cox(cox1,"M1_PRIMARY_zBray"),
  tidy_cox(cox2,"M2_PLUS_ADMISSION_ANTIBIOTICS"),
  tidy_cox(cox3,"M3_PLUS_ANTIBIOTICS_AND_SEPSIS_ADMISSION")
)

write_csv(
  cox_results,
  file.path(OUT, "03_cox_landmark_models.csv")
)

# ------------------------------------------------------------
# PH diagnostics
# ------------------------------------------------------------
zph_to_df <- function(z, fit_name) {
  x <- as.data.frame(z$table)
  tibble(
    model=fit_name,
    term=rownames(x),
    chisq=x[,1],
    p=x[,3]
  )
}

ph <- bind_rows(
  zph_to_df(cox.zph(cox1),"M1_PRIMARY_zBray"),
  zph_to_df(cox.zph(cox2),"M2_PLUS_ADMISSION_ANTIBIOTICS"),
  zph_to_df(cox.zph(cox3),"M3_PLUS_ANTIBIOTICS_AND_SEPSIS_ADMISSION")
)

write_csv(
  ph,
  file.path(OUT, "04_cox_PH_diagnostics.csv")
)

# ------------------------------------------------------------
# Logistic horizon sensitivity
# ------------------------------------------------------------
log1 <- glm(event ~ z_bray, data=d, family=binomial())
log2 <- glm(event ~ z_bray + abx, data=d, family=binomial())
log3 <- glm(event ~ z_bray + abx + sepsis, data=d, family=binomial())

log_results <- bind_rows(
  tidy_logit(log1,"L1_PRIMARY_zBray"),
  tidy_logit(log2,"L2_PLUS_ADMISSION_ANTIBIOTICS"),
  tidy_logit(log3,"L3_PLUS_ANTIBIOTICS_AND_SEPSIS_ADMISSION")
)

write_csv(
  log_results,
  file.path(OUT, "05_logistic_horizon_models.csv")
)

# ------------------------------------------------------------
# Bootstrap
# ------------------------------------------------------------
B <- 3000
boot_beta <- rep(NA_real_, B)

for (b in seq_len(B)) {

  idx <- sample(seq_len(nrow(d)), replace=TRUE)
  db <- d[idx,,drop=FALSE]

  fit <- try(
    coxph(
      Surv(followup,event) ~ z_bray,
      data=db,
      ties="efron"
    ),
    silent=TRUE
  )

  if (!inherits(fit,"try-error")) {
    cc <- coef(fit)

    if (
      "z_bray" %in% names(cc) &&
      is.finite(cc["z_bray"])
    ) {
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
  HR_CI_low=exp(unname(quantile(boot_beta,.025))),
  HR_CI_high=exp(unname(quantile(boot_beta,.975))),
  proportion_HR_gt_1=mean(boot_beta > 0)
)

write_csv(
  boot_summary,
  file.path(OUT, "06_bootstrap_primary_cox.csv")
)

# ------------------------------------------------------------
# Permutation
# ------------------------------------------------------------
obs_z <- summary(cox1)$coefficients["z_bray","z"]

P <- 5000
perm_z <- rep(NA_real_, P)

for (i in seq_len(P)) {

  dp <- d
  dp$z_bray <- sample(dp$z_bray, replace=FALSE)

  fit <- try(
    coxph(
      Surv(followup,event) ~ z_bray,
      data=dp,
      ties="efron"
    ),
    silent=TRUE
  )

  if (!inherits(fit,"try-error")) {
    zz <- summary(fit)$coefficients["z_bray","z"]
    if (is.finite(zz)) perm_z[i] <- zz
  }
}

perm_z <- perm_z[is.finite(perm_z)]

perm_p <- (
  1 + sum(abs(perm_z) >= abs(obs_z))
) / (
  1 + length(perm_z)
)

perm_summary <- tibble(
  permutations_requested=P,
  permutations_valid=length(perm_z),
  observed_z=obs_z,
  permutation_p_two_sided=perm_p
)

write_csv(
  perm_summary,
  file.path(OUT, "07_permutation_primary_cox.csv")
)

# ------------------------------------------------------------
# FIXED evidence gate
# ------------------------------------------------------------
get_model_term <- function(tbl, model_name, term_name) {

  out <- tbl %>%
    filter(
      .data$model == .env$model_name,
      .data$term == .env$term_name
    )

  if (nrow(out) != 1) {
    stop(
      "Evidence gate expected exactly one row for ",
      model_name, " / ", term_name,
      "; found ", nrow(out)
    )
  }

  out
}

p1 <- get_model_term(
  cox_results,
  "M1_PRIMARY_zBray",
  "z_bray"
)

p2 <- get_model_term(
  cox_results,
  "M2_PLUS_ADMISSION_ANTIBIOTICS",
  "z_bray"
)

p3 <- get_model_term(
  cox_results,
  "M3_PLUS_ANTIBIOTICS_AND_SEPSIS_ADMISSION",
  "z_bray"
)

same_direction <- (
  sign(p1$beta[[1]]) == sign(p2$beta[[1]]) &&
  sign(p1$beta[[1]]) == sign(p3$beta[[1]])
)

primary_significant <- p1$p[[1]] < 0.05
adjusted_supported <- p2$p[[1]] < 0.05

bootstrap_excludes_null <- (
  boot_summary$HR_CI_low[[1]] > 1 ||
  boot_summary$HR_CI_high[[1]] < 1
)

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

  primary_HR=p1$HR[[1]],
  primary_CI_low=p1$CI_low[[1]],
  primary_CI_high=p1$CI_high[[1]],
  primary_p=p1$p[[1]],

  antibiotic_adjusted_HR=p2$HR[[1]],
  antibiotic_adjusted_CI_low=p2$CI_low[[1]],
  antibiotic_adjusted_CI_high=p2$CI_high[[1]],
  antibiotic_adjusted_p=p2$p[[1]],

  fully_adjusted_HR=p3$HR[[1]],
  fully_adjusted_CI_low=p3$CI_low[[1]],
  fully_adjusted_CI_high=p3$CI_high[[1]],
  fully_adjusted_p=p3$p[[1]],

  bootstrap_HR_low=boot_summary$HR_CI_low[[1]],
  bootstrap_HR_high=boot_summary$HR_CI_high[[1]],

  permutation_p=perm_p,
  same_direction_across_models=same_direction,
  evidence_tier=tier
)

write_csv(
  decision,
  file.path(OUT, "08_EVIDENCE_DECISION.csv")
)

# ------------------------------------------------------------
# Figure
# ------------------------------------------------------------
plot_dat <- d %>%
  mutate(
    event_label=ifelse(
      event==1,
      "Subsequent infection/death",
      "No event through day 30"
    )
  )

p <- ggplot(
  plot_dat,
  aes(x=event_label,y=bray)
) +
  geom_boxplot(
    outlier.shape=NA,
    width=.5
  ) +
  geom_jitter(
    width=.10,
    height=0,
    size=2,
    alpha=.8
  ) +
  labs(
    x=NULL,
    y="Day-3 Bray displacement from patient baseline",
    title="PRJNA851469 Day-3 landmark analysis"
  ) +
  theme_classic(base_size=11)

ggsave(
  file.path(
    OUT,
    "09_Figure_candidate_day3_displacement_event.pdf"
  ),
  p,
  width=6.2,
  height=4.6
)

ggsave(
  file.path(
    OUT,
    "09_Figure_candidate_day3_displacement_event.png"
  ),
  p,
  width=6.2,
  height=4.6,
  dpi=300
)

# ------------------------------------------------------------
# Manuscript-safe summary
# ------------------------------------------------------------
summary_lines <- c(

  "V2 STEP36G2 PRJNA851469 DAY-3 LANDMARK CLINICAL OUTCOME ANALYSIS",

  "",

  paste0(
    "Landmark population: n=",
    nrow(d),
    "; subsequent infection/death events=",
    sum(d$event),
    "."
  ),

  paste0(
    "Primary Cox model, per 1-SD greater Day-3 Bray displacement: HR=",
    sprintf("%.3f",p1$HR[[1]]),
    " (95% CI ",
    sprintf("%.3f",p1$CI_low[[1]]),
    " to ",
    sprintf("%.3f",p1$CI_high[[1]]),
    "), p=",
    sprintf("%.4g",p1$p[[1]]),
    "."
  ),

  paste0(
    "Adjusted for antibiotics at ICU admission: HR=",
    sprintf("%.3f",p2$HR[[1]]),
    " (95% CI ",
    sprintf("%.3f",p2$CI_low[[1]]),
    " to ",
    sprintf("%.3f",p2$CI_high[[1]]),
    "), p=",
    sprintf("%.4g",p2$p[[1]]),
    "."
  ),

  paste0(
    "Further sensitivity including sepsis admission: HR=",
    sprintf("%.3f",p3$HR[[1]]),
    " (95% CI ",
    sprintf("%.3f",p3$CI_low[[1]]),
    " to ",
    sprintf("%.3f",p3$CI_high[[1]]),
    "), p=",
    sprintf("%.4g",p3$p[[1]]),
    "."
  ),

  paste0(
    "Bootstrap 95% HR interval: ",
    sprintf("%.3f",boot_summary$HR_CI_low[[1]]),
    " to ",
    sprintf("%.3f",boot_summary$HR_CI_high[[1]]),
    "."
  ),

  paste0(
    "Permutation p=",
    sprintf("%.4g",perm_p),
    "."
  ),

  "",

  paste0(
    "EVIDENCE TIER: ",
    tier
  ),

  "",

  "Interpretation guardrails:",

  "1. Exploratory/supportive clinical-context analysis in one ICU-background cohort.",
  "2. Do not describe ecological displacement as a validated prognostic biomarker.",
  "3. Admission-antibiotic status is binary and must not be interpreted causally.",
  "4. No individual-patient SOFA adjustment was possible from public data.",
  "5. No previously frozen V2 inference was changed."
)

writeLines(
  summary_lines,
  file.path(
    OUT,
    "10_MANUSCRIPT_SAFE_SUMMARY.txt"
  )
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP36G2 COMPLETE",
    paste0("Evidence tier: ",tier),
    "No previously frozen V2 result was modified."
  ),
  file.path(
    OUT,
    "_STEP36G2_COMPLETE.txt"
  )
)

cat("\nSTEP36G2 COMPLETE\n")
cat("Evidence tier:",tier,"\n")
