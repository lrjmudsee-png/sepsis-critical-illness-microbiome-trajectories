# ============================================================
# Sepsis V2 - Step 81D
# FINAL clinical repair for CRA002354 + re-freeze
#
# Why this step exists:
# Step81C correctly fixed sequence-layer modality and produced
# the correct 18-project / 1941-row primary 16S freeze.
# However, CRA002354 patient-level clinical fields were blank in
# the frozen master (outcome_28d, baseline/sample sepsis status,
# age/sex, SOFA/APACHE, antibiotics), although the supplement
# S1/S2 files contain them.
#
# Step81D repairs CRA002354 directly from the validated original
# supplementary workbooks and regenerates the final summaries.
#
# OFFLINE ONLY
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c("readr","readxl","dplyr","tidyr","stringr","purrr","tibble")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) install.packages(missing, repos="https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"

IN <- file.path(
  ROOT,
  "results",
  "V2_21C_FINAL_COHORT_FREEZE_CORRECTED"
)

CRA_SUPP <- file.path(
  ROOT,
  "data",
  "_step78A3_US_fast_download",
  "CRA002354"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_21D_FINAL_FREEZE_CLINICAL_REPAIRED"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

LOG <- file.path(OUT, "_STEP81D_runtime_checkpoints.txt")
ERR <- file.path(OUT, "_STEP81D_FATAL_ERROR.txt")
if (file.exists(ERR)) unlink(ERR)

ck <- function(x) {
  cat(paste0(x, ": ", Sys.time(), "\n"), file=LOG, append=TRUE)
}

clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[
    is.na(x) | x=="" |
    tolower(x) %in% c("na","nan","n/a","null","none")
  ] <- NA_character_
  x
}

clean_num <- function(x) suppressWarnings(as.numeric(clean_chr(x)))

safe_csv <- function(p) {
  suppressMessages(
    read_csv(
      p,
      show_col_types=FALSE,
      progress=FALSE,
      name_repair="unique"
    )
  )
}

find_latest <- function(folder, pattern) {
  x <- list.files(folder, pattern=pattern, full.names=TRUE)
  if (!length(x)) return(NA_character_)
  x[which.max(file.info(x)$mtime)]
}

first_existing <- function(paths) {
  z <- paths[file.exists(paths)]
  if (!length(z)) return(NA_character_)
  z[1]
}

find_regex_col <- function(nms, patterns) {
  pat <- paste(patterns, collapse="|")
  idx <- which(str_detect(nms, regex(pat, ignore_case=TRUE)))
  if (!length(idx)) return(NA_character_)
  nms[idx[1]]
}

read_source_value <- function(df, col) {
  if (is.na(col) || !(col %in% names(df))) {
    rep(NA_character_, nrow(df))
  } else {
    clean_chr(df[[col]])
  }
}

combine_severity <- function(sofa, apache) {
  sofa <- clean_chr(sofa)
  apache <- clean_chr(apache)

  out <- rep(NA_character_, length(sofa))

  for (i in seq_along(out)) {
    z <- character()

    if (!is.na(sofa[i])) {
      z <- c(z, paste0("SOFA=", sofa[i]))
    }

    if (!is.na(apache[i])) {
      z <- c(z, paste0("APACHEII=", apache[i]))
    }

    if (length(z)) out[i] <- paste(z, collapse=";")
  }

  out
}

main <- function() {

  ck("STEP81D STARTED")

  # ==========================================================
  # 1. Read Step81C master
  # ==========================================================

  master_path <- find_latest(
    IN,
    "^V2_FINAL_FROZEN_master_metadata_[0-9]{8}_[0-9]{6}\\.csv$"
  )

  if (is.na(master_path) || !file.exists(master_path)) {
    stop("Step81C final frozen master was not found.")
  }

  full <- safe_csv(master_path)

  if (nrow(full) != 2693 || dplyr::n_distinct(full$project) != 18) {
    stop("Step81C master guard failed: expected 2693 rows / 18 projects.")
  }

  if (sum(full$primary_16s_sequence_include, na.rm=TRUE) != 1941) {
    stop("Step81C master guard failed: expected 1941 primary 16S rows.")
  }

  ck("STEP81C MASTER READ AND GUARDED")

  # ==========================================================
  # 2. Locate validated CRA supplementary tables
  # ==========================================================

  s1_path <- first_existing(c(
    file.path(CRA_SUPP, "CRA002354_Table_S1_64_patients.xlsx"),
    file.path(CRA_SUPP, "mmc5.xlsx")
  ))

  s2_path <- first_existing(c(
    file.path(CRA_SUPP, "CRA002354_Table_S2_131_samples.xlsx"),
    file.path(CRA_SUPP, "mmc6.xlsx")
  ))

  if (is.na(s1_path) || is.na(s2_path)) {
    stop("CRA002354 S1/S2 supplementary XLSX files were not found.")
  }

  # True headers are Excel row 3.
  s1 <- read_excel(
    s1_path,
    skip=2,
    .name_repair="unique"
  ) |>
    mutate(across(everything(), clean_chr))

  s2 <- read_excel(
    s2_path,
    skip=2,
    .name_repair="unique"
  ) |>
    mutate(across(everything(), clean_chr))

  if (nrow(s1) != 64) {
    stop(paste0("CRA S1 expected 64 patient rows; observed ", nrow(s1)))
  }

  if (nrow(s2) != 131) {
    stop(paste0("CRA S2 expected 131 sample rows; observed ", nrow(s2)))
  }

  ck("CRA S1/S2 READ")

  # ==========================================================
  # 3. Identify source fields
  # ==========================================================

  s1_pid <- find_regex_col(names(s1), c("^Patient ID$"))
  s1_sid <- find_regex_col(names(s1), c("^Sample ID$"))

  s2_pid <- find_regex_col(names(s2), c("^Patient ID$"))
  s2_sid <- find_regex_col(names(s2), c("^Sample ID$"))

  s1_status <- find_regex_col(
    names(s1),
    c(
      "Sepsis/Septic shock",
      "sepsis.*septic shock"
    )
  )

  s2_status <- find_regex_col(
    names(s2),
    c(
      "^Sepsis\\s*=\\s*1.*septic shock\\s*=\\s*2",
      "sepsis.*septic shock"
    )
  )

  s1_survival <- find_regex_col(
    names(s1),
    c(
      "28.*day.*survival",
      "survival.*28"
    )
  )

  s1_infection <- find_regex_col(
    names(s1),
    c("^Site of infection")
  )

  s1_gender <- find_regex_col(
    names(s1),
    c("^Gender$", "^Sex$")
  )

  s1_age <- find_regex_col(
    names(s1),
    c("^Age")
  )

  s1_sofa <- find_regex_col(
    names(s1),
    c("^SOFA", "SOFA.*Sequential")
  )

  s1_apache <- find_regex_col(
    names(s1),
    c("APACHE")
  )

  s1_lactate <- find_regex_col(
    names(s1),
    c("lactate")
  )

  s1_antibiotics <- find_regex_col(
    names(s1),
    c(
      "Antibiotics use.*ICU",
      "Antibiotic"
    )
  )

  required_cols <- c(
    s1_pid, s1_sid, s2_pid, s2_sid,
    s1_status, s2_status, s1_survival, s1_infection
  )

  if (any(is.na(required_cols))) {
    stop(
      paste0(
        "Required CRA clinical columns were not all found.",
        "\nS1 columns: ", paste(names(s1), collapse=" | "),
        "\nS2 columns: ", paste(names(s2), collapse=" | ")
      )
    )
  }

  field_audit <- tibble(
    standardized_field=c(
      "patient_id",
      "first_sample_id",
      "baseline_sepsis_status",
      "outcome_28d",
      "infection_source",
      "sex_raw",
      "age",
      "baseline_sofa",
      "baseline_apache_ii",
      "baseline_lactate",
      "antibiotics",
      "sample_id",
      "sample_sepsis_status"
    ),

    source_table=c(
      rep("S1_patient_level",11),
      rep("S2_sample_level",2)
    ),

    source_column=c(
      s1_pid,
      s1_sid,
      s1_status,
      s1_survival,
      s1_infection,
      s1_gender,
      s1_age,
      s1_sofa,
      s1_apache,
      s1_lactate,
      s1_antibiotics,
      s2_sid,
      s2_status
    )
  )

  write_excel_csv(
    field_audit,
    file.path(OUT, "CRA002354_clinical_source_field_audit.csv"),
    na=""
  )

  # ==========================================================
  # 4. Patient-level map
  # ==========================================================

  infection_codes <- c(
    "0"="other",
    "1"="lung",
    "2"="intestinal",
    "3"="abdominal",
    "4"="blood",
    "5"="urinary",
    "6"="brain",
    "7"="surgical_site"
  )

  patient_map <- tibble(
    patient_id=clean_chr(s1[[s1_pid]]),

    baseline_sepsis_status=case_when(
      clean_chr(s1[[s1_status]])=="1" ~ "sepsis",
      clean_chr(s1[[s1_status]])=="2" ~ "septic_shock",
      TRUE ~ NA_character_
    ),

    outcome_28d=case_when(
      clean_chr(s1[[s1_survival]])=="1" ~ "survived",
      clean_chr(s1[[s1_survival]])=="2" ~ "dead",
      TRUE ~ NA_character_
    ),

    infection_source_standard=unname(
      infection_codes[
        clean_chr(s1[[s1_infection]])
      ]
    ),

    sex_raw=read_source_value(s1, s1_gender),

    age_repaired=if (!is.na(s1_age)) {
      clean_num(s1[[s1_age]])
    } else {
      rep(NA_real_, nrow(s1))
    },

    baseline_sofa=if (!is.na(s1_sofa)) {
      clean_num(s1[[s1_sofa]])
    } else {
      rep(NA_real_, nrow(s1))
    },

    baseline_apache_ii=if (!is.na(s1_apache)) {
      clean_num(s1[[s1_apache]])
    } else {
      rep(NA_real_, nrow(s1))
    },

    baseline_lactate=if (!is.na(s1_lactate)) {
      clean_num(s1[[s1_lactate]])
    } else {
      rep(NA_real_, nrow(s1))
    },

    antibiotics_repaired=read_source_value(
      s1,
      s1_antibiotics
    )
  ) |>
    distinct()

  if (
    nrow(patient_map) != 64 ||
    n_distinct(patient_map$patient_id) != 64
  ) {
    stop("CRA patient map is not 64 unique patients.")
  }

  # Hard clinical guards from the validated supplement.
  if (
    sum(patient_map$baseline_sepsis_status=="sepsis",na.rm=TRUE) != 46 ||
    sum(patient_map$baseline_sepsis_status=="septic_shock",na.rm=TRUE) != 18
  ) {
    stop("CRA baseline sepsis/shock guard failed: expected 46/18.")
  }

  if (
    sum(patient_map$outcome_28d=="survived",na.rm=TRUE) != 44 ||
    sum(patient_map$outcome_28d=="dead",na.rm=TRUE) != 20
  ) {
    stop("CRA 28-day outcome guard failed: expected 44 survived / 20 dead.")
  }

  # ==========================================================
  # 5. Sample-level status map
  # ==========================================================

  sample_map <- tibble(
    sample_id=clean_chr(s2[[s2_sid]]),

    source_patient_id=clean_chr(s2[[s2_pid]]),

    sample_sepsis_status=case_when(
      clean_chr(s2[[s2_status]])=="1" ~ "sepsis",
      clean_chr(s2[[s2_status]])=="2" ~ "septic_shock",
      TRUE ~ NA_character_
    )
  ) |>
    distinct()

  if (
    nrow(sample_map) != 131 ||
    n_distinct(sample_map$sample_id) != 131
  ) {
    stop("CRA sample map is not 131 unique samples.")
  }

  # ==========================================================
  # 6. Repair CRA rows in frozen master
  # ==========================================================

  cra_old <- full |>
    filter(project=="CRA002354")

  other <- full |>
    filter(project!="CRA002354")

  if (
    nrow(cra_old) != 131 ||
    n_distinct(cra_old$run_id, na.rm=TRUE) != 131
  ) {
    stop("Frozen CRA rows are not the expected 131/131 exact technical map.")
  }

  # Preserve existing exact technical identity fields.
  old_run_signature <- cra_old |>
    select(sample_id, experiment_id, run_id) |>
    arrange(sample_id)

  cra_new <- cra_old |>
    select(
      -any_of(c(
        "baseline_sepsis_status",
        "sample_sepsis_status",
        "outcome_28d",
        "infection_source_standard",
        "outcome",
        "infection_group",
        "sex",
        "age",
        "severity",
        "antibiotics",
        "baseline_sofa",
        "baseline_apache_ii",
        "baseline_lactate",
        "sex_raw"
      ))
    ) |>
    left_join(
      patient_map,
      by="patient_id"
    ) |>
    left_join(
      sample_map |>
        select(
          sample_id,
          source_patient_id,
          sample_sepsis_status
        ),
      by="sample_id"
    ) |>
    mutate(
      # Explicitly enforce patient/sample linkage agreement.
      patient_sample_link_match=
        patient_id==source_patient_id,

      outcome=outcome_28d,

      infection_group=
        infection_source_standard,

      # Keep raw source sex value without inventing code meaning.
      sex=sex_raw,

      age=age_repaired,

      severity=combine_severity(
        baseline_sofa,
        baseline_apache_ii
      ),

      antibiotics=antibiotics_repaired,

      clinical_metadata_source=
        "CRA002354 Supplementary Table S1 patient-level + S2 sample-level",

      clinical_mapping_status=
        "EXACT_PATIENT_ID_AND_SAMPLE_ID_JOIN"
    ) |>
    select(
      -source_patient_id,
      -age_repaired,
      -antibiotics_repaired
    )

  if (!all(cra_new$patient_sample_link_match, na.rm=TRUE)) {
    stop("CRA S1/S2 patient-sample linkage mismatch detected.")
  }

  # Ensure no technical identifiers changed.
  new_run_signature <- cra_new |>
    select(sample_id, experiment_id, run_id) |>
    arrange(sample_id)

  if (!identical(
    as.data.frame(old_run_signature),
    as.data.frame(new_run_signature)
  )) {
    stop("CRA technical Run/Experiment mapping changed during clinical repair.")
  }

  if (
    sum(!is.na(cra_new$outcome_28d)) != 131 ||
    sum(!is.na(cra_new$baseline_sepsis_status)) != 131
  ) {
    stop("CRA repaired patient-level clinical fields are not complete across 131 samples.")
  }

  if (
    !all(
      cra_new$phenotype ==
      cra_new$baseline_sepsis_status
    )
  ) {
    stop("CRA phenotype disagrees with repaired baseline sepsis status.")
  }

  ck("CRA CLINICAL REPAIR COMPLETE")

  # ==========================================================
  # 7. Recombine master
  # ==========================================================

  # Character-align all common columns before bind, then restore
  # numeric/logical fields below. This avoids all-NA type clashes.
  all_cols <- union(names(other), names(cra_new))

  add_missing <- function(df, cols) {
    miss <- setdiff(cols,names(df))
    for (nm in miss) df[[nm]] <- NA
    df[,cols,drop=FALSE]
  }

  other2 <- add_missing(other, all_cols)
  cra2 <- add_missing(cra_new, all_cols)

  # For bind safety, cast only columns with mixed types to character.
  common_numeric <- c(
    "time_day","time_order","age","gestational_age","birthweight",
    "baseline_sofa","baseline_apache_ii","baseline_lactate",
    "n_distinct_timepoints","sequence_n_distinct_timepoints"
  )

  common_logical <- c(
    "sequence_analysis_include",
    "freeze_metadata_include",
    "primary_16s_sequence_include",
    "shotgun_external_include",
    "longitudinal_ge2","longitudinal_ge3",
    "sequence_longitudinal_ge2","sequence_longitudinal_ge3",
    "patient_sample_link_match"
  )

  for (nm in intersect(common_numeric,all_cols)) {
    other2[[nm]] <- clean_num(other2[[nm]])
    cra2[[nm]] <- clean_num(cra2[[nm]])
  }

  for (nm in setdiff(all_cols, c(common_numeric,common_logical))) {
    other2[[nm]] <- as.character(other2[[nm]])
    cra2[[nm]] <- as.character(cra2[[nm]])
  }

  # Explicit logical conversion helper
  as_lgl <- function(x) {
    if (is.logical(x)) return(x)
    z <- tolower(clean_chr(x))
    case_when(
      z %in% c("true","t","1","yes","y") ~ TRUE,
      z %in% c("false","f","0","no","n") ~ FALSE,
      TRUE ~ NA
    )
  }

  for (nm in intersect(common_logical,all_cols)) {
    other2[[nm]] <- as_lgl(other2[[nm]])
    cra2[[nm]] <- as_lgl(cra2[[nm]])
  }

  full2 <- bind_rows(other2,cra2)

  # Preserve deterministic project/sample order.
  full2 <- full2 |>
    arrange(
      project,
      patient_id,
      time_order,
      time_day,
      sample_id
    )

  # ==========================================================
  # 8. Re-run global/run guards
  # ==========================================================

  if (nrow(full2) != 2693 || n_distinct(full2$project) != 18) {
    stop("Post-repair global row/project count changed.")
  }

  if (sum(full2$primary_16s_sequence_include,na.rm=TRUE) != 1941) {
    stop("Post-repair primary 16S row count changed.")
  }

  dup_run <- full2 |>
    filter(!is.na(run_uid)) |>
    count(run_uid,name="n") |>
    filter(n>1)

  run_conflicts <- full2 |>
    filter(!is.na(run_uid)) |>
    group_by(run_uid) |>
    summarise(
      n_rows=n(),
      patients=n_distinct(patient_id,na.rm=TRUE),
      samples=n_distinct(sample_id,na.rm=TRUE),
      timepoints=n_distinct(time_key,na.rm=TRUE),
      .groups="drop"
    ) |>
    filter(
      patients>1 |
      samples>1 |
      timepoints>1
    )

  if (nrow(dup_run)>0 || nrow(run_conflicts)>0) {
    stop("Duplicate/conflicting Run appeared after CRA clinical repair.")
  }

  # ==========================================================
  # 9. Clinical repair QC
  # ==========================================================

  cra_patient_qc <- cra_new |>
    distinct(
      patient_id,
      baseline_sepsis_status,
      outcome_28d,
      infection_source_standard,
      sex,
      age,
      baseline_sofa,
      baseline_apache_ii,
      baseline_lactate,
      antibiotics
    )

  cra_qc <- tibble(
    metric=c(
      "samples",
      "patients",
      "mapped_unique_runs",
      "baseline_sepsis_patients",
      "baseline_septic_shock_patients",
      "28d_survived_patients",
      "28d_dead_patients",
      "sample_status_nonmissing_rows",
      "infection_source_nonmissing_patients",
      "sex_nonmissing_patients",
      "age_nonmissing_patients",
      "SOFA_nonmissing_patients",
      "APACHEII_nonmissing_patients",
      "lactate_nonmissing_patients",
      "antibiotics_nonmissing_patients",
      "patient_sample_link_mismatches"
    ),
    value=c(
      nrow(cra_new),
      n_distinct(cra_new$patient_id),
      n_distinct(cra_new$run_id),
      sum(cra_patient_qc$baseline_sepsis_status=="sepsis",na.rm=TRUE),
      sum(cra_patient_qc$baseline_sepsis_status=="septic_shock",na.rm=TRUE),
      sum(cra_patient_qc$outcome_28d=="survived",na.rm=TRUE),
      sum(cra_patient_qc$outcome_28d=="dead",na.rm=TRUE),
      sum(!is.na(cra_new$sample_sepsis_status)),
      sum(!is.na(cra_patient_qc$infection_source_standard)),
      sum(!is.na(cra_patient_qc$sex)),
      sum(!is.na(cra_patient_qc$age)),
      sum(!is.na(cra_patient_qc$baseline_sofa)),
      sum(!is.na(cra_patient_qc$baseline_apache_ii)),
      sum(!is.na(cra_patient_qc$baseline_lactate)),
      sum(!is.na(cra_patient_qc$antibiotics)),
      sum(!cra_new$patient_sample_link_match,na.rm=TRUE)
    )
  )

  write_excel_csv(
    cra_qc,
    file.path(OUT,"CRA002354_FINAL_CLINICAL_REPAIR_QC.csv"),
    na=""
  )

  write_excel_csv(
    cra_patient_qc,
    file.path(OUT,"CRA002354_FINAL_PATIENT_LEVEL_CLINICAL.csv"),
    na=""
  )

  write_excel_csv(
    cra_new |>
      select(
        project,patient_id,sample_id,experiment_id,run_id,
        time_raw,time_day,time_order,
        baseline_sepsis_status,
        sample_sepsis_status,
        outcome_28d,
        infection_source_standard,
        baseline_sofa,
        baseline_apache_ii,
        baseline_lactate,
        antibiotics,
        clinical_mapping_status
      ),
    file.path(OUT,"CRA002354_FINAL_SAMPLE_LEVEL_CLINICAL_MAP.csv"),
    na=""
  )

  # ==========================================================
  # 10. Regenerate summaries/manifests
  # ==========================================================

  patient_summary <- full2 |>
    filter(!is.na(patient_uid)) |>
    group_by(project,patient_uid,patient_id,analysis_module) |>
    summarise(
      n_rows=n(),
      n_samples=n_distinct(sample_uid,na.rm=TRUE),
      n_runs=n_distinct(run_uid,na.rm=TRUE),
      n_timepoints=n_distinct(time_key[!is.na(time_key)]),
      n_sequence_timepoints=n_distinct(
        time_key[
          primary_16s_sequence_include &
          !is.na(time_key)
        ]
      ),
      metadata_longitudinal_ge2=n_timepoints>=2,
      metadata_longitudinal_ge3=n_timepoints>=3,
      sequence_longitudinal_ge2=n_sequence_timepoints>=2,
      sequence_longitudinal_ge3=n_sequence_timepoints>=3,
      .groups="drop"
    )

  cohort_summary <- full2 |>
    group_by(project,analysis_module,data_modality) |>
    summarise(
      rows=n(),
      patients=n_distinct(patient_uid,na.rm=TRUE),
      samples=n_distinct(sample_uid,na.rm=TRUE),
      mapped_runs=n_distinct(run_uid,na.rm=TRUE),
      missing_run_rows=sum(is.na(run_id)),
      metadata_patients_ge2=n_distinct(
        patient_uid[longitudinal_ge2],
        na.rm=TRUE
      ),
      metadata_patients_ge3=n_distinct(
        patient_uid[longitudinal_ge3],
        na.rm=TRUE
      ),
      sequence_patients_ge2=n_distinct(
        patient_uid[sequence_longitudinal_ge2],
        na.rm=TRUE
      ),
      sequence_patients_ge3=n_distinct(
        patient_uid[sequence_longitudinal_ge3],
        na.rm=TRUE
      ),
      primary_16s_rows=sum(primary_16s_sequence_include,na.rm=TRUE),
      shotgun_rows=sum(shotgun_external_include,na.rm=TRUE),
      outcome_nonmissing=sum(!is.na(outcome)),
      infection_group_nonmissing=sum(!is.na(infection_group)),
      .groups="drop"
    )

  module_summary <- full2 |>
    group_by(analysis_module) |>
    summarise(
      projects=n_distinct(project),
      patients=n_distinct(patient_uid,na.rm=TRUE),
      samples=n_distinct(sample_uid,na.rm=TRUE),
      mapped_runs=n_distinct(run_uid,na.rm=TRUE),
      metadata_patients_ge2=n_distinct(
        patient_uid[longitudinal_ge2],
        na.rm=TRUE
      ),
      sequence_patients_ge2=n_distinct(
        patient_uid[sequence_longitudinal_ge2],
        na.rm=TRUE
      ),
      primary_16s_rows=sum(primary_16s_sequence_include,na.rm=TRUE),
      shotgun_rows=sum(shotgun_external_include,na.rm=TRUE),
      .groups="drop"
    )

  registry <- cohort_summary |>
    mutate(
      freeze_decision=case_when(
        project=="PRJNA1125274" ~ "METADATA_ONLY_DEFERRED",
        data_modality=="SHOTGUN" ~ "SHOTGUN_EXTERNAL_SEPARATE_LAYER",
        analysis_module=="STATIC_SUPPORT" ~ "STATIC_SUPPORT",
        TRUE ~ "INCLUDE"
      ),
      next_sequence_action=case_when(
        project=="PRJNA1125274" ~
          "Finish correct 289-Run raw download/mapping before DADA2",
        project=="PRJNA884103" ~
          "Shotgun external layer; do not enter DADA2/16S workflow",
        project=="CRA002354" ~
          "Clinical metadata repaired; process cohort-specific sequence layer",
        project=="PRJNA851469" ~
          "Use only 119 strict exact-mapped samples in 16S sequence layer",
        TRUE ~
          "Audit existing ASV and raw sequence readiness"
      )
    )

  manifest16s <- full2 |>
    filter(primary_16s_sequence_include) |>
    select(
      project,patient_id,patient_uid,sample_id,sample_uid,
      biosample,experiment_id,run_id,run_uid,
      time_raw,time_class,time_day,time_order,
      analysis_module,cohort_role,sequence_status,freeze_status
    ) |>
    arrange(project,patient_id,time_order,time_day,sample_id)

  manifest16s_long <- full2 |>
    filter(
      primary_16s_sequence_include,
      sequence_longitudinal_ge2
    ) |>
    select(
      project,patient_id,patient_uid,sample_id,sample_uid,
      run_id,run_uid,time_raw,time_class,time_day,time_order,
      analysis_module,cohort_role,sequence_n_distinct_timepoints
    ) |>
    arrange(project,patient_id,time_order,time_day,sample_id)

  shotgun_manifest <- full2 |>
    filter(shotgun_external_include) |>
    select(
      project,patient_id,patient_uid,sample_id,sample_uid,
      biosample,run_id,run_uid,time_raw,time_class,time_day,time_order,
      analysis_module,cohort_role,sequence_status,freeze_status
    ) |>
    arrange(project,patient_id,time_order,time_day,sample_id)

  deferred <- full2 |>
    filter(
      project=="PRJNA1125274" |
      (
        data_modality=="16S" &
        !primary_16s_sequence_include
      ) |
      data_modality=="SHOTGUN"
    ) |>
    select(
      project,patient_id,sample_id,run_id,time_raw,
      analysis_module,data_modality,sequence_status,
      freeze_status,notes
    )

  global_summary <- tibble(
    metric=c(
      "Rows","Projects","Patients","Samples","Mapped_Runs",
      "Metadata_Patients_GE2","Metadata_Patients_GE3",
      "Primary_16S_Sequence_Rows",
      "Primary_16S_Patients_GE2","Primary_16S_Patients_GE3",
      "Shotgun_Metadata_Rows",
      "Deferred_PRJNA1125274_Rows",
      "Strict_Unmapped_16S_Rows",
      "Missing_Run_Rows",
      "Duplicate_Run_UID",
      "Run_Conflict_Rows",
      "CRA002354_Outcome_Nonmissing_Rows"
    ),
    value=c(
      nrow(full2),
      n_distinct(full2$project),
      n_distinct(full2$patient_uid,na.rm=TRUE),
      n_distinct(full2$sample_uid,na.rm=TRUE),
      n_distinct(full2$run_uid,na.rm=TRUE),
      n_distinct(full2$patient_uid[full2$longitudinal_ge2],na.rm=TRUE),
      n_distinct(full2$patient_uid[full2$longitudinal_ge3],na.rm=TRUE),
      sum(full2$primary_16s_sequence_include,na.rm=TRUE),
      n_distinct(full2$patient_uid[full2$sequence_longitudinal_ge2],na.rm=TRUE),
      n_distinct(full2$patient_uid[full2$sequence_longitudinal_ge3],na.rm=TRUE),
      sum(full2$shotgun_external_include,na.rm=TRUE),
      sum(full2$project=="PRJNA1125274"),
      sum(
        full2$data_modality=="16S" &
        full2$project!="PRJNA1125274" &
        !full2$primary_16s_sequence_include
      ),
      sum(is.na(full2$run_id)),
      nrow(dup_run),
      nrow(run_conflicts),
      sum(
        full2$project=="CRA002354" &
        !is.na(full2$outcome_28d)
      )
    )
  )

  # ==========================================================
  # 11. Final hard guards
  # ==========================================================

  expected <- c(
    Rows=2693,
    Projects=18,
    Primary_16S_Sequence_Rows=1941,
    Primary_16S_Patients_GE2=548,
    Primary_16S_Patients_GE3=210,
    Shotgun_Metadata_Rows=462,
    Deferred_PRJNA1125274_Rows=286,
    Strict_Unmapped_16S_Rows=4,
    Duplicate_Run_UID=0,
    Run_Conflict_Rows=0,
    CRA002354_Outcome_Nonmissing_Rows=131
  )

  for (nm in names(expected)) {
    observed <- global_summary$value[
      global_summary$metric==nm
    ]

    if (
      length(observed)!=1 ||
      observed != expected[[nm]]
    ) {
      stop(
        paste0(
          "Step81D final guard failed for ",
          nm,
          ": observed=",
          paste(observed,collapse=","),
          " expected=",
          expected[[nm]]
        )
      )
    }
  }

  ck("FINAL GLOBAL AND CRA CLINICAL GUARDS PASSED")

  # ==========================================================
  # 12. Write final freeze
  # ==========================================================

  stamp <- format(Sys.time(),"%Y%m%d_%H%M%S")
  version <- paste0("STEP81D_CLINICAL_REPAIRED_",stamp)

  full2$freeze_version <- version

  master_out <- file.path(
    OUT,
    paste0(
      "V2_FINAL_FROZEN_master_metadata_CLINICAL_REPAIRED_",
      stamp,
      ".csv"
    )
  )

  write_excel_csv(full2,master_out,na="")

  write_excel_csv(
    patient_summary,
    file.path(
      OUT,
      paste0(
        "V2_FINAL_FROZEN_patient_summary_",
        stamp,
        ".csv"
      )
    ),
    na=""
  )

  write_excel_csv(
    cohort_summary,
    file.path(
      OUT,
      paste0(
        "V2_FINAL_FROZEN_cohort_summary_",
        stamp,
        ".csv"
      )
    ),
    na=""
  )

  write_excel_csv(
    global_summary,
    file.path(
      OUT,
      paste0(
        "V2_FINAL_FROZEN_global_summary_",
        stamp,
        ".csv"
      )
    ),
    na=""
  )

  write_excel_csv(
    module_summary,
    file.path(
      OUT,
      "V2_FINAL_FROZEN_analysis_module_summary.csv"
    ),
    na=""
  )

  write_excel_csv(
    registry,
    file.path(
      OUT,
      "V2_FINAL_FROZEN_project_registry.csv"
    ),
    na=""
  )

  write_excel_csv(
    manifest16s,
    file.path(
      OUT,
      "V2_FINAL_FROZEN_16S_sequence_manifest.csv"
    ),
    na=""
  )

  write_excel_csv(
    manifest16s_long,
    file.path(
      OUT,
      "V2_FINAL_FROZEN_16S_longitudinal_manifest.csv"
    ),
    na=""
  )

  write_excel_csv(
    shotgun_manifest,
    file.path(
      OUT,
      "V2_FINAL_FROZEN_SHOTGUN_external_manifest.csv"
    ),
    na=""
  )

  write_excel_csv(
    deferred,
    file.path(
      OUT,
      "V2_FINAL_FROZEN_deferred_registry.csv"
    ),
    na=""
  )

  write_excel_csv(
    dup_run,
    file.path(
      OUT,
      "V2_FINAL_FROZEN_duplicate_run_QC.csv"
    ),
    na=""
  )

  write_excel_csv(
    run_conflicts,
    file.path(
      OUT,
      "V2_FINAL_FROZEN_run_conflict_QC.csv"
    ),
    na=""
  )

  readme <- c(
    "SEPSIS V2 - STEP81D FINAL FREEZE WITH CRA002354 CLINICAL REPAIR",
    paste0("Freeze version: ",version),
    paste0("Created: ",Sys.time()),
    "",
    "WHY STEP81D WAS REQUIRED",
    "Step81C correctly fixed the 16S/SHOTGUN modality and sequence-layer inclusion.",
    "During final audit, CRA002354 clinical outcome/severity fields in the Step81C master were found blank.",
    "Step81D re-read the validated original CRA002354 Supplementary Tables S1/S2 and repaired those fields by exact Patient ID / Sample ID joins.",
    "",
    "CRA002354 VERIFIED PATIENT-LEVEL COUNTS",
    "- Patients: 64",
    "- Baseline sepsis: 46",
    "- Baseline septic shock: 18",
    "- 28-day survived: 44",
    "- 28-day dead: 20",
    "- Technical Run map remains 131/131 exact and unchanged.",
    "",
    "GLOBAL FREEZE COUNTS UNCHANGED",
    "- Rows: 2693",
    "- Projects: 18",
    "- Primary 16S sequence rows: 1941",
    "- Primary 16S patients >=2 timepoints: 548",
    "- Primary 16S patients >=3 timepoints: 210",
    "- Shotgun metadata rows: 462",
    "- PRJNA1125274 deferred rows: 286",
    "- Strict unmapped eligible 16S rows: 4",
    "- Duplicate Run UID: 0",
    "- Run conflicts: 0",
    "",
    "UNIQUE DOWNSTREAM METADATA ENTRY POINT",
    master_out,
    "",
    "NEXT STEP",
    "Step82 may now audit ASV/raw-sequence readiness by cohort.",
    "Use only the Step81D master/manifests for downstream work."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP81D_FINAL_FREEZE.txt"
    ),
    useBytes=TRUE
  )

  ck("STEP81D COMPLETE")

  cat("\n============================================================\n")
  cat("SEPSIS V2 - STEP81D FINAL FREEZE COMPLETE\n")
  cat("============================================================\n\n")
  cat("CRA002354 clinical QC:\n")
  print(cra_qc,n=Inf,width=Inf)
  cat("\nGlobal summary:\n")
  print(global_summary,n=Inf,width=Inf)
  cat("\nOutput:\n",OUT,"\n",sep="")
}

tryCatch(
  main(),
  error=function(e) {

    msg <- c(
      paste0("STEP81D FATAL ERROR: ",Sys.time()),
      paste0("Message: ",conditionMessage(e)),
      paste0(
        "Call: ",
        paste(
          deparse(conditionCall(e)),
          collapse=" "
        )
      )
    )

    writeLines(msg,ERR,useBytes=TRUE)
    ck("STEP81D FAILED")
    message(paste(msg,collapse="\n"))

    quit(
      save="no",
      status=1,
      runLast=FALSE
    )
  }
)
