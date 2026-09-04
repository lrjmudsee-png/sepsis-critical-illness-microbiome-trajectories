# ============================================================
# Sepsis V2 - Step 80B
# Robust OFFLINE core metadata finalization
# CRA002354 + PRJNA851469
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c("readr","readxl","dplyr","stringr","tibble")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) install.packages(missing, repos="https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(readr); library(readxl); library(dplyr)
  library(stringr); library(tibble)
})

ROOT <- "E:/sepsis_project"
OUT <- file.path(ROOT, "results", "V2_20B_core_metadata_finalize")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

LOG <- file.path(OUT, "_STEP80B_runtime_checkpoints.txt")
ERR <- file.path(OUT, "_STEP80B_FATAL_ERROR.txt")
if (file.exists(ERR)) unlink(ERR)

ck <- function(x) cat(paste0(x, ": ", Sys.time(), "\n"), file=LOG, append=TRUE)

clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x=="" | tolower(x) %in% c("na","nan","n/a","null","none")] <- NA_character_
  x
}

norm_name <- function(x) toupper(str_replace_all(x, "[^A-Za-z0-9]", ""))

find_col <- function(nms, target) {
  z <- which(norm_name(nms) == norm_name(target))
  if (!length(z)) NA_character_ else nms[z[1]]
}

find_regex_col <- function(nms, pat) {
  z <- which(str_detect(nms, regex(pat, ignore_case=TRUE)))
  if (!length(z)) NA_character_ else nms[z[1]]
}

first_existing <- function(x) {
  z <- x[file.exists(x)]
  if (!length(z)) NA_character_ else z[1]
}

valid_xlsx <- function(path) {
  if (is.na(path) || !file.exists(path) || file.info(path)$size < 1000) return(FALSE)
  tryCatch(length(excel_sheets(path)) > 0, error=function(e) FALSE)
}

main <- function() {

  ck("STEP80B STARTED")

  cra_dir <- file.path(ROOT, "data", "_step78A3_US_fast_download", "CRA002354")
  cra_s1 <- first_existing(c(
    file.path(cra_dir, "CRA002354_Table_S1_64_patients.xlsx"),
    file.path(cra_dir, "mmc5.xlsx")
  ))
  cra_s2 <- first_existing(c(
    file.path(cra_dir, "CRA002354_Table_S2_131_samples.xlsx"),
    file.path(cra_dir, "mmc6.xlsx")
  ))
  p851 <- file.path(
    ROOT, "data", "_step78A3_US_fast_download", "PRJNA851469",
    "PRJNA851469_Supplementary_Tables_2_17.xlsx"
  )

  input_qc <- tibble(
    project=c("CRA002354","CRA002354","PRJNA851469"),
    source=c("S1_patient","S2_sample","Nature_Tables_2_17"),
    path=c(cra_s1,cra_s2,p851),
    valid_xlsx=c(valid_xlsx(cra_s1),valid_xlsx(cra_s2),valid_xlsx(p851))
  )

  write_excel_csv(input_qc, file.path(OUT,"STEP80B_input_validation.csv"), na="")
  if (!all(input_qc$valid_xlsx)) stop("One or more XLSX inputs are missing or invalid.")
  ck("INPUT VALIDATION COMPLETE")

  # ---------------- CRA002354 ----------------
  # True headers are on Excel row 3 -> skip first 2 rows.
  s1 <- read_excel(cra_s1, skip=2, .name_repair="unique") |>
    mutate(across(everything(), clean_chr))
  s2 <- read_excel(cra_s2, skip=2, .name_repair="unique") |>
    mutate(across(everything(), clean_chr))

  write_excel_csv(s1, file.path(OUT,"CRA002354_S1_TRUE_HEADER.csv"), na="")
  write_excel_csv(s2, file.path(OUT,"CRA002354_S2_TRUE_HEADER.csv"), na="")

  s1_pid <- find_col(names(s1),"Patient ID")
  s1_day <- find_col(names(s1),"Collected day")
  s1_sid <- find_col(names(s1),"Sample ID")
  s2_pid <- find_col(names(s2),"Patient ID")
  s2_day <- find_col(names(s2),"Collected day")
  s2_sid <- find_col(names(s2),"Sample ID")

  if (any(is.na(c(s1_pid,s1_day,s1_sid,s2_pid,s2_day,s2_sid)))) {
    stop(paste0(
      "CRA required columns missing. S1=",paste(names(s1),collapse=" | "),
      " ; S2=",paste(names(s2),collapse=" | ")
    ))
  }

  if (nrow(s1) != 64) stop(paste0("CRA S1 expected 64 rows; got ",nrow(s1)))
  if (nrow(s2) != 131) stop(paste0("CRA S2 expected 131 rows; got ",nrow(s2)))

  patients <- s1 |>
    rename(
      patient_id=all_of(s1_pid),
      first_time_raw=all_of(s1_day),
      first_sample_id=all_of(s1_sid)
    ) |>
    mutate(project="CRA002354", .before=1)

  samples <- s2 |>
    rename(
      patient_id=all_of(s2_pid),
      time_raw=all_of(s2_day),
      sample_id=all_of(s2_sid)
    ) |>
    mutate(
      project="CRA002354",
      time_day=suppressWarnings(as.integer(str_extract(time_raw,"[0-9]+"))),
      time_order=time_day,
      .before=1
    )

  cra <- samples |>
    left_join(patients, by=c("project","patient_id"), suffix=c("_sample","_patient"))

  # Preserve clinical columns and add standardized fields where public codes exist.
  sepsis_col <- find_regex_col(names(cra),"Sepsis.*Septic shock")
  survival_col <- find_regex_col(names(cra),"28.*day.*survival")
  infection_col <- find_regex_col(names(cra),"^Site of infection")

  infection_codes <- c(
    "0"="other","1"="lung","2"="intestinal","3"="abdominal",
    "4"="blood","5"="urinary","6"="brain","7"="surgical_site"
  )

  cra <- cra |>
    mutate(
      sepsis_status = if (!is.na(sepsis_col)) {
        case_when(
          .data[[sepsis_col]]=="1" ~ "sepsis",
          .data[[sepsis_col]]=="2" ~ "septic_shock",
          TRUE ~ NA_character_
        )
      } else NA_character_,
      outcome_28d = if (!is.na(survival_col)) {
        case_when(
          .data[[survival_col]]=="1" ~ "survived",
          .data[[survival_col]]=="2" ~ "dead",
          TRUE ~ NA_character_
        )
      } else NA_character_,
      infection_source = if (!is.na(infection_col)) {
        unname(infection_codes[clean_chr(.data[[infection_col]])])
      } else NA_character_,
      cohort_role="CORE_SEPSIS_LONGITUDINAL",
      patient_mapping_status="EXPLICIT_PUBLIC_PATIENT_ID",
      time_mapping_status="EXPLICIT_PUBLIC_COLLECTED_DAY",
      metadata_source="CRA002354 Supplementary Tables S1/S2",
      run_id=NA_character_
    )

  pt_cra <- cra |> distinct(patient_id,sample_id) |> count(patient_id,name="n")
  cra_qc <- tibble(
    metric=c("samples","patients","patients_GE2","patients_GE3","max_timepoints",
             "sepsis_patients","septic_shock_patients",
             "28d_survived_patients","28d_dead_patients"),
    value=c(
      nrow(cra), n_distinct(cra$patient_id),
      sum(pt_cra$n>=2), sum(pt_cra$n>=3), max(pt_cra$n),
      n_distinct(cra$patient_id[cra$sepsis_status=="sepsis"]),
      n_distinct(cra$patient_id[cra$sepsis_status=="septic_shock"]),
      n_distinct(cra$patient_id[cra$outcome_28d=="survived"]),
      n_distinct(cra$patient_id[cra$outcome_28d=="dead"])
    )
  )

  write_excel_csv(cra, file.path(OUT,"CRA002354_MASTER_READY_core_metadata.csv"), na="")
  write_excel_csv(patients, file.path(OUT,"CRA002354_patient_level_clinical.csv"), na="")
  write_excel_csv(cra_qc, file.path(OUT,"CRA002354_FINAL_QC.csv"), na="")
  ck("CRA002354 COMPLETE")

  # ---------------- PRJNA851469 ----------------
  sh <- excel_sheets(p851)
  hit <- sh[str_detect(sh, regex("^\\s*Supplementary\\s+table\\s+2\\s*$",ignore_case=TRUE))]
  sheet2 <- if (length(hit)) hit[1] else sh[2]

  t2 <- read_excel(p851, sheet=sheet2, .name_repair="unique") |>
    mutate(across(everything(), clean_chr))

  sid_col <- find_col(names(t2),"Sample ID")
  if (is.na(sid_col)) stop(paste0("PRJNA851469 Sample ID missing. Columns=",paste(names(t2),collapse=" | ")))

  pmap <- t2 |>
    distinct(sample_id=.data[[sid_col]]) |>
    mutate(
      sample_class=case_when(
        str_detect(sample_id,regex("^HealthyVolunteer-",ignore_case=TRUE)) ~ "HEALTHY_CONTROL",
        str_detect(sample_id,regex("^Patient[0-9]+-Day-[0-9]+$",ignore_case=TRUE)) ~ "ICU_PATIENT",
        TRUE ~ "OTHER"
      ),
      patient_id=case_when(
        sample_class=="ICU_PATIENT" ~ str_extract(sample_id,regex("^Patient[0-9]+",ignore_case=TRUE)),
        sample_class=="HEALTHY_CONTROL" ~ sample_id,
        TRUE ~ NA_character_
      ),
      time_raw=case_when(
        sample_class=="ICU_PATIENT" ~ str_extract(sample_id,regex("Day-[0-9]+",ignore_case=TRUE)),
        TRUE ~ NA_character_
      ),
      time_day=suppressWarnings(as.integer(str_extract(time_raw,"[0-9]+"))),
      time_order=time_day,
      cohort_role=case_when(
        sample_class=="ICU_PATIENT" ~ "ICU_BACKGROUND_LONGITUDINAL",
        sample_class=="HEALTHY_CONTROL" ~ "HEALTHY_STATIC_CONTROL",
        TRUE ~ "UNCLASSIFIED"
      ),
      patient_mapping_status="EXPLICIT_PUBLIC_SAMPLE_ID",
      time_mapping_status=ifelse(
        sample_class=="ICU_PATIENT",
        "EXPLICIT_PUBLIC_SAMPLE_ID_DAY",
        "NOT_APPLICABLE"
      ),
      metadata_source="Nature Medicine Supplementary Table 2",
      run_id=NA_character_
    )

  icu <- pmap |> filter(sample_class=="ICU_PATIENT")
  pt851 <- icu |> distinct(patient_id,time_day) |> count(patient_id,name="n")

  q851 <- tibble(
    metric=c("all_microbiome_samples","healthy_control_samples","ICU_microbiome_samples",
             "ICU_patients","Day1","Day3","Day7","patients_GE2","patients_GE3"),
    value=c(
      nrow(pmap),
      sum(pmap$sample_class=="HEALTHY_CONTROL"),
      nrow(icu),
      n_distinct(icu$patient_id),
      sum(icu$time_day==1,na.rm=TRUE),
      sum(icu$time_day==3,na.rm=TRUE),
      sum(icu$time_day==7,na.rm=TRUE),
      sum(pt851$n>=2),
      sum(pt851$n>=3)
    )
  )

  write_excel_csv(pmap, file.path(OUT,"PRJNA851469_MASTER_READY_core_metadata.csv"), na="")
  write_excel_csv(icu, file.path(OUT,"PRJNA851469_ICU_LONGITUDINAL_ONLY.csv"), na="")
  write_excel_csv(q851, file.path(OUT,"PRJNA851469_FINAL_QC.csv"), na="")
  ck("PRJNA851469 COMPLETE")

  decision <- tibble(
    project=c("CRA002354","PRJNA851469"),
    patient_time_mapping=c("EXPLICIT_CONFIRMED","EXPLICIT_CONFIRMED"),
    longitudinal_patients=c(sum(pt_cra$n>=2),sum(pt851$n>=2)),
    run_mapping_status=c("PENDING_STEP80C","PENDING_STEP80C"),
    recommended_role=c(
      "CORE_SEPSIS_LONGITUDINAL_WITH_CLINICAL_OUTCOMES",
      "ICU_BACKGROUND_LONGITUDINAL"
    ),
    metadata_decision=c("INCLUDE","INCLUDE")
  )

  write_excel_csv(decision,file.path(OUT,"V2_step80B_cohort_decision.csv"),na="")
  ck("STEP80B COMPLETE")

  cat("\nSTEP80B COMPLETE\n")
  print(cra_qc,n=Inf)
  print(q851,n=Inf)
  print(decision,n=Inf)
  cat("\nOutput: ",OUT,"\n",sep="")
}

tryCatch(
  main(),
  error=function(e) {
    msg <- c(
      paste0("STEP80B FATAL ERROR: ",Sys.time()),
      paste0("Message: ",conditionMessage(e)),
      paste0("Call: ",paste(deparse(conditionCall(e)),collapse=" "))
    )
    writeLines(msg,ERR,useBytes=TRUE)
    ck("STEP80B FAILED")
    message(paste(msg,collapse="\n"))
    quit(save="no",status=1,runLast=FALSE)
  }
)
