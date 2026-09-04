# ============================================================
# Sepsis V2 - Step 80F
# Exact CRA002354 two-hop technical mapping
# + strict PRJNA851469 exact mapping freeze
#
# OFFLINE ONLY
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c("readr","readxl","dplyr","stringr","tibble")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) install.packages(missing, repos="https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"
IN80C <- file.path(ROOT, "results", "V2_20C_run_mapping_and_clinical_fix")
OUT <- file.path(ROOT, "results", "V2_20F_exact_final_mapping")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

LOG <- file.path(OUT, "_STEP80F_runtime_checkpoints.txt")
ERR <- file.path(OUT, "_STEP80F_FATAL_ERROR.txt")
if (file.exists(ERR)) unlink(ERR)

ck <- function(x) {
  cat(paste0(x, ": ", Sys.time(), "\n"), file=LOG, append=TRUE)
}

clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x=="" | tolower(x) %in% c("na","nan","n/a","null","none")] <- NA_character_
  x
}

norm_id <- function(x) {
  x <- toupper(clean_chr(x))
  str_replace_all(x, "[^A-Z0-9]", "")
}

safe_csv <- function(p) {
  tryCatch(
    suppressMessages(read_csv(p, show_col_types=FALSE, progress=FALSE)),
    error=function(e) NULL
  )
}

find_exact_col <- function(nms, target) {
  z <- which(toupper(gsub("[^A-Za-z0-9]","",nms)) ==
               toupper(gsub("[^A-Za-z0-9]","",target)))
  if (!length(z)) NA_character_ else nms[z[1]]
}

main <- function() {

  ck("STEP80F STARTED")

  # ==========================================================
  # Inputs from Step80C
  # ==========================================================
  cra_path <- file.path(IN80C, "CRA002354_MASTER_READY_CORRECTED.csv")
  p851_path <- file.path(IN80C, "PRJNA851469_MASTER_READY_WITH_RUN_SEARCH.csv")

  if (!file.exists(cra_path) || !file.exists(p851_path)) {
    stop("Required Step80C files are missing.")
  }

  cra <- safe_csv(cra_path)
  p851 <- safe_csv(p851_path)

  if (is.null(cra) || is.null(p851)) stop("Failed to read Step80C metadata.")

  # ==========================================================
  # A. Locate official CRA002354.xlsx
  # ==========================================================
  candidates <- c(
    file.path(Sys.getenv("USERPROFILE"), "Downloads", "CRA002354.xlsx"),
    file.path(ROOT, "data", "CRA002354", "CRA002354.xlsx")
  )

  candidates <- c(
    candidates,
    list.files(
      ROOT,
      recursive=TRUE,
      full.names=TRUE,
      pattern="^CRA002354\\.xlsx$",
      ignore.case=TRUE
    )
  )

  candidates <- unique(candidates[file.exists(candidates)])

  if (!length(candidates)) {
    stop("Official CRA002354.xlsx was not found.")
  }

  # Prefer the workbook that contains Experiment + Run sheets.
  wb <- NA_character_

  for (p in candidates) {
    sh <- tryCatch(excel_sheets(p), error=function(e) character())
    if (all(c("Experiment","Run") %in% sh)) {
      wb <- p
      break
    }
  }

  if (is.na(wb)) {
    stop("A CRA002354.xlsx was found, but it lacks Experiment and Run sheets.")
  }

  ck(paste0("CRA WORKBOOK FOUND = ", wb))

  # ==========================================================
  # B. Exact two-hop CRA mapping
  #
  # Supplement sample_id
  #   -> Experiment sheet: BioSample name
  #   -> Experiment accession
  #   -> Run sheet: Experiment accession
  #   -> CRR Run accession
  # ==========================================================
  exp_tbl <- read_excel(wb, sheet="Experiment", .name_repair="unique") |>
    mutate(across(everything(), clean_chr))

  run_tbl <- read_excel(wb, sheet="Run", .name_repair="unique") |>
    mutate(across(everything(), clean_chr))

  exp_sample_col <- find_exact_col(names(exp_tbl), "BioSample name")
  exp_acc_col <- find_exact_col(names(exp_tbl), "Accession")
  exp_biosample_acc_col <- find_exact_col(names(exp_tbl), "BioSample accession")

  run_exp_col <- find_exact_col(names(run_tbl), "Experiment accession")
  run_acc_col <- find_exact_col(names(run_tbl), "Accession")

  needed <- c(exp_sample_col, exp_acc_col, run_exp_col, run_acc_col)

  if (any(is.na(needed))) {
    stop(
      paste0(
        "Required CRA technical columns missing. Experiment cols: ",
        paste(names(exp_tbl), collapse=" | "),
        " ; Run cols: ",
        paste(names(run_tbl), collapse=" | ")
      )
    )
  }

  exp_map <- tibble(
    sample_id_technical=clean_chr(exp_tbl[[exp_sample_col]]),
    sample_norm=norm_id(exp_tbl[[exp_sample_col]]),
    experiment_id=clean_chr(exp_tbl[[exp_acc_col]]),
    archive_sample_accession=
      if (!is.na(exp_biosample_acc_col)) {
        clean_chr(exp_tbl[[exp_biosample_acc_col]])
      } else {
        NA_character_
      }
  ) |>
    filter(
      !is.na(sample_norm),
      !is.na(experiment_id)
    ) |>
    distinct()

  run_map <- tibble(
    experiment_id=clean_chr(run_tbl[[run_exp_col]]),
    run_id=clean_chr(run_tbl[[run_acc_col]])
  ) |>
    filter(
      !is.na(experiment_id),
      !is.na(run_id),
      str_detect(run_id, regex("^CRR[0-9]+$", ignore_case=TRUE))
    ) |>
    distinct()

  technical_map <- exp_map |>
    left_join(run_map, by="experiment_id")

  # Strict QC: 131 unique samples -> 131 unique experiments -> 131 unique runs
  tech_qc <- tibble(
    metric=c(
      "technical_rows",
      "unique_sample_ids",
      "unique_experiments",
      "run_nonmissing",
      "unique_runs",
      "duplicate_sample_ids",
      "duplicate_experiment_ids",
      "duplicate_run_ids"
    ),
    value=c(
      nrow(technical_map),
      n_distinct(technical_map$sample_norm),
      n_distinct(technical_map$experiment_id),
      sum(!is.na(technical_map$run_id)),
      n_distinct(technical_map$run_id, na.rm=TRUE),
      sum(duplicated(technical_map$sample_norm)),
      sum(duplicated(technical_map$experiment_id)),
      sum(duplicated(na.omit(technical_map$run_id)))
    )
  )

  write_excel_csv(
    technical_map,
    file.path(OUT, "CRA002354_EXACT_technical_map.csv"),
    na=""
  )

  write_excel_csv(
    tech_qc,
    file.path(OUT, "CRA002354_technical_map_QC.csv"),
    na=""
  )

  cra_final <- cra |>
    mutate(sample_norm=norm_id(sample_id)) |>
    select(
      -any_of(c(
        "run_id",
        "run_mapping_source_file",
        "run_mapping_source_sheet"
      ))
    ) |>
    left_join(
      technical_map |>
        select(
          sample_norm,
          experiment_id,
          archive_sample_accession,
          run_id
        ),
      by="sample_norm"
    ) |>
    mutate(
      run_mapping_status=ifelse(
        !is.na(run_id),
        "EXACT_TWO_HOP_GSA_WORKBOOK",
        "UNMAPPED"
      ),
      run_mapping_source=wb
    )

  cra_qc <- tibble(
    metric=c(
      "supplement_samples",
      "mapped_runs",
      "missing_runs",
      "unique_runs",
      "duplicate_run_assignments"
    ),
    value=c(
      nrow(cra_final),
      sum(!is.na(cra_final$run_id)),
      sum(is.na(cra_final$run_id)),
      n_distinct(cra_final$run_id, na.rm=TRUE),
      sum(duplicated(na.omit(cra_final$run_id)))
    )
  )

  write_excel_csv(
    cra_final,
    file.path(OUT, "CRA002354_MASTER_READY_FINAL.csv"),
    na=""
  )

  write_excel_csv(
    cra_qc,
    file.path(OUT, "CRA002354_FINAL_QC.csv"),
    na=""
  )

  ck("CRA002354 EXACT TWO-HOP MAPPING COMPLETE")

  # ==========================================================
  # C. PRJNA851469 strict exact freeze
  #
  # Step80C already contains the safe exact mappings.
  # Do not use Step80D substring containment.
  # ==========================================================
  p851_final <- p851 |>
    mutate(
      run_mapping_status=ifelse(
        !is.na(run_id),
        "EXACT_CONFIRMED",
        "UNMAPPED"
      ),
      sequence_available=!is.na(run_id),
      sequence_analysis_include=!is.na(run_id)
    )

  p851_unmatched <- p851_final |>
    filter(is.na(run_id)) |>
    select(
      sample_id,
      sample_class,
      patient_id,
      time_raw,
      time_day,
      run_mapping_status
    )

  icu_seq <- p851_final |>
    filter(
      sample_class=="ICU_PATIENT",
      sequence_analysis_include
    )

  pt851 <- icu_seq |>
    distinct(patient_id,time_day) |>
    count(patient_id,name="n_timepoints")

  p851_qc <- tibble(
    metric=c(
      "supplement_samples",
      "exact_mapped_runs",
      "unmapped_samples",
      "unique_runs",
      "ICU_sequence_samples",
      "ICU_sequence_patients",
      "ICU_sequence_patients_GE2",
      "ICU_sequence_patients_GE3"
    ),
    value=c(
      nrow(p851_final),
      sum(!is.na(p851_final$run_id)),
      sum(is.na(p851_final$run_id)),
      n_distinct(p851_final$run_id,na.rm=TRUE),
      nrow(icu_seq),
      n_distinct(icu_seq$patient_id),
      sum(pt851$n_timepoints>=2),
      sum(pt851$n_timepoints>=3)
    )
  )

  write_excel_csv(
    p851_final,
    file.path(OUT, "PRJNA851469_MASTER_READY_FINAL.csv"),
    na=""
  )

  write_excel_csv(
    p851_unmatched,
    file.path(OUT, "PRJNA851469_STRICT_UNMATCHED_SAMPLES.csv"),
    na=""
  )

  write_excel_csv(
    p851_qc,
    file.path(OUT, "PRJNA851469_FINAL_QC.csv"),
    na=""
  )

  ck("PRJNA851469 STRICT EXACT FREEZE COMPLETE")

  # ==========================================================
  # Final decision
  # ==========================================================
  decision <- tibble(
    project=c("CRA002354","PRJNA851469"),
    patient_time_metadata_ready=c(TRUE,TRUE),
    mapped_runs=c(
      sum(!is.na(cra_final$run_id)),
      sum(!is.na(p851_final$run_id))
    ),
    total_samples=c(
      nrow(cra_final),
      nrow(p851_final)
    ),
    run_mapping_status=c(
      ifelse(
        all(!is.na(cra_final$run_id)),
        "COMPLETE_EXACT",
        "PARTIAL"
      ),
      "STRICT_EXACT_ONLY"
    ),
    recommended_role=c(
      "CORE_SEPSIS_LONGITUDINAL_WITH_CLINICAL_OUTCOMES",
      "ICU_BACKGROUND_LONGITUDINAL"
    ),
    ready_for_step81=c(
      all(!is.na(cra_final$run_id)),
      TRUE
    )
  )

  write_excel_csv(
    decision,
    file.path(OUT, "V2_step80F_decision.csv"),
    na=""
  )

  ck("STEP80F COMPLETE")

  cat("\n============================================\n")
  cat("STEP80F COMPLETE\n")
  cat("============================================\n")
  cat("\nCRA002354:\n")
  print(cra_qc,n=Inf)
  cat("\nPRJNA851469:\n")
  print(p851_qc,n=Inf)
  cat("\nDecision:\n")
  print(decision,n=Inf)
  cat("\nOutput: ",OUT,"\n",sep="")
}

tryCatch(
  main(),
  error=function(e) {
    msg <- c(
      paste0("STEP80F FATAL ERROR: ",Sys.time()),
      paste0("Message: ",conditionMessage(e)),
      paste0("Call: ",paste(deparse(conditionCall(e)),collapse=" "))
    )
    writeLines(msg,ERR,useBytes=TRUE)
    ck("STEP80F FAILED")
    message(paste(msg,collapse="\n"))
    quit(save="no",status=1,runLast=FALSE)
  }
)
