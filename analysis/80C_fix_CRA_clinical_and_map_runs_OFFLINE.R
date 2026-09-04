# ============================================================
# Sepsis V2 - Step 80C
# Correct CRA clinical coding + recover sample -> Run mapping
# for CRA002354 and PRJNA851469.
#
# OFFLINE ONLY
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c("readr","readxl","dplyr","stringr","purrr","tibble")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) install.packages(missing, repos="https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

ROOT <- "E:/sepsis_project"
IN <- file.path(ROOT, "results", "V2_20B_core_metadata_finalize")
OUT <- file.path(ROOT, "results", "V2_20C_run_mapping_and_clinical_fix")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

LOG <- file.path(OUT, "_STEP80C_runtime_checkpoints.txt")
ERR <- file.path(OUT, "_STEP80C_FATAL_ERROR.txt")
if (file.exists(ERR)) unlink(ERR)

ck <- function(x) cat(paste0(x, ": ", Sys.time(), "\n"), file=LOG, append=TRUE)

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

safe_tsv <- function(p) {
  tryCatch(
    suppressMessages(read_tsv(p, show_col_types=FALSE, progress=FALSE)),
    error=function(e) NULL
  )
}

run_count <- function(x) {
  x <- clean_chr(x)
  sum(str_detect(x, regex("^(SRR|ERR|DRR|CRR)[0-9]+$", ignore_case=TRUE)), na.rm=TRUE)
}

# Search local CSV/TSV/XLSX tables for a sample-ID column
# overlapping the target sample IDs AND a Run accession column.
search_mapping <- function(project, target_ids, roots) {

  target_norm <- unique(norm_id(target_ids))
  target_norm <- target_norm[!is.na(target_norm)]

  files <- unique(unlist(lapply(roots, function(r) {
    if (!dir.exists(r)) return(character())
    list.files(
      r,
      recursive=TRUE,
      full.names=TRUE,
      pattern="\\.(csv|tsv|xlsx)$",
      ignore.case=TRUE
    )
  })))

  files <- files[
    file.exists(files) &
    file.info(files)$size <= 50*1024^2
  ]

  # Prefer project-related files and exclude current finalization outputs.
  project_hits <- files[
    str_detect(files, regex(project, ignore_case=TRUE))
  ]
  if (length(project_hits)) files <- project_hits

  files <- files[
    !str_detect(
      files,
      regex(
        "V2_20B_core_metadata_finalize|V2_20C_run_mapping_and_clinical_fix|Supplementary_Tables_2_17|mmc[0-9]+\\.xlsx",
        ignore_case=TRUE
      )
    )
  ]

  scores <- list()
  maps <- list()

  inspect <- function(dat, file, sheet=NA_character_) {
    if (is.null(dat) || !nrow(dat) || !ncol(dat)) return(NULL)
    dat <- as_tibble(dat)

    overlap <- map_int(dat, function(v) {
      vv <- unique(norm_id(v))
      vv <- vv[!is.na(vv)]
      length(intersect(vv, target_norm))
    })

    runs <- map_int(dat, run_count)

    if (max(overlap, na.rm=TRUE) <= 0 || max(runs, na.rm=TRUE) <= 0) return(NULL)

    sample_col <- names(dat)[which.max(overlap)]
    run_col <- names(dat)[which.max(runs)]

    score <- tibble(
      project=project,
      source_file=file,
      source_sheet=sheet,
      rows=nrow(dat),
      cols=ncol(dat),
      sample_column=sample_col,
      exact_sample_overlap=max(overlap,na.rm=TRUE),
      run_column=run_col,
      run_accession_count=max(runs,na.rm=TRUE)
    )

    map <- tibble(
      sample_norm=norm_id(dat[[sample_col]]),
      run_id=clean_chr(dat[[run_col]]),
      run_mapping_source_file=file,
      run_mapping_source_sheet=sheet
    ) |>
      filter(
        !is.na(sample_norm),
        !is.na(run_id),
        str_detect(run_id, regex("^(SRR|ERR|DRR|CRR)[0-9]+$", ignore_case=TRUE))
      ) |>
      distinct()

    list(score=score,map=map)
  }

  for (p in files) {
    ext <- tolower(tools::file_ext(p))

    if (ext=="csv") {
      z <- inspect(safe_csv(p),p)
      if (!is.null(z)) {
        scores[[length(scores)+1]] <- z$score
        maps[[length(maps)+1]] <- z$map
      }
    }

    if (ext=="tsv") {
      z <- inspect(safe_tsv(p),p)
      if (!is.null(z)) {
        scores[[length(scores)+1]] <- z$score
        maps[[length(maps)+1]] <- z$map
      }
    }

    if (ext=="xlsx") {
      shs <- tryCatch(excel_sheets(p), error=function(e) character())
      for (sh in shs) {
        dat <- tryCatch(
          read_excel(p,sheet=sh,.name_repair="unique"),
          error=function(e) NULL
        )
        z <- inspect(dat,p,sh)
        if (!is.null(z)) {
          scores[[length(scores)+1]] <- z$score
          maps[[length(maps)+1]] <- z$map
        }
      }
    }
  }

  score_df <- if (!length(scores)) tibble() else bind_rows(scores) |>
    arrange(desc(exact_sample_overlap),desc(run_accession_count))

  best_map <- tibble(
    sample_norm=character(),
    run_id=character(),
    run_mapping_source_file=character(),
    run_mapping_source_sheet=character()
  )

  if (nrow(score_df)) {
    b <- score_df[1,]
    idx <- which(map_lgl(scores,function(s) {
      identical(s$source_file[1],b$source_file[1]) &&
      identical(s$source_sheet[1],b$source_sheet[1])
    }))
    if (length(idx)) best_map <- maps[[idx[1]]]
  }

  list(scores=score_df,map=best_map)
}

main <- function() {

  ck("STEP80C STARTED")

  cra_path <- file.path(IN,"CRA002354_MASTER_READY_core_metadata.csv")
  cra_pat_path <- file.path(IN,"CRA002354_patient_level_clinical.csv")
  p851_path <- file.path(IN,"PRJNA851469_MASTER_READY_core_metadata.csv")

  if (!all(file.exists(c(cra_path,cra_pat_path,p851_path)))) {
    stop("Required Step80B files are missing.")
  }

  cra <- safe_csv(cra_path)
  cra_pat <- safe_csv(cra_pat_path)
  p851 <- safe_csv(p851_path)

  if (is.null(cra) || is.null(cra_pat) || is.null(p851)) {
    stop("Failed to read Step80B core metadata.")
  }

  ck("INPUTS READ")

  # ==========================================================
  # CRA002354: correct patient-level and sample-level status
  # ==========================================================

  patient_status_col <- "Sepsis/Septic shock (1 = sepsis, 2 = septic shock)"
  sample_status_col <- "Sepsis = 1, septic shock = 2"
  survival_col <- "28 days survival (1 = survive, 2 = dead)"
  infection_col <- "Site of infection ( 0 = other, 1 = lung, 2 = intestinal, 3 = abdominal, 4 = blood, 5 = urinary, 6 = brain, 7 = surgical site)"

  needed <- c(patient_status_col,sample_status_col,survival_col,infection_col)
  if (!all(needed %in% names(cra))) {
    stop(
      paste0(
        "Expected CRA clinical columns not found. Missing: ",
        paste(setdiff(needed,names(cra)),collapse=" | ")
      )
    )
  }

  infection_codes <- c(
    "0"="other","1"="lung","2"="intestinal","3"="abdominal",
    "4"="blood","5"="urinary","6"="brain","7"="surgical_site"
  )

  cra <- cra |>
    mutate(
      baseline_sepsis_status=case_when(
        .data[[patient_status_col]]=="1" ~ "sepsis",
        .data[[patient_status_col]]=="2" ~ "septic_shock",
        TRUE ~ NA_character_
      ),
      sample_sepsis_status=case_when(
        .data[[sample_status_col]]=="1" ~ "sepsis",
        .data[[sample_status_col]]=="2" ~ "septic_shock",
        TRUE ~ NA_character_
      ),
      outcome_28d=case_when(
        .data[[survival_col]]=="1" ~ "survived",
        .data[[survival_col]]=="2" ~ "dead",
        TRUE ~ NA_character_
      ),
      infection_source=unname(
        infection_codes[clean_chr(.data[[infection_col]])]
      ),
      sepsis_status=baseline_sepsis_status
    )

  # Patient-level QC must be deduplicated by patient.
  cra_patient_qc <- cra |>
    distinct(
      patient_id,
      baseline_sepsis_status,
      outcome_28d,
      infection_source
    )

  cra_transitions <- cra |>
    distinct(
      patient_id,
      time_day,
      sample_sepsis_status
    ) |>
    arrange(patient_id,time_day) |>
    group_by(patient_id) |>
    summarise(
      n_timepoints=n(),
      statuses=paste(unique(na.omit(sample_sepsis_status)),collapse=" -> "),
      ever_septic_shock=any(sample_sepsis_status=="septic_shock",na.rm=TRUE),
      .groups="drop"
    )

  # ==========================================================
  # Run mappings
  # ==========================================================

  cra_search <- search_mapping(
    "CRA002354",
    cra$sample_id,
    roots=c(
      file.path(ROOT,"data","CRA002354"),
      file.path(ROOT,"results")
    )
  )

  write_excel_csv(
    cra_search$scores,
    file.path(OUT,"CRA002354_run_mapping_candidates.csv"),
    na=""
  )

  cra$sample_norm <- norm_id(cra$sample_id)

  # Remove old blank Run column.
  if ("run_id" %in% names(cra)) cra <- cra |> select(-run_id)

  cra <- cra |>
    left_join(cra_search$map,by="sample_norm")

  if (!("run_id" %in% names(cra))) cra$run_id <- NA_character_

  ck("CRA CLINICAL FIX AND RUN SEARCH COMPLETE")

  p851_search <- search_mapping(
    "PRJNA851469",
    p851$sample_id,
    roots=c(
      file.path(ROOT,"data","PRJNA851469"),
      file.path(ROOT,"results")
    )
  )

  write_excel_csv(
    p851_search$scores,
    file.path(OUT,"PRJNA851469_run_mapping_candidates.csv"),
    na=""
  )

  p851$sample_norm <- norm_id(p851$sample_id)
  if ("run_id" %in% names(p851)) p851 <- p851 |> select(-run_id)

  p851 <- p851 |>
    left_join(p851_search$map,by="sample_norm")

  if (!("run_id" %in% names(p851))) p851$run_id <- NA_character_

  ck("PRJNA851469 RUN SEARCH COMPLETE")

  # ==========================================================
  # QC
  # ==========================================================

  cra_qc <- tibble(
    metric=c(
      "samples","patients","patients_GE2","patients_GE3",
      "baseline_sepsis_patients","baseline_septic_shock_patients",
      "28d_survived_patients","28d_dead_patients",
      "patients_ever_sampled_in_septic_shock",
      "run_id_nonmissing","run_id_missing"
    ),
    value=c(
      nrow(cra),
      n_distinct(cra$patient_id),
      sum((cra |> distinct(patient_id,sample_id) |> count(patient_id))$n>=2),
      sum((cra |> distinct(patient_id,sample_id) |> count(patient_id))$n>=3),
      sum(cra_patient_qc$baseline_sepsis_status=="sepsis",na.rm=TRUE),
      sum(cra_patient_qc$baseline_sepsis_status=="septic_shock",na.rm=TRUE),
      sum(cra_patient_qc$outcome_28d=="survived",na.rm=TRUE),
      sum(cra_patient_qc$outcome_28d=="dead",na.rm=TRUE),
      sum(cra_transitions$ever_septic_shock,na.rm=TRUE),
      sum(!is.na(cra$run_id)),
      sum(is.na(cra$run_id))
    )
  )

  icu851 <- p851 |> filter(sample_class=="ICU_PATIENT")
  pt851 <- icu851 |> distinct(patient_id,time_day) |> count(patient_id)

  p851_qc <- tibble(
    metric=c(
      "all_samples","healthy_controls","ICU_samples","ICU_patients",
      "Day1","Day3","Day7","patients_GE2","patients_GE3",
      "run_id_nonmissing","run_id_missing"
    ),
    value=c(
      nrow(p851),
      sum(p851$sample_class=="HEALTHY_CONTROL"),
      nrow(icu851),
      n_distinct(icu851$patient_id),
      sum(icu851$time_day==1,na.rm=TRUE),
      sum(icu851$time_day==3,na.rm=TRUE),
      sum(icu851$time_day==7,na.rm=TRUE),
      sum(pt851$n>=2),
      sum(pt851$n>=3),
      sum(!is.na(p851$run_id)),
      sum(is.na(p851$run_id))
    )
  )

  decision <- tibble(
    project=c("CRA002354","PRJNA851469"),
    patient_time_mapping=c("EXPLICIT_CONFIRMED","EXPLICIT_CONFIRMED"),
    clinical_mapping=c(
      "EXPLICIT_PATIENT_LEVEL_PLUS_SAMPLE_LEVEL_STATUS",
      "PUBLIC_SAMPLE_LEVEL_PATIENT_DAY"
    ),
    run_mapping=c(
      paste0(sum(!is.na(cra$run_id)),"/",nrow(cra)),
      paste0(sum(!is.na(p851$run_id)),"/",nrow(p851))
    ),
    recommended_role=c(
      "CORE_SEPSIS_LONGITUDINAL_WITH_OUTCOME_INFECTION_SOURCE",
      "ICU_BACKGROUND_LONGITUDINAL"
    ),
    ready_for_master=c(
      TRUE,
      TRUE
    )
  )

  # ==========================================================
  # Write
  # ==========================================================

  write_excel_csv(
    cra,
    file.path(OUT,"CRA002354_MASTER_READY_CORRECTED.csv"),
    na=""
  )

  write_excel_csv(
    cra_patient_qc,
    file.path(OUT,"CRA002354_PATIENT_LEVEL_FINAL.csv"),
    na=""
  )

  write_excel_csv(
    cra_transitions,
    file.path(OUT,"CRA002354_SAMPLE_STATUS_TRANSITIONS.csv"),
    na=""
  )

  write_excel_csv(
    cra_qc,
    file.path(OUT,"CRA002354_STEP80C_QC.csv"),
    na=""
  )

  write_excel_csv(
    p851,
    file.path(OUT,"PRJNA851469_MASTER_READY_WITH_RUN_SEARCH.csv"),
    na=""
  )

  write_excel_csv(
    p851_qc,
    file.path(OUT,"PRJNA851469_STEP80C_QC.csv"),
    na=""
  )

  write_excel_csv(
    decision,
    file.path(OUT,"V2_step80C_decision.csv"),
    na=""
  )

  ck("STEP80C COMPLETE")

  cat("\n============================================\n")
  cat("STEP80C COMPLETE\n")
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
      paste0("STEP80C FATAL ERROR: ",Sys.time()),
      paste0("Message: ",conditionMessage(e)),
      paste0("Call: ",paste(deparse(conditionCall(e)),collapse=" "))
    )
    writeLines(msg,ERR,useBytes=TRUE)
    ck("STEP80C FAILED")
    message(paste(msg,collapse="\n"))
    quit(save="no",status=1,runLast=FALSE)
  }
)
