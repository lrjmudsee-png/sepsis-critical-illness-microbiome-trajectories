# ============================================================
# Sepsis V2 - Step 80E
# Strict final technical mapping
#
# 1) CRA002354:
#    Search the ENTIRE project + user's Downloads for the
#    official CRA002354.xlsx / GSA metadata workbook.
#    Try skip=0:10 and map supplement Sample ID -> CRR exactly.
#
# 2) PRJNA851469:
#    Roll back unsafe substring mappings from Step80D.
#    Keep ONLY exact Step80C mappings.
#    Document 4 supplement-only samples + ENA-only runs.
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
IN80C <- file.path(ROOT,"results","V2_20C_run_mapping_and_clinical_fix")
OUT <- file.path(ROOT,"results","V2_20E_strict_final_technical_mapping")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

LOG <- file.path(OUT,"_STEP80E_runtime_checkpoints.txt")
ERR <- file.path(OUT,"_STEP80E_FATAL_ERROR.txt")
if (file.exists(ERR)) unlink(ERR)

ck <- function(x) cat(paste0(x,": ",Sys.time(),"\n"),file=LOG,append=TRUE)

clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)|x==""|tolower(x)%in%c("na","nan","n/a","null","none")] <- NA_character_
  x
}

norm_id <- function(x) {
  x <- toupper(clean_chr(x))
  str_replace_all(x,"[^A-Z0-9]","")
}

safe_csv <- function(p) {
  tryCatch(
    suppressMessages(read_csv(p,show_col_types=FALSE,progress=FALSE)),
    error=function(e) NULL
  )
}

run_n <- function(x) {
  x <- clean_chr(x)
  sum(
    str_detect(
      x,
      regex("^(CRR|SRR|ERR|DRR)[0-9]+$",ignore_case=TRUE)
    ),
    na.rm=TRUE
  )
}

exp_n <- function(x) {
  x <- clean_chr(x)
  sum(
    str_detect(
      x,
      regex("^(CRX|SRX|ERX|DRX)[0-9]+$",ignore_case=TRUE)
    ),
    na.rm=TRUE
  )
}

sample_acc_n <- function(x) {
  x <- clean_chr(x)
  sum(
    str_detect(
      x,
      regex("^(SAMC|SAMN|SAMEA|SAMD|SRS|ERS|DRS)[A-Z0-9]+$",ignore_case=TRUE)
    ),
    na.rm=TRUE
  )
}

main <- function() {

  ck("STEP80E STARTED")

  cra_path <- file.path(IN80C,"CRA002354_MASTER_READY_CORRECTED.csv")
  p851_path <- file.path(IN80C,"PRJNA851469_MASTER_READY_WITH_RUN_SEARCH.csv")

  if (!file.exists(cra_path) || !file.exists(p851_path)) {
    stop("Required Step80C master-ready files are missing.")
  }

  cra <- safe_csv(cra_path)
  p851 <- safe_csv(p851_path)

  if (is.null(cra) || is.null(p851)) stop("Failed to read Step80C files.")

  # ==========================================================
  # A. CRA002354 - exhaustive workbook search
  # ==========================================================

  userprofile <- Sys.getenv("USERPROFILE")

  roots <- unique(c(
    ROOT,
    file.path(userprofile,"Downloads"),
    file.path(userprofile,"Desktop")
  ))

  xlsx <- unique(unlist(lapply(roots,function(r) {
    if (!dir.exists(r)) return(character())
    list.files(
      r,
      recursive=TRUE,
      full.names=TRUE,
      pattern="\\.xlsx$",
      ignore.case=TRUE
    )
  })))

  xlsx <- xlsx[file.exists(xlsx)]

  # Prefer filenames explicitly referring to CRA002354.
  named <- xlsx[
    str_detect(
      basename(xlsx),
      regex("CRA002354",ignore_case=TRUE)
    )
  ]

  # Scan named candidates first; if none, scan reasonably small xlsx
  # under the project's CRA-related paths/results.
  if (length(named)) {
    scan_files <- named
  } else {
    scan_files <- xlsx[
      file.info(xlsx)$size < 10*1024^2 &
      (
        str_detect(xlsx,regex("CRA002354|sepsis",ignore_case=TRUE))
      )
    ]
  }

  target_norm <- unique(norm_id(cra$sample_id))
  target_norm <- target_norm[!is.na(target_norm)]

  candidate_scores <- list()
  candidate_maps <- list()
  technical_exports <- list()

  inspect <- function(dat,file,sheet,skip) {

    if (is.null(dat) || nrow(dat)==0 || ncol(dat)==0) return(NULL)
    dat <- as_tibble(dat)

    overlaps <- map_int(dat,function(v) {
      vv <- unique(norm_id(v))
      vv <- vv[!is.na(vv)]
      length(intersect(vv,target_norm))
    })

    runs <- map_int(dat,run_n)
    exps <- map_int(dat,exp_n)
    sacc <- map_int(dat,sample_acc_n)

    best_overlap <- max(overlaps,na.rm=TRUE)
    best_runs <- max(runs,na.rm=TRUE)

    # Keep candidate if it either has direct supplement overlap,
    # OR has a strong technical GSA structure (~131 Runs).
    technical_like <- best_runs >= 100

    if (
      (!is.finite(best_overlap) || best_overlap==0) &&
      !technical_like
    ) return(NULL)

    sample_col <- if (best_overlap>0) names(dat)[which.max(overlaps)] else NA_character_
    run_col <- if (best_runs>0) names(dat)[which.max(runs)] else NA_character_
    exp_col <- if (max(exps,na.rm=TRUE)>0) names(dat)[which.max(exps)] else NA_character_
    sample_acc_col <- if (max(sacc,na.rm=TRUE)>0) names(dat)[which.max(sacc)] else NA_character_

    score <- tibble(
      source_file=file,
      source_sheet=sheet,
      skip_rows=skip,
      rows=nrow(dat),
      columns=ncol(dat),
      supplement_sample_column=sample_col,
      exact_supplement_overlap=best_overlap,
      run_column=run_col,
      run_accession_count=best_runs,
      experiment_column=exp_col,
      experiment_accession_count=max(exps,na.rm=TRUE),
      sample_accession_column=sample_acc_col,
      sample_accession_count=max(sacc,na.rm=TRUE)
    )

    map <- tibble()

    if (!is.na(sample_col) && !is.na(run_col)) {
      map <- tibble(
        sample_norm=norm_id(dat[[sample_col]]),
        run_id=clean_chr(dat[[run_col]]),
        experiment_id=if (!is.na(exp_col)) clean_chr(dat[[exp_col]]) else NA_character_,
        archive_sample_accession=if (!is.na(sample_acc_col)) clean_chr(dat[[sample_acc_col]]) else NA_character_,
        source_file=file,
        source_sheet=sheet,
        skip_rows=skip
      ) |>
        filter(
          !is.na(sample_norm),
          !is.na(run_id),
          str_detect(run_id,regex("^CRR[0-9]+$",ignore_case=TRUE))
        ) |>
        distinct()
    }

    # Export technical-looking table even if supplement IDs don't overlap,
    # so the exact obstacle is visible.
    if (technical_like) {
      keycols <- unique(na.omit(c(sample_col,sample_acc_col,exp_col,run_col)))
      if (length(keycols)) {
        technical <- dat |>
          select(all_of(keycols))
      } else {
        technical <- dat
      }
    } else {
      technical <- NULL
    }

    list(score=score,map=map,technical=technical)
  }

  for (f in scan_files) {

    if (file.info(f)$size > 20*1024^2) next

    sheets <- tryCatch(excel_sheets(f),error=function(e) character())

    for (sh in sheets) {
      for (skip in 0:10) {

        dat <- tryCatch(
          read_excel(
            f,
            sheet=sh,
            skip=skip,
            .name_repair="unique"
          ),
          error=function(e) NULL
        )

        z <- inspect(dat,f,sh,skip)

        if (!is.null(z)) {
          candidate_scores[[length(candidate_scores)+1]] <- z$score
          candidate_maps[[length(candidate_maps)+1]] <- z$map

          if (!is.null(z$technical)) {
            technical_exports[[length(technical_exports)+1]] <- list(
              score=z$score,
              data=z$technical
            )
          }
        }
      }
    }
  }

  score_df <- if (!length(candidate_scores)) {
    tibble()
  } else {
    bind_rows(candidate_scores) |>
      arrange(
        desc(exact_supplement_overlap),
        desc(run_accession_count),
        desc(experiment_accession_count),
        skip_rows
      )
  }

  write_excel_csv(
    score_df,
    file.path(OUT,"CRA002354_EXHAUSTIVE_workbook_search.csv"),
    na=""
  )

  cra_map <- tibble(
    sample_norm=character(),
    run_id=character(),
    experiment_id=character(),
    archive_sample_accession=character(),
    source_file=character(),
    source_sheet=character(),
    skip_rows=integer()
  )

  if (nrow(score_df)>0 && score_df$exact_supplement_overlap[1] > 0) {

    b <- score_df[1,]

    idx <- which(map_lgl(candidate_scores,function(s) {
      identical(s$source_file[1],b$source_file[1]) &&
      identical(s$source_sheet[1],b$source_sheet[1]) &&
      identical(s$skip_rows[1],b$skip_rows[1])
    }))

    if (length(idx)) cra_map <- candidate_maps[[idx[1]]]
  }

  # If no direct map, export the best 131-Run technical table for review.
  if (nrow(cra_map)==0 && length(technical_exports)) {

    t_scores <- bind_rows(lapply(technical_exports,`[[`,"score")) |>
      arrange(
        desc(run_accession_count),
        desc(experiment_accession_count),
        desc(sample_accession_count)
      )

    best <- t_scores[1,]

    idx <- which(map_lgl(technical_exports,function(z) {
      s <- z$score
      identical(s$source_file[1],best$source_file[1]) &&
      identical(s$source_sheet[1],best$source_sheet[1]) &&
      identical(s$skip_rows[1],best$skip_rows[1])
    }))

    if (length(idx)) {
      write_excel_csv(
        technical_exports[[idx[1]]]$data,
        file.path(OUT,"CRA002354_BEST_GSA_TECHNICAL_TABLE_REVIEW.csv"),
        na=""
      )
    }
  }

  cra_final <- cra |>
    mutate(sample_norm=norm_id(sample_id)) |>
    select(
      -any_of(c(
        "run_id",
        "run_mapping_source_file",
        "run_mapping_source_sheet"
      ))
    ) |>
    left_join(cra_map,by="sample_norm")

  if (!("run_id" %in% names(cra_final))) cra_final$run_id <- NA_character_

  write_excel_csv(
    cra_final,
    file.path(OUT,"CRA002354_MASTER_READY_STRICT_FINAL.csv"),
    na=""
  )

  cra_qc <- tibble(
    metric=c(
      "samples",
      "mapped_runs",
      "missing_runs",
      "unique_runs",
      "candidate_workbooks_found",
      "best_exact_sample_overlap"
    ),
    value=c(
      nrow(cra_final),
      sum(!is.na(cra_final$run_id)),
      sum(is.na(cra_final$run_id)),
      n_distinct(cra_final$run_id,na.rm=TRUE),
      length(scan_files),
      ifelse(nrow(score_df),score_df$exact_supplement_overlap[1],0)
    )
  )

  write_excel_csv(
    cra_qc,
    file.path(OUT,"CRA002354_STEP80E_QC.csv"),
    na=""
  )

  ck("CRA002354 EXHAUSTIVE SEARCH COMPLETE")

  # ==========================================================
  # B. PRJNA851469 - STRICT rollback to Step80C exact mapping
  # ==========================================================

  # Step80C mapping was exact-match based. Preserve it as truth.
  p851_final <- p851 |>
    mutate(
      run_mapping_status=case_when(
        !is.na(run_id) ~ "EXACT_CONFIRMED",
        TRUE ~ "UNAVAILABLE_OR_UNMAPPED"
      )
    )

  # Never use substring containment such as HealthyVolunteer-1
  # matching HealthyVolunteer-14/15.
  strict_missing <- p851_final |>
    filter(is.na(run_id)) |>
    select(
      sample_id,
      sample_class,
      patient_id,
      time_raw,
      time_day,
      run_mapping_status
    )

  write_excel_csv(
    strict_missing,
    file.path(OUT,"PRJNA851469_STRICT_UNMATCHED_4_SAMPLES.csv"),
    na=""
  )

  # Identify sequence-analysis availability.
  p851_final <- p851_final |>
    mutate(
      sequence_available=!is.na(run_id),
      sequence_analysis_include=!is.na(run_id)
    )

  write_excel_csv(
    p851_final,
    file.path(OUT,"PRJNA851469_MASTER_READY_STRICT_FINAL.csv"),
    na=""
  )

  icu_seq <- p851_final |>
    filter(
      sample_class=="ICU_PATIENT",
      sequence_analysis_include
    )

  pt_seq <- icu_seq |>
    distinct(patient_id,time_day) |>
    count(patient_id,name="n_timepoints")

  q851 <- tibble(
    metric=c(
      "supplement_samples",
      "exact_mapped_runs",
      "unmapped_samples",
      "unique_mapped_runs",
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
      sum(pt_seq$n_timepoints>=2),
      sum(pt_seq$n_timepoints>=3)
    )
  )

  write_excel_csv(
    q851,
    file.path(OUT,"PRJNA851469_STEP80E_QC.csv"),
    na=""
  )

  ck("PRJNA851469 STRICT MAPPING COMPLETE")

  # ==========================================================
  # Decision
  # ==========================================================

  decision <- tibble(
    project=c("CRA002354","PRJNA851469"),
    patient_time_metadata_ready=c(TRUE,TRUE),
    run_mapping=c(
      paste0(sum(!is.na(cra_final$run_id)),"/",nrow(cra_final)),
      paste0(sum(!is.na(p851_final$run_id)),"/",nrow(p851_final))
    ),
    sequence_analysis_policy=c(
      ifelse(
        all(!is.na(cra_final$run_id)),
        "READY",
        "HOLD_SEQUENCE_LAYER_UNTIL_GSA_MAP_RESOLVED"
      ),
      "USE_119_EXACT_MAPPED_SAMPLES_ONLY"
    ),
    ready_for_master_metadata=c(TRUE,TRUE)
  )

  write_excel_csv(
    decision,
    file.path(OUT,"V2_step80E_decision.csv"),
    na=""
  )

  ck("STEP80E COMPLETE")

  cat("\n============================================\n")
  cat("STEP80E COMPLETE\n")
  cat("============================================\n")
  cat("\nCRA002354:\n")
  print(cra_qc,n=Inf)
  cat("\nPRJNA851469:\n")
  print(q851,n=Inf)
  cat("\nDecision:\n")
  print(decision,n=Inf)
  cat("\nOutput: ",OUT,"\n",sep="")
}

tryCatch(
  main(),
  error=function(e) {
    msg <- c(
      paste0("STEP80E FATAL ERROR: ",Sys.time()),
      paste0("Message: ",conditionMessage(e)),
      paste0("Call: ",paste(deparse(conditionCall(e)),collapse=" "))
    )
    writeLines(msg,ERR,useBytes=TRUE)
    ck("STEP80E FAILED")
    message(paste(msg,collapse="\n"))
    quit(save="no",status=1,runLast=FALSE)
  }
)
