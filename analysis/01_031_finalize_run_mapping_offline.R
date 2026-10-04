# ============================================================
# Sepsis V2 - Step 80D
# FINAL technical sample -> Run mapping repair
#
# Targets:
#   CRA002354      : recover 131 Sample -> Run mappings
#   PRJNA851469    : adjudicate remaining 4 unmatched samples
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
IN <- file.path(ROOT, "results", "V2_20C_run_mapping_and_clinical_fix")
OUT <- file.path(ROOT, "results", "V2_20D_final_run_mapping")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

LOG <- file.path(OUT, "_STEP80D_runtime_checkpoints.txt")
ERR <- file.path(OUT, "_STEP80D_FATAL_ERROR.txt")
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

count_runs <- function(x) {
  x <- clean_chr(x)
  sum(
    str_detect(
      x,
      regex("^(SRR|ERR|DRR|CRR)[0-9]+$", ignore_case=TRUE)
    ),
    na.rm=TRUE
  )
}

# ============================================================
# CRA002354
# Search local xlsx with skip=0:6 because the technical
# workbook may have title/metadata rows above the real header.
# ============================================================

recover_cra <- function(cra) {

  project_dir <- file.path(ROOT, "data", "CRA002354")

  xlsx <- list.files(
    project_dir,
    recursive=TRUE,
    full.names=TRUE,
    pattern="\\.xlsx$",
    ignore.case=TRUE
  )

  # Prefer technical workbook; exclude downloaded supplements.
  xlsx <- xlsx[
    !str_detect(
      xlsx,
      regex(
        "_step78A3_US_fast_download|mmc[0-9]+|Table_S1|Table_S2",
        ignore_case=TRUE
      )
    )
  ]

  target <- unique(norm_id(cra$sample_id))
  target <- target[!is.na(target)]

  candidates <- list()
  maps <- list()

  inspect <- function(dat, file, sheet, skip) {

    if (is.null(dat) || nrow(dat)==0 || ncol(dat)==0) return(NULL)

    dat <- as_tibble(dat)

    overlap <- map_int(
      dat,
      function(v) {
        vv <- unique(norm_id(v))
        vv <- vv[!is.na(vv)]
        length(intersect(vv,target))
      }
    )

    runs <- map_int(dat,count_runs)

    best_overlap <- max(overlap,na.rm=TRUE)
    best_runs <- max(runs,na.rm=TRUE)

    if (
      !is.finite(best_overlap) ||
      !is.finite(best_runs) ||
      best_overlap==0 ||
      best_runs==0
    ) return(NULL)

    sample_col <- names(dat)[which.max(overlap)]
    run_col <- names(dat)[which.max(runs)]

    score <- tibble(
      source_file=file,
      source_sheet=sheet,
      skip_rows=skip,
      rows=nrow(dat),
      columns=ncol(dat),
      sample_column=sample_col,
      exact_sample_overlap=best_overlap,
      run_column=run_col,
      run_accession_count=best_runs
    )

    map <- tibble(
      sample_norm=norm_id(dat[[sample_col]]),
      run_id_new=clean_chr(dat[[run_col]]),
      run_mapping_source_file=file,
      run_mapping_source_sheet=sheet,
      run_mapping_skip_rows=skip
    ) |>
      filter(
        !is.na(sample_norm),
        !is.na(run_id_new),
        str_detect(
          run_id_new,
          regex("^(SRR|ERR|DRR|CRR)[0-9]+$",ignore_case=TRUE)
        )
      ) |>
      distinct()

    list(score=score,map=map)
  }

  for (p in xlsx) {

    shs <- tryCatch(excel_sheets(p),error=function(e) character())

    for (sh in shs) {
      for (skip in 0:6) {

        dat <- tryCatch(
          read_excel(
            p,
            sheet=sh,
            skip=skip,
            .name_repair="unique"
          ),
          error=function(e) NULL
        )

        z <- inspect(dat,p,sh,skip)

        if (!is.null(z)) {
          candidates[[length(candidates)+1]] <- z$score
          maps[[length(maps)+1]] <- z$map
        }
      }
    }
  }

  score_df <- if (!length(candidates)) {
    tibble()
  } else {
    bind_rows(candidates) |>
      arrange(
        desc(exact_sample_overlap),
        desc(run_accession_count),
        skip_rows
      )
  }

  write_excel_csv(
    score_df,
    file.path(OUT,"CRA002354_technical_workbook_candidates.csv"),
    na=""
  )

  best_map <- tibble(
    sample_norm=character(),
    run_id_new=character(),
    run_mapping_source_file=character(),
    run_mapping_source_sheet=character(),
    run_mapping_skip_rows=integer()
  )

  if (nrow(score_df)) {

    b <- score_df[1,]

    idx <- which(
      map_lgl(
        candidates,
        function(s) {
          identical(s$source_file[1],b$source_file[1]) &&
          identical(s$source_sheet[1],b$source_sheet[1]) &&
          identical(s$skip_rows[1],b$skip_rows[1])
        }
      )
    )

    if (length(idx)) best_map <- maps[[idx[1]]]
  }

  cra2 <- cra |>
    mutate(sample_norm=norm_id(sample_id)) |>
    select(
      -any_of(c(
        "run_id",
        "run_mapping_source_file",
        "run_mapping_source_sheet"
      ))
    ) |>
    left_join(best_map,by="sample_norm") |>
    rename(run_id=run_id_new)

  if (!("run_id" %in% names(cra2))) cra2$run_id <- NA_character_

  list(
    metadata=cra2,
    scores=score_df,
    best_map=best_map
  )
}

# ============================================================
# PRJNA851469
# Use the exact local ENA read_run table.
#
# Auto-map:
# - exact normalized value matches only
# - unique containment matches only
#
# Fuzzy/nearest matches are OUTPUT FOR REVIEW ONLY and are
# NEVER automatically accepted.
# ============================================================

recover_851 <- function(p851) {

  ena_candidates <- list.files(
    file.path(ROOT,"data","PRJNA851469"),
    recursive=TRUE,
    full.names=TRUE,
    pattern="\\.(tsv|csv)$",
    ignore.case=TRUE
  )

  ena_candidates <- ena_candidates[
    str_detect(
      ena_candidates,
      regex("ENA.*read_run|read_run.*ENA",ignore_case=TRUE)
    )
  ]

  if (!length(ena_candidates)) {
    return(list(
      metadata=p851,
      audit=tibble(),
      unmatched_ena=tibble()
    ))
  }

  # Prefer Step75 ENA table explicitly.
  pref <- ena_candidates[
    str_detect(
      ena_candidates,
      regex("step75",ignore_case=TRUE)
    )
  ]

  ena_path <- if (length(pref)) pref[1] else ena_candidates[1]

  ena <- if (tolower(tools::file_ext(ena_path))=="tsv") {
    safe_tsv(ena_path)
  } else {
    safe_csv(ena_path)
  }

  if (is.null(ena)) stop("Failed reading PRJNA851469 ENA read_run table.")

  run_col <- names(ena)[
    which.max(map_int(ena,count_runs))
  ]

  if (count_runs(ena[[run_col]])==0) stop("No Run accession column in ENA table.")

  # Build long index of every textual ENA field.
  ena_index <- bind_rows(
    lapply(
      names(ena),
      function(col) {
        tibble(
          ena_row=seq_len(nrow(ena)),
          field=col,
          raw_value=clean_chr(ena[[col]]),
          norm_value=norm_id(ena[[col]]),
          run_id_ena=clean_chr(ena[[run_col]])
        )
      }
    )
  ) |>
    filter(
      !is.na(norm_value),
      !is.na(run_id_ena)
    )

  p <- p851 |>
    mutate(sample_norm=norm_id(sample_id))

  already_runs <- clean_chr(p$run_id)

  # Exact all-field match first.
  exact <- ena_index |>
    inner_join(
      p |>
        select(sample_id,sample_norm),
      by=c("norm_value"="sample_norm")
    ) |>
    distinct(sample_id,run_id_ena,field,raw_value)

  exact_unique <- exact |>
    group_by(sample_id) |>
    summarise(
      n_runs=n_distinct(run_id_ena),
      run_id_exact=ifelse(n_runs==1,first(run_id_ena),NA_character_),
      exact_fields=paste(unique(field),collapse=" | "),
      .groups="drop"
    )

  p <- p |>
    left_join(exact_unique,by="sample_id") |>
    mutate(
      run_id=coalesce(clean_chr(run_id),run_id_exact),
      mapping_method=case_when(
        !is.na(run_id_exact) ~ "EXACT_ANY_ENA_FIELD",
        !is.na(run_id) ~ "EXISTING_STEP80C_EXACT",
        TRUE ~ NA_character_
      )
    )

  # For remaining unmatched samples, search unique containment.
  unmatched_ids <- p |>
    filter(is.na(run_id)) |>
    pull(sample_id)

  containment_rows <- list()

  for (sid in unmatched_ids) {

    sn <- norm_id(sid)

    hits <- ena_index |>
      filter(
        str_detect(norm_value,fixed(sn)) |
        str_detect(sn,fixed(norm_value))
      ) |>
      distinct(run_id_ena,field,raw_value,norm_value)

    if (nrow(hits)) {
      containment_rows[[length(containment_rows)+1]] <- hits |>
        mutate(sample_id=sid,.before=1)
    }
  }

  containment <- if (!length(containment_rows)) {
    tibble()
  } else {
    bind_rows(containment_rows)
  }

  containment_unique <- if (!nrow(containment)) {
    tibble(
      sample_id=character(),
      n_runs=integer(),
      run_id_containment=character()
    )
  } else {
    containment |>
      group_by(sample_id) |>
      summarise(
        n_runs=n_distinct(run_id_ena),
        run_id_containment=ifelse(n_runs==1,first(run_id_ena),NA_character_),
        containment_evidence=paste(
          unique(paste(field,raw_value,sep="=")),
          collapse=" | "
        ),
        .groups="drop"
      )
  }

  p <- p |>
    left_join(containment_unique,by="sample_id") |>
    mutate(
      run_id=coalesce(run_id,run_id_containment),
      mapping_method=case_when(
        !is.na(mapping_method) ~ mapping_method,
        !is.na(run_id_containment) ~ "UNIQUE_CONTAINMENT_ENA_FIELD",
        TRUE ~ NA_character_
      )
    )

  # Remaining unmatched: nearest candidates for human review only.
  remaining <- p |>
    filter(is.na(run_id)) |>
    select(sample_id,sample_norm)

  audit <- list()

  candidate_strings <- ena_index |>
    filter(
      field %in% c(
        "sample_description",
        "sample_alias",
        "sample_title",
        "sample_accession",
        "secondary_sample_accession"
      )
    ) |>
    distinct(run_id_ena,field,raw_value,norm_value)

  if (!nrow(candidate_strings)) {
    candidate_strings <- ena_index |>
      distinct(run_id_ena,field,raw_value,norm_value)
  }

  for (i in seq_len(nrow(remaining))) {

    sid <- remaining$sample_id[i]
    sn <- remaining$sample_norm[i]

    d <- adist(sn,candidate_strings$norm_value,ignore.case=TRUE)
    ord <- order(d)[seq_len(min(10,length(d)))]

    audit[[length(audit)+1]] <- candidate_strings[ord,] |>
      mutate(
        sample_id=sid,
        edit_distance=as.integer(d[ord]),
        .before=1
      )
  }

  audit_df <- if (!length(audit)) tibble() else bind_rows(audit)

  matched_runs <- unique(clean_chr(p$run_id))
  matched_runs <- matched_runs[!is.na(matched_runs)]

  unmatched_ena <- ena |>
    filter(
      !(.data[[run_col]] %in% matched_runs)
    )

  list(
    metadata=p,
    audit=audit_df,
    unmatched_ena=unmatched_ena,
    ena_path=ena_path
  )
}

main <- function() {

  ck("STEP80D STARTED")

  cra_path <- file.path(IN,"CRA002354_MASTER_READY_CORRECTED.csv")
  p851_path <- file.path(IN,"PRJNA851469_MASTER_READY_WITH_RUN_SEARCH.csv")

  if (!file.exists(cra_path) || !file.exists(p851_path)) {
    stop("Step80C master-ready files are missing.")
  }

  cra <- safe_csv(cra_path)
  p851 <- safe_csv(p851_path)

  if (is.null(cra) || is.null(p851)) stop("Failed to read Step80C files.")

  # CRA
  rcra <- recover_cra(cra)

  cra_final <- rcra$metadata

  cra_qc <- tibble(
    metric=c(
      "samples",
      "run_id_nonmissing",
      "run_id_missing",
      "unique_runs",
      "duplicate_run_assignments"
    ),
    value=c(
      nrow(cra_final),
      sum(!is.na(cra_final$run_id)),
      sum(is.na(cra_final$run_id)),
      n_distinct(cra_final$run_id,na.rm=TRUE),
      sum(duplicated(na.omit(cra_final$run_id)))
    )
  )

  write_excel_csv(
    cra_final,
    file.path(OUT,"CRA002354_MASTER_READY_FINAL.csv"),
    na=""
  )

  write_excel_csv(
    cra_qc,
    file.path(OUT,"CRA002354_FINAL_RUN_QC.csv"),
    na=""
  )

  ck("CRA002354 TECHNICAL MAPPING COMPLETE")

  # PRJNA851469
  r851 <- recover_851(p851)
  p851_final <- r851$metadata

  write_excel_csv(
    r851$audit,
    file.path(OUT,"PRJNA851469_UNMATCHED_nearest_candidates_REVIEW_ONLY.csv"),
    na=""
  )

  write_excel_csv(
    r851$unmatched_ena,
    file.path(OUT,"PRJNA851469_UNMATCHED_ENA_runs.csv"),
    na=""
  )

  p851_qc <- tibble(
    metric=c(
      "samples",
      "run_id_nonmissing",
      "run_id_missing",
      "unique_runs",
      "duplicate_run_assignments"
    ),
    value=c(
      nrow(p851_final),
      sum(!is.na(p851_final$run_id)),
      sum(is.na(p851_final$run_id)),
      n_distinct(p851_final$run_id,na.rm=TRUE),
      sum(duplicated(na.omit(p851_final$run_id)))
    )
  )

  write_excel_csv(
    p851_final,
    file.path(OUT,"PRJNA851469_MASTER_READY_FINAL.csv"),
    na=""
  )

  write_excel_csv(
    p851_qc,
    file.path(OUT,"PRJNA851469_FINAL_RUN_QC.csv"),
    na=""
  )

  ck("PRJNA851469 TECHNICAL MAPPING COMPLETE")

  decision <- tibble(
    project=c("CRA002354","PRJNA851469"),
    total_samples=c(nrow(cra_final),nrow(p851_final)),
    mapped_runs=c(
      sum(!is.na(cra_final$run_id)),
      sum(!is.na(p851_final$run_id))
    ),
    missing_runs=c(
      sum(is.na(cra_final$run_id)),
      sum(is.na(p851_final$run_id))
    ),
    technical_mapping_status=c(
      ifelse(all(!is.na(cra_final$run_id)),"COMPLETE","PARTIAL"),
      ifelse(all(!is.na(p851_final$run_id)),"COMPLETE","PARTIAL")
    ),
    ready_for_step81=c(
      all(!is.na(cra_final$run_id)),
      sum(is.na(p851_final$run_id)) <= 4
    )
  )

  write_excel_csv(
    decision,
    file.path(OUT,"V2_step80D_decision.csv"),
    na=""
  )

  ck("STEP80D COMPLETE")

  cat("\n============================================\n")
  cat("STEP80D COMPLETE\n")
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
      paste0("STEP80D FATAL ERROR: ",Sys.time()),
      paste0("Message: ",conditionMessage(e)),
      paste0("Call: ",paste(deparse(conditionCall(e)),collapse=" "))
    )
    writeLines(msg,ERR,useBytes=TRUE)
    ck("STEP80D FAILED")
    message(paste(msg,collapse="\n"))
    quit(save="no",status=1,runLast=FALSE)
  }
)
