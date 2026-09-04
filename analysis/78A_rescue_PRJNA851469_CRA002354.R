# ============================================================
# Sepsis V2 - Step 78A
# Rescue PRJNA851469 and CRA002354 public metadata
# R 4.4.0 / Windows
#
# Independent project sections:
# failure in one project does NOT stop the other.
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr","readxl","dplyr","tidyr","stringr",
  "purrr","tibble","httr2"
)

missing <- pkgs[
  !vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)
]

if (length(missing) > 0) {
  install.packages(
    missing,
    repos="https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(httr2)
})

PROJECT_ROOT <- "E:/sepsis_project"
DATA_ROOT <- file.path(PROJECT_ROOT, "data")
OUT_ROOT <- file.path(
  PROJECT_ROOT,
  "results",
  "V2_18A_851469_CRA002354_rescue"
)

dir.create(
  OUT_ROOT,
  recursive=TRUE,
  showWarnings=FALSE
)

LOG <- file.path(
  OUT_ROOT,
  "_STEP78A_log.txt"
)

cat(
  paste0("STEP78A START: ", Sys.time(), "\n"),
  file=LOG,
  append=TRUE
)

UA <- "SepsisV2-Step78A-R44/1.0"

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[
    is.na(x) |
    x == "" |
    tolower(x) %in% c("na","nan","n/a","null","none")
  ] <- NA_character_
  x
}

safe_download <- function(url, dest) {

  dir.create(
    dirname(dest),
    recursive=TRUE,
    showWarnings=FALSE
  )

  if (
    file.exists(dest) &&
    file.info(dest)$size > 1000
  ) {
    return("EXISTS")
  }

  tryCatch(
    {
      resp <- request(url) |>
        req_user_agent(UA) |>
        req_timeout(300) |>
        req_perform(path=dest)

      paste0("HTTP_", resp_status(resp))
    },
    error=function(e) {
      paste0("FAILED:", conditionMessage(e))
    }
  )
}

profile_excel <- function(path, project) {

  if (
    !file.exists(path) ||
    file.info(path)$size < 1000
  ) return(tibble())

  sheets <- tryCatch(
    excel_sheets(path),
    error=function(e) character()
  )

  out <- list()

  for (sh in sheets) {

    x <- tryCatch(
      read_excel(
        path,
        sheet=sh,
        n_max=10000,
        .name_repair="unique"
      ),
      error=function(e) NULL
    )

    if (is.null(x)) next

    nms <- names(x)

    out[[length(out)+1]] <- tibble(
      project=project,
      file=basename(path),
      sheet=sh,
      rows=nrow(x),
      columns=ncol(x),
      column_names=paste(nms, collapse=";"),

      patient_field=any(
        str_detect(
          nms,
          regex(
            "patient|subject|participant|individual",
            ignore_case=TRUE
          )
        )
      ),

      sample_field=any(
        str_detect(
          nms,
          regex(
            "sample|biosample|specimen|swab|run",
            ignore_case=TRUE
          )
        )
      ),

      time_field=any(
        str_detect(
          nms,
          regex(
            "day|time|visit|admission|discharge",
            ignore_case=TRUE
          )
        )
      ),

      infection_field=any(
        str_detect(
          nms,
          regex(
            "infection|nosocomial|sepsis|culture|pathogen",
            ignore_case=TRUE
          )
        )
      ),

      outcome_field=any(
        str_detect(
          nms,
          regex(
            "mortality|death|survival|outcome",
            ignore_case=TRUE
          )
        )
      ),

      antibiotics_field=any(
        str_detect(
          nms,
          regex(
            "antibiotic|abx|antimicrobial",
            ignore_case=TRUE
          )
        )
      )
    )
  }

  if (length(out)==0) tibble() else bind_rows(out)
}

export_sheet_previews <- function(path, project, prefix) {

  if (
    !file.exists(path) ||
    file.info(path)$size < 1000
  ) return(invisible(NULL))

  sheets <- tryCatch(
    excel_sheets(path),
    error=function(e) character()
  )

  for (i in seq_along(sheets)) {

    sh <- sheets[i]

    x <- tryCatch(
      read_excel(
        path,
        sheet=sh,
        n_max=10000,
        .name_repair="unique"
      ),
      error=function(e) NULL
    )

    if (is.null(x)) next

    safe_sh <- str_replace_all(
      sh,
      "[^A-Za-z0-9]+",
      "_"
    )

    write_excel_csv(
      x,
      file.path(
        OUT_ROOT,
        paste0(
          prefix,
          "_sheet_",
          i,
          "_",
          safe_sh,
          ".csv"
        )
      ),
      na=""
    )
  }
}

# ============================================================
# A. PRJNA851469
# ============================================================
status851 <- tryCatch({

  project <- "PRJNA851469"

  dir851 <- file.path(
    DATA_ROOT,
    project,
    "00_metadata",
    "step78A_public_supplement"
  )

  dir.create(
    dir851,
    recursive=TRUE,
    showWarnings=FALSE
  )

  # Official Nature Medicine supplementary workbook,
  # Supplementary Tables 2-17.
  url851 <- paste0(
    "https://media.springernature.com/original/",
    "springer-static/esm/art%3A10.1038%2Fs41591-023-02243-5/",
    "MediaObjects/41591_2023_2243_MOESM3_ESM.xlsx"
  )

  path851 <- file.path(
    dir851,
    "PRJNA851469_Supplementary_Tables_2_17.xlsx"
  )

  st <- safe_download(
    url851,
    path851
  )

  prof <- profile_excel(
    path851,
    project
  )

  write_excel_csv(
    prof,
    file.path(
      OUT_ROOT,
      "PRJNA851469_supplement_profile.csv"
    ),
    na=""
  )

  export_sheet_previews(
    path851,
    project,
    "PRJNA851469"
  )

  # Search values, not just headers, for public identifiers.
  candidate_rows <- list()

  if (
    file.exists(path851) &&
    file.info(path851)$size > 1000
  ) {

    sheets <- excel_sheets(path851)

    for (sh in sheets) {

      x <- tryCatch(
        read_excel(
          path851,
          sheet=sh,
          n_max=10000,
          .name_repair="unique"
        ),
        error=function(e) NULL
      )

      if (is.null(x)) next

      row_text <- apply(
        as.data.frame(x),
        1,
        function(z) {
          paste(
            clean_chr(z),
            collapse=" | "
          )
        }
      )

      hit <- str_detect(
        row_text,
        regex(
          "patient|subject|day ?1|day ?3|day ?7|infection|nosocomial|sample",
          ignore_case=TRUE
        )
      )

      if (any(hit, na.rm=TRUE)) {
        candidate_rows[[length(candidate_rows)+1]] <-
          x[which(hit), , drop=FALSE] |>
          mutate(
            Source_Sheet=sh,
            .before=1
          )
      }
    }
  }

  cand <- if (length(candidate_rows)>0) {
    bind_rows(candidate_rows)
  } else {
    tibble()
  }

  if (ncol(cand)>0) {
    write_excel_csv(
      cand,
      file.path(
        OUT_ROOT,
        "PRJNA851469_candidate_public_metadata_rows.csv"
      ),
      na=""
    )
  }

  cat(
    paste0(
      "PRJNA851469 COMPLETE: ",
      st,
      " ",
      Sys.time(),
      "\n"
    ),
    file=LOG,
    append=TRUE
  )

  tibble(
    Project=project,
    Download_Status=st,
    Workbook_Exists=file.exists(path851),
    Workbook_Bytes=if(file.exists(path851)) file.info(path851)$size else 0,
    Sheets_Profiled=nrow(prof),
    Status=if(nrow(prof)>0) "PUBLIC_SUPPLEMENT_RECOVERED" else "SUPPLEMENT_RECOVERY_FAILED"
  )

}, error=function(e) {

  cat(
    paste0(
      "PRJNA851469 ERROR: ",
      conditionMessage(e),
      "\n"
    ),
    file=LOG,
    append=TRUE
  )

  tibble(
    Project="PRJNA851469",
    Download_Status="ERROR",
    Workbook_Exists=FALSE,
    Workbook_Bytes=0,
    Sheets_Profiled=0,
    Status=paste0("ERROR:", conditionMessage(e))
  )
})

# ============================================================
# B. CRA002354
# ============================================================
statusCRA <- tryCatch({

  project <- "CRA002354"

  dirCRA <- file.path(
    DATA_ROOT,
    project,
    "00_metadata",
    "step78A_public_supplement"
  )

  dir.create(
    dirCRA,
    recursive=TRUE,
    showWarnings=FALSE
  )

  # Official supplementary tables described by the PMC article.
  resources <- tibble(
    Resource=c(
      "Supplementary_Table_S1_64_patients",
      "Supplementary_Table_S2_131_samples"
    ),

    File=c(
      "CRA002354_Table_S1_64_patients.xlsx",
      "CRA002354_Table_S2_131_samples.xlsx"
    ),

    URL1=c(
      "https://pmc.ncbi.nlm.nih.gov/articles/PMC8377022/bin/mmc5.xlsx",
      "https://pmc.ncbi.nlm.nih.gov/articles/PMC8377022/bin/mmc6.xlsx"
    ),

    URL2=c(
      "https://pmc.ncbi.nlm.nih.gov/articles/instance/8377022/bin/mmc5.xlsx",
      "https://pmc.ncbi.nlm.nih.gov/articles/instance/8377022/bin/mmc6.xlsx"
    )
  )

  logs <- list()
  profiles <- list()

  for (i in seq_len(nrow(resources))) {

    dest <- file.path(
      dirCRA,
      resources$File[i]
    )

    st <- safe_download(
      resources$URL1[i],
      dest
    )

    if (
      !file.exists(dest) ||
      file.info(dest)$size < 1000
    ) {
      st <- safe_download(
        resources$URL2[i],
        dest
      )
    }

    logs[[i]] <- tibble(
      Resource=resources$Resource[i],
      Local_Path=dest,
      Status=st,
      Exists=file.exists(dest),
      Bytes=if(file.exists(dest)) file.info(dest)$size else 0
    )

    p <- profile_excel(
      dest,
      project
    )

    if (nrow(p)>0) {
      profiles[[length(profiles)+1]] <- p
    }

    export_sheet_previews(
      dest,
      project,
      paste0(
        "CRA002354_",
        ifelse(i==1, "S1", "S2")
      )
    )
  }

  logdf <- bind_rows(logs)

  profdf <- if (length(profiles)>0) {
    bind_rows(profiles)
  } else {
    tibble()
  }

  write_excel_csv(
    logdf,
    file.path(
      OUT_ROOT,
      "CRA002354_supplement_download_log.csv"
    ),
    na=""
  )

  write_excel_csv(
    profdf,
    file.path(
      OUT_ROOT,
      "CRA002354_supplement_profile.csv"
    ),
    na=""
  )

  # ----------------------------------------------------------
  # Locate already-local CRA002354 workbook from GSA.
  # ----------------------------------------------------------
  local_xlsx <- list.files(
    file.path(DATA_ROOT, project),
    pattern="\\.xlsx$",
    recursive=TRUE,
    full.names=TRUE,
    ignore.case=TRUE
  )

  local_xlsx <- local_xlsx[
    !str_detect(
      local_xlsx,
      regex(
        "Table_S1|Table_S2",
        ignore_case=TRUE
      )
    )
  ]

  local_profile <- list()

  for (p in local_xlsx) {

    q <- profile_excel(
      p,
      project
    )

    if (nrow(q)>0) {
      local_profile[[length(local_profile)+1]] <- q
    }
  }

  local_profile_df <- if (length(local_profile)>0) {
    bind_rows(local_profile)
  } else {
    tibble()
  }

  write_excel_csv(
    local_profile_df,
    file.path(
      OUT_ROOT,
      "CRA002354_existing_local_workbook_profile.csv"
    ),
    na=""
  )

  cat(
    paste0(
      "CRA002354 COMPLETE: ",
      Sys.time(),
      "\n"
    ),
    file=LOG,
    append=TRUE
  )

  tibble(
    Project=project,
    Download_Status=paste(logdf$Status, collapse=";"),
    Workbook_Exists=all(logdf$Exists),
    Workbook_Bytes=sum(logdf$Bytes),
    Sheets_Profiled=nrow(profdf),
    Status=if(all(logdf$Exists)) "S1_S2_SUPPLEMENTS_RECOVERED" else "PARTIAL_SUPPLEMENT_RECOVERY"
  )

}, error=function(e) {

  cat(
    paste0(
      "CRA002354 ERROR: ",
      conditionMessage(e),
      "\n"
    ),
    file=LOG,
    append=TRUE
  )

  tibble(
    Project="CRA002354",
    Download_Status="ERROR",
    Workbook_Exists=FALSE,
    Workbook_Bytes=0,
    Sheets_Profiled=0,
    Status=paste0("ERROR:", conditionMessage(e))
  )
})

# ============================================================
# Final summary
# ============================================================
summary <- bind_rows(
  status851,
  statusCRA
)

write_excel_csv(
  summary,
  file.path(
    OUT_ROOT,
    "V2_step78A_status.csv"
  ),
  na=""
)

cat(
  paste0(
    "STEP78A COMPLETE: ",
    Sys.time(),
    "\n"
  ),
  file=LOG,
  append=TRUE
)

cat("\n=============================================\n")
cat("SEPSIS V2 STEP 78A COMPLETE\n")
cat("=============================================\n")
print(summary, n=Inf, width=Inf)
cat("\nOutput: ", OUT_ROOT, "\n", sep="")
