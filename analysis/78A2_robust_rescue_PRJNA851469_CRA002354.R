# ============================================================
# Sepsis V2 - Step 78A2
# Robust rescue for PRJNA851469 + CRA002354
# R 4.4.0 / Windows
#
# Key fixes:
# - validate XLSX ZIP magic + readxl readability
# - delete stale/corrupt pseudo-xlsx files
# - PRJNA851469: redownload official Nature MOESM3 XLSX
# - CRA002354: obtain official PMC OA package and extract
#   mmc5.xlsx (Table S1) + mmc6.xlsx (Table S2)
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr","readxl","dplyr","tidyr","stringr",
  "purrr","tibble","httr2","xml2"
)

missing <- pkgs[
  !vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)
]

if (length(missing) > 0) {
  install.packages(missing, repos="https://cloud.r-project.org")
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
  library(xml2)
})

PROJECT_ROOT <- "E:/sepsis_project"
DATA_ROOT <- file.path(PROJECT_ROOT, "data")
OUT_ROOT <- file.path(
  PROJECT_ROOT,
  "results",
  "V2_18A2_robust_supplement_rescue"
)

dir.create(
  OUT_ROOT,
  recursive=TRUE,
  showWarnings=FALSE
)

LOG <- file.path(
  OUT_ROOT,
  "_STEP78A2_log.txt"
)

cat(
  paste0("STEP78A2 START: ", Sys.time(), "\n"),
  file=LOG,
  append=TRUE
)

UA <- "SepsisV2-Step78A2-R44/1.0"

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

is_valid_xlsx <- function(path) {

  if (
    !file.exists(path) ||
    file.info(path)$size < 100
  ) return(FALSE)

  # XLSX is a ZIP archive and should start with PK.
  con <- file(path, "rb")
  on.exit(close(con), add=TRUE)

  magic <- readBin(
    con,
    what="raw",
    n=4
  )

  if (length(magic) < 2) return(FALSE)

  pk_ok <- (
    as.integer(magic[1]) == 0x50 &&
    as.integer(magic[2]) == 0x4B
  )

  if (!pk_ok) return(FALSE)

  # Also require readxl to recognize workbook structure.
  tryCatch(
    {
      sh <- excel_sheets(path)
      length(sh) > 0
    },
    error=function(e) FALSE
  )
}

download_binary <- function(url, dest) {

  dir.create(
    dirname(dest),
    recursive=TRUE,
    showWarnings=FALSE
  )

  if (file.exists(dest)) {
    unlink(dest, force=TRUE)
  }

  tryCatch(
    {
      resp <- request(url) |>
        req_user_agent(UA) |>
        req_headers(
          Accept="application/octet-stream,*/*"
        ) |>
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

  if (!is_valid_xlsx(path)) {
    return(tibble())
  }

  sheets <- excel_sheets(path)
  out <- list()

  for (sh in sheets) {

    x <- tryCatch(
      read_excel(
        path,
        sheet=sh,
        n_max=20000,
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
            "antibiotic|abx|antimicrobial|carbapenem",
            ignore_case=TRUE
          )
        )
      )
    )
  }

  if (length(out)==0) tibble() else bind_rows(out)
}

export_all_sheets <- function(path, prefix) {

  if (!is_valid_xlsx(path)) return(invisible(NULL))

  sheets <- excel_sheets(path)

  for (i in seq_along(sheets)) {

    sh <- sheets[i]

    x <- tryCatch(
      read_excel(
        path,
        sheet=sh,
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
    "step78A2_public_supplement"
  )

  dir.create(
    dir851,
    recursive=TRUE,
    showWarnings=FALSE
  )

  url851 <- paste0(
    "https://media.springernature.com/original/",
    "springer-static/esm/art%3A10.1038%2Fs41591-023-02243-5/",
    "MediaObjects/41591_2023_2243_MOESM3_ESM.xlsx"
  )

  path851 <- file.path(
    dir851,
    "PRJNA851469_Supplementary_Tables_2_17.xlsx"
  )

  st <- download_binary(
    url851,
    path851
  )

  valid <- is_valid_xlsx(
    path851
  )

  if (!valid) {
    stop(
      paste0(
        "Nature file downloaded but failed XLSX validation. Status=",
        st,
        ", bytes=",
        if(file.exists(path851)) file.info(path851)$size else 0
      )
    )
  }

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

  export_all_sheets(
    path851,
    "PRJNA851469"
  )

  cat(
    paste0(
      "PRJNA851469 COMPLETE: ",
      st,
      " bytes=",
      file.info(path851)$size,
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
    Valid_XLSX=valid,
    Bytes=file.info(path851)$size,
    Sheets_Profiled=nrow(prof),
    Status="PUBLIC_SUPPLEMENT_RECOVERED"
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
    Valid_XLSX=FALSE,
    Bytes=0,
    Sheets_Profiled=0,
    Status=paste0("ERROR:", conditionMessage(e))
  )
})

# ============================================================
# B. CRA002354 via official NCBI PMC OA package
# ============================================================
statusCRA <- tryCatch({

  project <- "CRA002354"

  dirCRA <- file.path(
    DATA_ROOT,
    project,
    "00_metadata",
    "step78A2_public_supplement"
  )

  dir.create(
    dirCRA,
    recursive=TRUE,
    showWarnings=FALSE
  )

  oa_url <- paste0(
    "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?",
    "id=PMC8377022"
  )

  oa_xml_path <- file.path(
    dirCRA,
    "PMC8377022_oa.xml"
  )

  oa_status <- download_binary(
    oa_url,
    oa_xml_path
  )

  if (
    !file.exists(oa_xml_path) ||
    file.info(oa_xml_path)$size < 50
  ) {
    stop("Failed to retrieve PMC OA XML.")
  }

  doc <- read_xml(
    oa_xml_path
  )

  links <- xml_find_all(
    doc,
    "//link"
  )

  href <- xml_attr(
    links,
    "href"
  )

  format <- xml_attr(
    links,
    "format"
  )

  tgz <- href[
    format == "tgz"
  ]

  if (length(tgz)==0) {
    tgz <- href[
      str_detect(
        href,
        regex("\\.tar\\.gz$", ignore_case=TRUE)
      )
    ]
  }

  if (length(tgz)==0) {
    stop("PMC OA XML contains no tgz package link.")
  }

  package_url <- tgz[1]

  package_url <- sub(
    "^ftp://ftp\\.ncbi\\.nlm\\.nih\\.gov",
    "https://ftp.ncbi.nlm.nih.gov",
    package_url
  )

  tgz_path <- file.path(
    dirCRA,
    "PMC8377022_oa_package.tar.gz"
  )

  tgz_status <- download_binary(
    package_url,
    tgz_path
  )

  if (
    !file.exists(tgz_path) ||
    file.info(tgz_path)$size < 1000
  ) {
    stop("Failed to download PMC OA package.")
  }

  extract_dir <- file.path(
    dirCRA,
    "PMC8377022_extracted"
  )

  if (dir.exists(extract_dir)) {
    unlink(
      extract_dir,
      recursive=TRUE,
      force=TRUE
    )
  }

  dir.create(
    extract_dir,
    recursive=TRUE,
    showWarnings=FALSE
  )

  untar(
    tgz_path,
    exdir=extract_dir
  )

  all_files <- list.files(
    extract_dir,
    recursive=TRUE,
    full.names=TRUE
  )

  mmc5 <- all_files[
    basename(all_files) == "mmc5.xlsx"
  ]

  mmc6 <- all_files[
    basename(all_files) == "mmc6.xlsx"
  ]

  if (length(mmc5)==0 || length(mmc6)==0) {
    stop(
      paste0(
        "OA package extracted but mmc5/mmc6 not found. XLSX files found: ",
        paste(
          basename(
            all_files[
              str_detect(
                all_files,
                regex("\\.xlsx$", ignore_case=TRUE)
              )
            ]
          ),
          collapse=";"
        )
      )
    )
  }

  s1 <- file.path(
    dirCRA,
    "CRA002354_Table_S1_64_patients.xlsx"
  )

  s2 <- file.path(
    dirCRA,
    "CRA002354_Table_S2_131_samples.xlsx"
  )

  file.copy(
    mmc5[1],
    s1,
    overwrite=TRUE
  )

  file.copy(
    mmc6[1],
    s2,
    overwrite=TRUE
  )

  valid_s1 <- is_valid_xlsx(s1)
  valid_s2 <- is_valid_xlsx(s2)

  if (!valid_s1 || !valid_s2) {
    stop(
      paste0(
        "Extracted CRA supplements failed validation: S1=",
        valid_s1,
        ", S2=",
        valid_s2
      )
    )
  }

  prof1 <- profile_excel(
    s1,
    project
  )

  prof2 <- profile_excel(
    s2,
    project
  )

  prof <- bind_rows(
    prof1,
    prof2
  )

  write_excel_csv(
    prof,
    file.path(
      OUT_ROOT,
      "CRA002354_supplement_profile.csv"
    ),
    na=""
  )

  export_all_sheets(
    s1,
    "CRA002354_S1"
  )

  export_all_sheets(
    s2,
    "CRA002354_S2"
  )

  detail <- tibble(
    Resource=c(
      "Table_S1_64_patients",
      "Table_S2_131_samples"
    ),
    File=c(s1,s2),
    Bytes=c(
      file.info(s1)$size,
      file.info(s2)$size
    ),
    Valid_XLSX=c(
      valid_s1,
      valid_s2
    )
  )

  write_excel_csv(
    detail,
    file.path(
      OUT_ROOT,
      "CRA002354_supplement_validation.csv"
    ),
    na=""
  )

  cat(
    paste0(
      "CRA002354 COMPLETE: OA=",
      oa_status,
      " TGZ=",
      tgz_status,
      " S1=",
      file.info(s1)$size,
      " S2=",
      file.info(s2)$size,
      " ",
      Sys.time(),
      "\n"
    ),
    file=LOG,
    append=TRUE
  )

  tibble(
    Project=project,
    Download_Status=paste0(
      "OA:",
      oa_status,
      ";TGZ:",
      tgz_status
    ),
    Valid_XLSX=valid_s1 & valid_s2,
    Bytes=
      file.info(s1)$size +
      file.info(s2)$size,
    Sheets_Profiled=nrow(prof),
    Status="S1_S2_SUPPLEMENTS_RECOVERED_AND_VALIDATED"
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
    Valid_XLSX=FALSE,
    Bytes=0,
    Sheets_Profiled=0,
    Status=paste0("ERROR:", conditionMessage(e))
  )
})

# ============================================================
# Final
# ============================================================
summary <- bind_rows(
  status851,
  statusCRA
)

write_excel_csv(
  summary,
  file.path(
    OUT_ROOT,
    "V2_step78A2_status.csv"
  ),
  na=""
)

cat(
  paste0(
    "STEP78A2 COMPLETE: ",
    Sys.time(),
    "\n"
  ),
  file=LOG,
  append=TRUE
)

cat("\n=============================================\n")
cat("SEPSIS V2 STEP 78A2 COMPLETE\n")
cat("=============================================\n")
print(summary, n=Inf, width=Inf)
cat("\nOutput: ", OUT_ROOT, "\n", sep="")
