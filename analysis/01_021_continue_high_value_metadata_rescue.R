# ============================================================
# Sepsis V2 - Step 76B
# Continue/fix high-value missing cohort rescue
# R 4.4.0 / Windows
#
# Fixes:
# 1) PRJEB68229 pivot_wider name collision
# 2) distinguish primary vs secondary ENA sample accessions
# 3) build PRJEB68229 high-confidence patient/time candidate map
# 4) reconcile PRJNA1125274 local FASTQ vs current ENA
# 5) audit PRJNA1125274 A/B/C naming structure without forcing
#    patient/time interpretation
# 6) continue PRJNA851469 + CRA002354 public supplement rescue
#
# DOES NOT redownload already-complete sample XML.
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr","readxl","dplyr","tidyr","stringr",
  "purrr","tibble","httr2","xml2","rvest","officer"
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
  library(xml2)
  library(rvest)
  library(officer)
})

PROJECT_ROOT <- "E:/sepsis_project"
DATA_ROOT <- file.path(PROJECT_ROOT, "data")

STEP76_ROOT <- file.path(
  PROJECT_ROOT,
  "results",
  "V2_16_high_value_metadata_rescue"
)

OUT_ROOT <- file.path(
  PROJECT_ROOT,
  "results",
  "V2_16B_high_value_metadata_rescue_fixed"
)

dir.create(
  OUT_ROOT,
  recursive=TRUE,
  showWarnings=FALSE
)

UA <- "SepsisV2-Step76B-R44/1.0"


# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
clean_chr <- function(x) {

  x <- trimws(
    as.character(x)
  )

  x[
    is.na(x) |
    x == "" |
    tolower(x) %in% c(
      "na","nan","n/a","null","none"
    )
  ] <- NA_character_

  x
}


safe_csv <- function(path) {

  if (
    length(path)==0 ||
    is.na(path) ||
    !file.exists(path)
  ) return(NULL)

  tryCatch(
    suppressMessages(
      read_csv(
        path,
        show_col_types=FALSE,
        progress=FALSE,
        name_repair="unique"
      )
    ),
    error=function(e) NULL
  )
}


safe_tsv <- function(path) {

  if (
    length(path)==0 ||
    is.na(path) ||
    !file.exists(path)
  ) return(NULL)

  tryCatch(
    suppressMessages(
      read_tsv(
        path,
        show_col_types=FALSE,
        progress=FALSE,
        name_repair="unique"
      )
    ),
    error=function(e) NULL
  )
}


find_recursive <- function(folder, pattern) {

  if (!dir.exists(folder)) {
    return(character())
  }

  x <- list.files(
    folder,
    recursive=TRUE,
    full.names=TRUE
  )

  x[
    str_detect(
      basename(x),
      regex(
        pattern,
        ignore_case=TRUE
      )
    )
  ]
}


get_step75_ena <- function(project) {

  p <- find_recursive(
    file.path(
      DATA_ROOT,
      project
    ),
    paste0(
      "^",
      project,
      "_ENA_read_run_step75\\.tsv$"
    )
  )

  if (length(p)==0) return(NULL)

  safe_tsv(
    p[1]
  )
}


safe_download <- function(url, dest, force=FALSE) {

  dir.create(
    dirname(dest),
    recursive=TRUE,
    showWarnings=FALSE
  )

  if (
    file.exists(dest) &&
    file.info(dest)$size > 0 &&
    !force
  ) {
    return("EXISTS")
  }

  tryCatch(
    {

      resp <- request(url) |>
        req_user_agent(UA) |>
        req_headers(
          Accept="*/*"
        ) |>
        req_timeout(300) |>
        req_perform(
          path=dest
        )

      paste0(
        "HTTP_",
        resp_status(resp)
      )

    },
    error=function(e) {

      paste0(
        "FAILED:",
        conditionMessage(e)
      )
    }
  )
}


# ============================================================
# A. PRJEB68229 - FIX WIDE TABLE
# ============================================================
project <- "PRJEB68229"

long682_path <- file.path(
  STEP76_ROOT,
  "PRJEB68229_sample_attributes_LONG.csv"
)

long682 <- safe_csv(
  long682_path
)

if (is.null(long682)) {

  stop(
    "Cannot find PRJEB68229_sample_attributes_LONG.csv"
  )
}


# Prefix all attribute columns with attr_
# so an ENA tag such as sample_(Description cannot collide with
# the existing sample_description identifier column.
long682_fixed <- long682 |>
  mutate(

    attribute_tag_safe =
      attribute_tag |>
      coalesce("unknown_attribute") |>
      str_to_lower() |>
      str_replace_all(
        "[^a-z0-9]+",
        "_"
      ) |>
      str_replace_all(
        "^_|_$",
        ""
      ) |>
      paste0(
        "attr_",
        .
      )
  )


wide682 <- long682_fixed |>
  select(
    project,
    sample_accession,
    sample_alias,
    sample_title,
    sample_description,
    attribute_tag_safe,
    attribute_value
  ) |>
  distinct() |>
  pivot_wider(
    names_from=attribute_tag_safe,
    values_from=attribute_value,

    values_fn=function(x) {

      x <- sort(
        unique(
          clean_chr(x)
        )
      )

      x <- x[
        !is.na(x)
      ]

      paste(
        x,
        collapse=";"
      )
    },

    values_fill=""
  ) |>
  mutate(
    across(
      everything(),
      ~na_if(.x, "")
    )
  )


write_excel_csv(
  wide682,
  file.path(
    OUT_ROOT,
    "PRJEB68229_sample_attributes_WIDE_FIXED.csv"
  ),
  na=""
)


# ------------------------------------------------------------
# Logical sample deduplication
# ------------------------------------------------------------
# ENA exposes SAMEA and ERS representations of the same logical
# submitted sample. sample_alias is the stable submission alias.
logical682 <- wide682 |>
  mutate(
    accession_type = case_when(

      str_detect(
        sample_accession,
        "^SAMEA"
      ) ~ "SAMEA",

      str_detect(
        sample_accession,
        "^ERS"
      ) ~ "ERS",

      TRUE ~ "OTHER"
    )
  ) |>
  arrange(
    sample_alias,
    factor(
      accession_type,
      levels=c(
        "SAMEA",
        "ERS",
        "OTHER"
      )
    )
  ) |>
  group_by(
    sample_alias
  ) |>
  summarise(

    primary_sample_accession =
      first(
        sample_accession
      ),

    all_sample_accessions =
      paste(
        sort(
          unique(
            sample_accession
          )
        ),
        collapse=";"
      ),

    sample_title =
      first(
        clean_chr(
          sample_title
        )
      ),

    sample_description =
      first(
        clean_chr(
          sample_description
        )
      ),

    attr_sample_description =
      first(
        clean_chr(
          attr_sample_description
        )
      ),

    .groups="drop"
  )


write_excel_csv(
  logical682,
  file.path(
    OUT_ROOT,
    "PRJEB68229_logical_samples_deduplicated.csv"
  ),
  na=""
)


# ------------------------------------------------------------
# Parse the deposited alias structure
# ------------------------------------------------------------
parsed682 <- logical682 |>
  mutate(

    patient_candidate =
      str_match(
        sample_alias,
        "^([0-9]+)_(D1|D3)$"
      )[,2],

    time_raw =
      str_match(
        sample_alias,
        "^([0-9]+)_(D1|D3)$"
      )[,3],

    is_patient_sample =
      !is.na(
        patient_candidate
      ),

    is_NEC_control =
      str_detect(
        sample_alias,
        "^NEC_"
      )
  )


# ------------------------------------------------------------
# Join logical sample alias to ENA Run
# ------------------------------------------------------------
ena682 <- get_step75_ena(
  project
)

if (is.null(ena682)) {
  stop(
    "Cannot find Step75 ENA metadata for PRJEB68229."
  )
}


run682 <- ena682 |>
  transmute(
    Run_ID =
      clean_chr(
        run_accession
      ),

    Sample_Alias =
      clean_chr(
        sample_alias
      ),

    BioSample =
      clean_chr(
        sample_accession
      ),

    Secondary_BioSample =
      clean_chr(
        secondary_sample_accession
      ),

    Sample_Title =
      clean_chr(
        sample_title
      ),

    Collection_Date =
      clean_chr(
        collection_date
      )
  ) |>
  left_join(
    parsed682 |>
      select(
        sample_alias,
        patient_candidate,
        time_raw,
        is_patient_sample,
        is_NEC_control,
        attr_sample_description
      ),
    by=c(
      "Sample_Alias"="sample_alias"
    )
  )


# ------------------------------------------------------------
# Validation structure
# ------------------------------------------------------------
patient682 <- run682 |>
  filter(
    is_patient_sample
  ) |>
  distinct(
    patient_candidate,
    time_raw
  )


n_pat682 <- n_distinct(
  patient682$patient_candidate
)

n_d1_682 <- sum(
  patient682$time_raw == "D1"
)

n_d3_682 <- sum(
  patient682$time_raw == "D3"
)


d1_pat <- patient682 |>
  filter(
    time_raw == "D1"
  ) |>
  pull(
    patient_candidate
  )


d3_pat <- patient682 |>
  filter(
    time_raw == "D3"
  ) |>
  pull(
    patient_candidate
  )


structure_valid682 <- (
  n_pat682 == 96 &&
  n_d1_682 == 96 &&
  n_d3_682 == 94 &&
  all(
    d3_pat %in% d1_pat
  )
)


map682 <- run682 |>
  filter(
    is_patient_sample
  ) |>
  transmute(

    Project =
      "PRJEB68229",

    Patient_ID =
      patient_candidate,

    Run_ID,

    BioSample,

    Secondary_BioSample,

    Sample_ID =
      Sample_Alias,

    Time_Raw =
      time_raw,

    Sample_Description =
      attr_sample_description,

    Collection_Date,

    Mapping_Confidence = ifelse(
      structure_valid682,
      "HIGH_BY_ENA_SUBMISSION_STRUCTURE_AND_COHORT_COUNT",
      "REVIEW_REQUIRED"
    ),

    Mapping_Method =
      paste(
        "ENA sample alias numeric prefix grouped as patient;",
        "D1/D3 retained exactly as deposited time labels;",
        "validated against 96-patient cohort structure"
      )
  )


write_excel_csv(
  map682,
  file.path(
    OUT_ROOT,
    "PRJEB68229_patient_time_map_CANDIDATE_HIGH_CONFIDENCE.csv"
  ),
  na=""
)


qc682 <- tibble(

  Metric=c(
    "ENA_runs",
    "Logical_sample_aliases",
    "Patient_samples",
    "NEC_controls",
    "Candidate_patients",
    "Patients_D1",
    "Patients_D3",
    "D3_patients_not_in_D1",
    "Structure_validation_passed"
  ),

  Value=c(
    nrow(ena682),
    nrow(logical682),
    sum(
      parsed682$is_patient_sample
    ),
    sum(
      parsed682$is_NEC_control
    ),
    n_pat682,
    n_d1_682,
    n_d3_682,
    sum(
      !d3_pat %in% d1_pat
    ),
    structure_valid682
  )
)


write_excel_csv(
  qc682,
  file.path(
    OUT_ROOT,
    "PRJEB68229_mapping_structure_QC.csv"
  ),
  na=""
)


# ============================================================
# B. PRJEB68229 supplementary DOCX
# ============================================================
supp682_url <- paste0(
  "https://media.springernature.com/original/",
  "springer-static/esm/art%3A10.1186%2Fs13613-024-01407-x/",
  "MediaObjects/13613_2024_1407_MOESM1_ESM.docx"
)

supp682_path <- file.path(
  DATA_ROOT,
  project,
  "00_metadata",
  "step76B_public_supplement",
  "PRJEB68229_Additional_file_1.docx"
)

supp682_status <- safe_download(
  supp682_url,
  supp682_path
)


docx682 <- tibble()

if (
  file.exists(supp682_path) &&
  file.info(supp682_path)$size > 1000
) {

  docx682 <- tryCatch(

    officer::docx_summary(
      officer::read_docx(
        supp682_path
      )
    ) |>
      as_tibble(),

    error=function(e) tibble()
  )
}


if (nrow(docx682)>0) {

  write_excel_csv(
    docx682,
    file.path(
      OUT_ROOT,
      "PRJEB68229_supplement_DOCX_summary.csv"
    ),
    na=""
  )
}


# ============================================================
# C. PRJNA1125274 - clean primary 289-sample view
# ============================================================
project <- "PRJNA1125274"

wide112_path <- file.path(
  STEP76_ROOT,
  "PRJNA1125274_sample_attributes_WIDE.csv"
)

map112_path <- file.path(
  STEP76_ROOT,
  "PRJNA1125274_EXPLICIT_run_patient_time_candidate_map.csv"
)

wide112 <- safe_csv(
  wide112_path
)

map112 <- safe_csv(
  map112_path
)


if (
  is.null(wide112) ||
  is.null(map112)
) {

  stop(
    "Cannot find Step76 PRJNA1125274 files."
  )
}


# SAMN is the primary BioSample accession represented in the Run table.
primary112 <- wide112 |>
  filter(
    str_detect(
      sample_accession,
      "^SAMN"
    )
  )


write_excel_csv(
  primary112,
  file.path(
    OUT_ROOT,
    "PRJNA1125274_primary_BioSamples_289.csv"
  ),
  na=""
)


# ------------------------------------------------------------
# Naming-pattern audit only
# ------------------------------------------------------------
pattern112 <- map112 |>
  mutate(

    alias_number =
      str_match(
        sample_alias,
        "^([0-9]+)([ABC])-([A-Z]+)$"
      )[,2],

    alias_letter =
      str_match(
        sample_alias,
        "^([0-9]+)([ABC])-([A-Z]+)$"
      )[,3],

    alias_site =
      str_match(
        sample_alias,
        "^([0-9]+)([ABC])-([A-Z]+)$"
      )[,4],

    patient_candidate =
      ifelse(
        !is.na(alias_number) &
        !is.na(alias_site),

        paste0(
          alias_site,
          "_",
          alias_number
        ),

        NA_character_
      ),

    # Hypothesis only; not accepted into master yet.
    time_hypothesis = case_when(
      alias_letter == "A" ~ "baseline_candidate",
      alias_letter == "B" ~ "sepsis_onset_candidate",
      alias_letter == "C" ~ "discharge_candidate",
      TRUE ~ NA_character_
    ),

    mapping_status =
      "NAMING_PATTERN_HYPOTHESIS_NOT_FINAL"
  )


write_excel_csv(
  pattern112,
  file.path(
    OUT_ROOT,
    "PRJNA1125274_alias_pattern_audit_NOT_FINAL.csv"
  ),
  na=""
)


qc112_pattern <- pattern112 |>
  summarise(

    Runs =
      n(),

    A_count =
      sum(
        alias_letter == "A",
        na.rm=TRUE
      ),

    B_count =
      sum(
        alias_letter == "B",
        na.rm=TRUE
      ),

    C_count =
      sum(
        alias_letter == "C",
        na.rm=TRUE
      ),

    Unique_number_site_candidates =
      n_distinct(
        patient_candidate,
        na.rm=TRUE
      ),

    Explicit_ENA_patient_fields =
      sum(
        !is.na(
          patient_id
        )
      ),

    Explicit_ENA_time_fields =
      sum(
        !is.na(
          timepoint
        )
      )
  )


write_excel_csv(
  qc112_pattern,
  file.path(
    OUT_ROOT,
    "PRJNA1125274_alias_pattern_QC.csv"
  ),
  na=""
)


# ------------------------------------------------------------
# Current ENA vs local FASTQ reconciliation
# ------------------------------------------------------------
ena112 <- get_step75_ena(
  project
)


local_fastqs <- list.files(
  file.path(
    DATA_ROOT,
    project
  ),
  pattern="\\.fastq\\.gz$",
  recursive=TRUE,
  full.names=TRUE,
  ignore.case=TRUE
)


local_runs <- sort(
  unique(
    str_extract(
      basename(
        local_fastqs
      ),
      "SRR[0-9]+"
    )
  )
)

local_runs <- local_runs[
  !is.na(local_runs)
]


ena_runs <- if (
  !is.null(
    ena112
  )
) {

  sort(
    unique(
      clean_chr(
        ena112$run_accession
      )
    )
  )

} else {

  character()
}


srr_list_files <- find_recursive(
  file.path(
    DATA_ROOT,
    project
  ),
  "^SRR_Acc_List\\.txt$"
)


list_runs <- character()

if (length(srr_list_files)>0) {

  list_runs <- sort(
    unique(
      clean_chr(
        readLines(
          srr_list_files[1],
          warn=FALSE
        )
      )
    )
  )
}


all_runs112 <- sort(
  unique(
    c(
      local_runs,
      ena_runs,
      list_runs
    )
  )
)


reconcile112 <- tibble(
  Run_ID=all_runs112
) |>
  mutate(

    In_Local_FASTQ =
      Run_ID %in% local_runs,

    In_Current_ENA =
      Run_ID %in% ena_runs,

    In_SRR_Acc_List =
      Run_ID %in% list_runs,

    Status = case_when(

      In_Local_FASTQ &
      In_Current_ENA ~
        "LOCAL_AND_CURRENT_ENA",

      In_Local_FASTQ &
      !In_Current_ENA ~
        "LOCAL_NOT_IN_CURRENT_ENA",

      !In_Local_FASTQ &
      In_Current_ENA ~
        "CURRENT_ENA_NOT_LOCAL",

      TRUE ~
        "OTHER"
    )
  )


write_excel_csv(
  reconcile112,
  file.path(
    OUT_ROOT,
    "PRJNA1125274_run_reconciliation.csv"
  ),
  na=""
)


reconcile112_summary <- tibble(

  Metric=c(
    "FASTQ_files",
    "Unique_local_runs",
    "Current_ENA_runs",
    "SRR_Acc_List_runs",
    "Local_not_in_current_ENA",
    "Current_ENA_not_local"
  ),

  Value=c(
    length(local_fastqs),
    length(local_runs),
    length(ena_runs),
    length(list_runs),

    sum(
      reconcile112$Status ==
      "LOCAL_NOT_IN_CURRENT_ENA"
    ),

    sum(
      reconcile112$Status ==
      "CURRENT_ENA_NOT_LOCAL"
    )
  )
)


write_excel_csv(
  reconcile112_summary,
  file.path(
    OUT_ROOT,
    "PRJNA1125274_run_reconciliation_summary.csv"
  ),
  na=""
)


# ============================================================
# D. PRJNA851469 - continue Nature supplement rescue
# ============================================================
project <- "PRJNA851469"

nature_url <- "https://www.nature.com/articles/s41591-023-02243-5"

nature_links <- tibble(
  text=character(),
  href=character()
)


try({

  resp <- request(
    nature_url
  ) |>
    req_user_agent(
      UA
    ) |>
    req_timeout(
      120
    ) |>
    req_perform()


  page <- read_html(
    resp_body_string(
      resp
    )
  )


  a <- html_elements(
    page,
    "a"
  )


  nature_links <- tibble(

    text =
      html_text2(a),

    href =
      html_attr(
        a,
        "href"
      )

  ) |>
    filter(
      !is.na(
        href
      )
    ) |>
    filter(
      str_detect(
        paste(
          text,
          href
        ),
        regex(
          "supp|\\.xlsx|\\.xls|\\.csv",
          ignore_case=TRUE
        )
      )
    ) |>
    distinct()

}, silent=TRUE)


write_excel_csv(
  nature_links,
  file.path(
    OUT_ROOT,
    "PRJNA851469_Nature_supplement_links.csv"
  ),
  na=""
)


supp851_dir <- file.path(
  DATA_ROOT,
  project,
  "00_metadata",
  "step76B_public_supplement"
)

dir.create(
  supp851_dir,
  recursive=TRUE,
  showWarnings=FALSE
)


download851 <- list()


if (nrow(nature_links)>0) {

  for (i in seq_len(nrow(nature_links))) {

    url <- nature_links$href[i]


    if (
      str_starts(
        url,
        "/"
      )
    ) {

      url <- paste0(
        "https://www.nature.com",
        url
      )
    }


    if (
      !str_detect(
        url,
        regex(
          "\\.(xlsx|xls|csv)(\\?|$)",
          ignore_case=TRUE
        )
      )
    ) next


    clean_url <- sub(
      "\\?.*$",
      "",
      url
    )


    nm <- basename(
      clean_url
    )


    if (
      is.na(nm) ||
      nm == ""
    ) {

      nm <- paste0(
        "supplement_",
        i,
        ".xlsx"
      )
    }


    dest <- file.path(
      supp851_dir,
      nm
    )


    st <- safe_download(
      url,
      dest
    )


    download851[[
      length(download851)+1
    ]] <- tibble(

      text =
        nature_links$text[i],

      url =
        url,

      local_path =
        dest,

      status =
        st
    )
  }
}


download851_df <- if (
  length(download851)>0
) {

  bind_rows(
    download851
  )

} else {

  tibble(
    text=character(),
    url=character(),
    local_path=character(),
    status=character()
  )
}


write_excel_csv(
  download851_df,
  file.path(
    OUT_ROOT,
    "PRJNA851469_Nature_supplement_download_log.csv"
  ),
  na=""
)


profile851 <- list()


for (p in list.files(
  supp851_dir,
  pattern="\\.(xlsx|xls)$",
  full.names=TRUE,
  ignore.case=TRUE
)) {

  sheets <- tryCatch(
    excel_sheets(
      p
    ),
    error=function(e) character()
  )


  for (sh in sheets) {

    x <- tryCatch(

      read_excel(
        p,
        sheet=sh,
        n_max=5000,
        .name_repair="unique"
      ),

      error=function(e) NULL
    )


    if (is.null(x)) next


    profile851[[
      length(profile851)+1
    ]] <- tibble(

      file =
        basename(p),

      sheet =
        sh,

      rows =
        nrow(x),

      columns =
        ncol(x),

      column_names =
        paste(
          names(x),
          collapse=";"
        ),

      patient_field =
        any(
          str_detect(
            names(x),
            regex(
              "patient|subject|participant",
              ignore_case=TRUE
            )
          )
        ),

      sample_field =
        any(
          str_detect(
            names(x),
            regex(
              "sample|biosample|specimen",
              ignore_case=TRUE
            )
          )
        ),

      time_field =
        any(
          str_detect(
            names(x),
            regex(
              "day|time|visit",
              ignore_case=TRUE
            )
          )
        ),

      infection_field =
        any(
          str_detect(
            names(x),
            regex(
              "infection|nosocomial|culture",
              ignore_case=TRUE
            )
          )
        )
    )
  }
}


profile851_df <- if (
  length(profile851)>0
) {

  bind_rows(
    profile851
  )

} else {

  tibble()
}


write_excel_csv(
  profile851_df,
  file.path(
    OUT_ROOT,
    "PRJNA851469_Nature_supplement_profile.csv"
  ),
  na=""
)


# ============================================================
# E. CRA002354 - continue public metadata/supplement rescue
# ============================================================
project <- "CRA002354"

cra_dir <- file.path(
  DATA_ROOT,
  project,
  "00_metadata",
  "step76B_public_metadata"
)

dir.create(
  cra_dir,
  recursive=TRUE,
  showWarnings=FALSE
)


cra_page_url <- "https://ngdc.cncb.ac.cn/gsa/browse/CRA002354"


cra_links <- tibble(
  text=character(),
  href=character()
)


try({

  resp <- request(
    cra_page_url
  ) |>
    req_user_agent(
      UA
    ) |>
    req_timeout(
      120
    ) |>
    req_perform()


  page <- read_html(
    resp_body_string(
      resp
    )
  )


  a <- html_elements(
    page,
    "a"
  )


  cra_links <- tibble(

    text =
      html_text2(a),

    href =
      html_attr(
        a,
        "href"
      )

  ) |>
    filter(
      !is.na(
        href
      )
    ) |>
    distinct()

}, silent=TRUE)


write_excel_csv(
  cra_links,
  file.path(
    OUT_ROOT,
    "CRA002354_GSA_page_links.csv"
  ),
  na=""
)


# ------------------------------------------------------------
# Also inspect already-local CRA xlsx files
# ------------------------------------------------------------
cra_local_xlsx <- list.files(
  file.path(
    DATA_ROOT,
    project
  ),
  pattern="\\.xlsx$",
  recursive=TRUE,
  full.names=TRUE,
  ignore.case=TRUE
)


cra_profiles <- list()


for (p in cra_local_xlsx) {

  sheets <- tryCatch(
    excel_sheets(
      p
    ),
    error=function(e) character()
  )


  for (sh in sheets) {

    x <- tryCatch(

      read_excel(
        p,
        sheet=sh,
        n_max=5000,
        .name_repair="unique"
      ),

      error=function(e) NULL
    )


    if (is.null(x)) next


    cra_profiles[[
      length(cra_profiles)+1
    ]] <- tibble(

      file =
        basename(p),

      sheet =
        sh,

      rows =
        nrow(x),

      columns =
        ncol(x),

      column_names =
        paste(
          names(x),
          collapse=";"
        ),

      patient_field =
        any(
          str_detect(
            names(x),
            regex(
              "patient|subject|host",
              ignore_case=TRUE
            )
          )
        ),

      sample_field =
        any(
          str_detect(
            names(x),
            regex(
              "sample|biosample|run|specimen",
              ignore_case=TRUE
            )
          )
        ),

      time_field =
        any(
          str_detect(
            names(x),
            regex(
              "day|time|date",
              ignore_case=TRUE
            )
          )
        ),

      sepsis_field =
        any(
          str_detect(
            names(x),
            regex(
              "sepsis|septic|shock",
              ignore_case=TRUE
            )
          )
        ),

      mortality_field =
        any(
          str_detect(
            names(x),
            regex(
              "mortality|death|survival",
              ignore_case=TRUE
            )
          )
        )
    )
  }
}


cra_profiles_df <- if (
  length(cra_profiles)>0
) {

  bind_rows(
    cra_profiles
  )

} else {

  tibble()
}


write_excel_csv(
  cra_profiles_df,
  file.path(
    OUT_ROOT,
    "CRA002354_all_local_xlsx_profile.csv"
  ),
  na=""
)


# ============================================================
# F. Final status summary
# ============================================================
summary <- tibble(

  Project=c(
    "PRJEB68229",
    "PRJNA1125274",
    "PRJNA851469",
    "CRA002354"
  ),

  Status=c(

    if (
      structure_valid682
    ) {
      "PATIENT_TIME_STRUCTURE_RECOVERED_HIGH_CONFIDENCE"
    } else {
      "REVIEW_REQUIRED"
    },

    "289_RUN_SET_READY_BUT_PATIENT_TIME_NAMING_NOT_FINAL",

    if (
      nrow(
        profile851_df
      ) > 0
    ) {
      "PUBLIC_SUPPLEMENT_RECOVERED"
    } else {
      "PUBLIC_SUPPLEMENT_NOT_YET_RECOVERED"
    },

    if (
      nrow(
        cra_profiles_df
      ) > 0
    ) {
      "LOCAL_PUBLIC_METADATA_PROFILED"
    } else {
      "CRA_METADATA_STILL_NEEDS_RECOVERY"
    }
  ),

  Next_Action=c(

    "Review D1/D3 interpretation, then eligible for ICU outcome longitudinal module.",

    paste(
      "Use only current 289 ENA Runs;",
      "seek source validation for A/B/C and patient-center naming before master inclusion."
    ),

    "Inspect public supplement fields for exact patient/sample/time linkage.",

    "Inspect CRA local XLSX and public supplement/GSA links for exact patient/time mapping."
  )
)


write_excel_csv(
  summary,
  file.path(
    OUT_ROOT,
    "V2_step76B_final_status.csv"
  ),
  na=""
)


cat("\n============================================================\n")
cat("SEPSIS V2 - STEP 76B COMPLETE\n")
cat("R version: ", R.version.string, "\n", sep="")
cat("============================================================\n\n")

cat("PRJEB68229 QC:\n")
print(
  qc682,
  n=Inf,
  width=Inf
)

cat("\nPRJNA1125274 naming audit:\n")
print(
  qc112_pattern,
  n=Inf,
  width=Inf
)

cat("\nPRJNA1125274 Run reconciliation:\n")
print(
  reconcile112_summary,
  n=Inf,
  width=Inf
)

cat("\nFinal status:\n")
print(
  summary,
  n=Inf,
  width=Inf
)

cat("\nOutput folder:\n", OUT_ROOT, "\n", sep="")
cat("============================================================\n")
