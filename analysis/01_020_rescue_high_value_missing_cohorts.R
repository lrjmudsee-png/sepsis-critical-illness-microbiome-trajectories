# ============================================================
# Sepsis V2 - Step 76
# Rescue high-value missing longitudinal cohorts
# R 4.4.0 / Windows
#
# Targets:
#   PRJNA1125274
#   PRJEB68229
#   PRJNA851469
#   CRA002354
#
# Main goals:
# 1) fetch ENA SAMPLE XML and expand all sample attributes
# 2) recover explicit patient/time fields where publicly present
# 3) reconcile PRJNA1125274 local 308-ish paired runs vs 289 ENA runs
# 4) fetch PRJEB68229 official supplementary DOCX
# 5) fetch PRJNA851469 Nature supplementary XLSX links
# 6) fetch CRA002354 GSA metadata + Supplementary Tables S1/S2
#
# IMPORTANT:
# - no patient ID is inferred from opaque sample prefixes
# - no cohort is automatically added to master in this step
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr","readxl","dplyr","tidyr","stringr","purrr",
  "tibble","httr2","xml2","rvest","officer"
)

missing <- pkgs[
  !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
  message("Installing missing packages: ", paste(missing, collapse=", "))
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
  library(rvest)
  library(officer)
})

PROJECT_ROOT <- "E:/sepsis_project"
DATA_ROOT <- file.path(PROJECT_ROOT, "data")
OUT_ROOT <- file.path(
  PROJECT_ROOT,
  "results",
  "V2_16_high_value_metadata_rescue"
)

dir.create(
  OUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)

ENA_PROJECTS <- c(
  "PRJNA1125274",
  "PRJEB68229",
  "PRJNA851469"
)

UA <- "SepsisV2-Step76-R44/1.0"

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

norm_id <- function(x) {
  x <- clean_chr(x)
  ifelse(
    is.na(x),
    NA_character_,
    toupper(gsub("[^A-Za-z0-9._-]", "", x))
  )
}

find_recursive <- function(folder, pattern) {
  if (!dir.exists(folder)) return(character())
  x <- list.files(folder, recursive=TRUE, full.names=TRUE)
  x[
    str_detect(
      basename(x),
      regex(pattern, ignore_case=TRUE)
    )
  ]
}

safe_tsv <- function(path) {
  if (length(path)==0 || is.na(path) || !file.exists(path)) return(NULL)

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

safe_csv <- function(path) {
  if (length(path)==0 || is.na(path) || !file.exists(path)) return(NULL)

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
  ) return("EXISTS")

  tryCatch(
    {
      resp <- request(url) |>
        req_user_agent(UA) |>
        req_headers(Accept="*/*") |>
        req_timeout(300) |>
        req_perform(path=dest)

      paste0("HTTP_", resp_status(resp))
    },
    error=function(e) {
      paste0("FAILED:", conditionMessage(e))
    }
  )
}

get_step75_ena <- function(project) {

  p <- find_recursive(
    file.path(DATA_ROOT, project),
    paste0(
      "^",
      project,
      "_ENA_read_run_step75\\.tsv$"
    )
  )

  if (length(p)==0) return(NULL)

  safe_tsv(p[1])
}

# ------------------------------------------------------------
# ENA SAMPLE XML parser
# ------------------------------------------------------------
fetch_sample_xml <- function(accession, dest) {

  if (
    is.na(accession) ||
    accession == ""
  ) return("NO_ACCESSION")

  if (
    file.exists(dest) &&
    file.info(dest)$size > 100
  ) return("EXISTS")

  url <- paste0(
    "https://www.ebi.ac.uk/ena/browser/api/xml/",
    accession
  )

  safe_download(
    url,
    dest
  )
}

parse_sample_xml <- function(path, project, accession) {

  if (
    !file.exists(path) ||
    file.info(path)$size < 100
  ) return(tibble())

  doc <- tryCatch(
    read_xml(path),
    error=function(e) NULL
  )

  if (is.null(doc)) return(tibble())

  alias <- xml_attr(
    xml_find_first(doc, "//SAMPLE"),
    "alias"
  )

  title <- xml_text(
    xml_find_first(doc, "//TITLE")
  )

  description <- xml_text(
    xml_find_first(doc, "//DESCRIPTION")
  )

  attrs <- xml_find_all(
    doc,
    "//SAMPLE_ATTRIBUTE"
  )

  if (length(attrs)==0) {

    return(
      tibble(
        project=project,
        sample_accession=accession,
        sample_alias=clean_chr(alias),
        sample_title=clean_chr(title),
        sample_description=clean_chr(description),
        attribute_tag=NA_character_,
        attribute_value=NA_character_
      )
    )
  }

  map_dfr(attrs, function(node) {

    tibble(
      project=project,
      sample_accession=accession,
      sample_alias=clean_chr(alias),
      sample_title=clean_chr(title),
      sample_description=clean_chr(description),

      attribute_tag=clean_chr(
        xml_text(
          xml_find_first(node, "./TAG")
        )
      ),

      attribute_value=clean_chr(
        xml_text(
          xml_find_first(node, "./VALUE")
        )
      )
    )
  })
}

pick_explicit_field <- function(df, regex_pattern) {

  hit <- names(df)[
    str_detect(
      names(df),
      regex(regex_pattern, ignore_case=TRUE)
    )
  ]

  if (length(hit)==0) {
    return(
      rep(
        NA_character_,
        nrow(df)
      )
    )
  }

  out <- rep(
    NA_character_,
    nrow(df)
  )

  for (nm in hit) {

    v <- clean_chr(
      df[[nm]]
    )

    take <- is.na(out) & !is.na(v)

    out[take] <- v[take]
  }

  out
}

extract_explicit_title_patient <- function(x) {

  x <- clean_chr(x)

  out <- rep(
    NA_character_,
    length(x)
  )

  if (length(x)==0) return(out)

  # Only patterns containing explicit semantic labels.
  pats <- c(
    "(?i)patient[_ -]*([A-Za-z0-9.-]+)",
    "(?i)subject[_ -]*([A-Za-z0-9.-]+)",
    "(?i)participant[_ -]*([A-Za-z0-9.-]+)"
  )

  for (pat in pats) {

    m <- str_match(
      x,
      pat
    )[,2]

    take <- is.na(out) & !is.na(m)
    out[take] <- m[take]
  }

  clean_chr(out)
}

extract_explicit_title_time <- function(x) {

  x <- clean_chr(x)

  out <- rep(
    NA_character_,
    length(x)
  )

  if (length(x)==0) return(out)

  # Explicit clinical/time labels only.
  patterns <- c(
    "(?i)\\b(T[0-9]+)\\b",
    "(?i)\\b(S[12])\\b",
    "(?i)\\b(day[_ -]*[0-9]+)\\b",
    "(?i)\\b(D[0-9]+)\\b",
    "(?i)\\b(baseline)\\b",
    "(?i)\\b(discharge)\\b",
    "(?i)\\b(sepsis[_ -]*onset)\\b",
    "(?i)\\b(admission)\\b"
  )

  for (pat in patterns) {

    m <- str_match(
      x,
      pat
    )[,2]

    take <- is.na(out) & !is.na(m)
    out[take] <- m[take]
  }

  clean_chr(out)
}

# ------------------------------------------------------------
# A. ENA projects: fetch all sample XML
# ------------------------------------------------------------
ena_rescue_summary <- list()

for (project in ENA_PROJECTS) {

  message("\n==== ", project, " ====")

  ena <- get_step75_ena(project)

  if (is.null(ena)) {
    warning("Step75 ENA metadata not found for ", project)
    next
  }

  xml_dir <- file.path(
    DATA_ROOT,
    project,
    "00_metadata",
    "step76_sample_xml"
  )

  dir.create(
    xml_dir,
    recursive=TRUE,
    showWarnings=FALSE
  )

  sample_ids <- unique(
    clean_chr(
      c(
        ena$sample_accession,
        ena$secondary_sample_accession
      )
    )
  )

  sample_ids <- sample_ids[
    !is.na(sample_ids)
  ]

  # Prefer accessions that look like sample records.
  # Duplicates are harmless but removed.
  sample_ids <- unique(sample_ids)

  xml_rows <- list()
  xml_log <- list()

  for (i in seq_along(sample_ids)) {

    acc <- sample_ids[i]

    dest <- file.path(
      xml_dir,
      paste0(acc, ".xml")
    )

    status <- fetch_sample_xml(
      acc,
      dest
    )

    xml_log[[i]] <- tibble(
      project=project,
      sample_accession=acc,
      status=status,
      local_path=dest
    )

    parsed <- parse_sample_xml(
      dest,
      project,
      acc
    )

    if (nrow(parsed)>0) {
      xml_rows[[length(xml_rows)+1]] <- parsed
    }

    if (
      i %% 50 == 0 ||
      i == length(sample_ids)
    ) {
      message(
        "Sample XML ",
        i,
        "/",
        length(sample_ids)
      )
    }
  }

  xml_long <- if (length(xml_rows)>0) {
    bind_rows(xml_rows)
  } else {
    tibble()
  }

  xml_log_df <- bind_rows(xml_log)

  write_excel_csv(
    xml_log_df,
    file.path(
      OUT_ROOT,
      paste0(
        project,
        "_sample_XML_download_log.csv"
      )
    ),
    na=""
  )

  if (nrow(xml_long)==0) next

  write_excel_csv(
    xml_long,
    file.path(
      OUT_ROOT,
      paste0(
        project,
        "_sample_attributes_LONG.csv"
      )
    ),
    na=""
  )

  # Make attribute names safe.
  xml_long2 <- xml_long |>
    mutate(
      attribute_tag_safe = attribute_tag |>
        coalesce("attribute") |>
        str_to_lower() |>
        str_replace_all("[^a-z0-9]+", "_") |>
        str_replace_all("^_|_$", "")
    )

  xml_wide <- xml_long2 |>
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
        paste(
          sort(
            unique(
              clean_chr(x)[!is.na(clean_chr(x))]
            )
          ),
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

  # Explicit structured patient/time fields.
  patient_structured <- pick_explicit_field(
    xml_wide,
    paste(
      c(
        "^patient",
        "^subject",
        "^participant",
        "host.*subject",
        "host.*patient",
        "individual.*id",
        "donor.*id"
      ),
      collapse="|"
    )
  )

  time_structured <- pick_explicit_field(
    xml_wide,
    paste(
      c(
        "timepoint",
        "time_point",
        "^visit",
        "^day$",
        "sampling.*day",
        "collection.*day",
        "icu.*day",
        "study.*day"
      ),
      collapse="|"
    )
  )

  patient_title <- extract_explicit_title_patient(
    paste(
      xml_wide$sample_title,
      xml_wide$sample_description
    )
  )

  time_title <- extract_explicit_title_time(
    paste(
      xml_wide$sample_title,
      xml_wide$sample_description
    )
  )

  xml_wide <- xml_wide |>
    mutate(
      patient_id_explicit = coalesce(
        clean_chr(patient_structured),
        clean_chr(patient_title)
      ),

      timepoint_explicit = coalesce(
        clean_chr(time_structured),
        clean_chr(time_title)
      ),

      patient_mapping_source = case_when(
        !is.na(clean_chr(patient_structured)) ~
          "ENA_SAMPLE_ATTRIBUTE",

        !is.na(clean_chr(patient_title)) ~
          "EXPLICIT_PATIENT_LABEL_IN_TITLE_OR_DESCRIPTION",

        TRUE ~ NA_character_
      ),

      time_mapping_source = case_when(
        !is.na(clean_chr(time_structured)) ~
          "ENA_SAMPLE_ATTRIBUTE",

        !is.na(clean_chr(time_title)) ~
          "EXPLICIT_TIME_LABEL_IN_TITLE_OR_DESCRIPTION",

        TRUE ~ NA_character_
      )
    )

  write_excel_csv(
    xml_wide,
    file.path(
      OUT_ROOT,
      paste0(
        project,
        "_sample_attributes_WIDE.csv"
      )
    ),
    na=""
  )

  # Run-level exact join through sample accessions.
  sample_lookup <- xml_wide |>
    select(
      sample_accession,
      patient_id_explicit,
      timepoint_explicit,
      patient_mapping_source,
      time_mapping_source,
      everything()
    )

  runmap_primary <- ena |>
    left_join(
      sample_lookup |>
        select(
          sample_accession,
          patient_id_explicit,
          timepoint_explicit,
          patient_mapping_source,
          time_mapping_source
        ),
      by="sample_accession"
    )

  # Secondary sample accession fallback.
  miss <- is.na(runmap_primary$patient_id_explicit) &
          is.na(runmap_primary$timepoint_explicit)

  if (any(miss)) {

    secondary_lookup <- sample_lookup |>
      rename(
        secondary_sample_accession=sample_accession
      )

    sec <- ena[miss, , drop=FALSE] |>
      left_join(
        secondary_lookup |>
          select(
            secondary_sample_accession,
            patient_id_explicit,
            timepoint_explicit,
            patient_mapping_source,
            time_mapping_source
          ),
        by="secondary_sample_accession"
      )

    runmap_primary[
      miss,
      c(
        "patient_id_explicit",
        "timepoint_explicit",
        "patient_mapping_source",
        "time_mapping_source"
      )
    ] <- sec[
      ,
      c(
        "patient_id_explicit",
        "timepoint_explicit",
        "patient_mapping_source",
        "time_mapping_source"
      )
    ]
  }

  runmap <- runmap_primary |>
    transmute(
      project=project,

      run_id=clean_chr(run_accession),

      sample_accession=clean_chr(sample_accession),

      secondary_sample_accession=
        clean_chr(secondary_sample_accession),

      sample_alias=clean_chr(sample_alias),

      sample_title=clean_chr(sample_title),

      sample_description=
        clean_chr(sample_description),

      collection_date=
        clean_chr(collection_date),

      body_site=
        clean_chr(host_body_site),

      patient_id=
        clean_chr(patient_id_explicit),

      timepoint=
        clean_chr(timepoint_explicit),

      patient_mapping_source,

      time_mapping_source
    )

  write_excel_csv(
    runmap,
    file.path(
      OUT_ROOT,
      paste0(
        project,
        "_EXPLICIT_run_patient_time_candidate_map.csv"
      )
    ),
    na=""
  )

  ena_rescue_summary[[project]] <- tibble(
    project=project,
    ena_runs=nrow(ena),
    unique_samples=nrow(xml_wide),
    runs_with_patient=sum(!is.na(runmap$patient_id)),
    runs_with_time=sum(!is.na(runmap$timepoint)),
    runs_with_patient_and_time=sum(
      !is.na(runmap$patient_id) &
      !is.na(runmap$timepoint)
    )
  )
}

# ------------------------------------------------------------
# B. PRJNA1125274 run reconciliation
# ------------------------------------------------------------
project <- "PRJNA1125274"

ena112 <- get_step75_ena(project)

local_fastqs <- list.files(
  file.path(DATA_ROOT, project),
  pattern="\\.fastq\\.gz$",
  recursive=TRUE,
  full.names=TRUE,
  ignore.case=TRUE
)

local_runs <- unique(
  str_extract(
    basename(local_fastqs),
    "SRR[0-9]+"
  )
)

local_runs <- sort(
  local_runs[
    !is.na(local_runs)
  ]
)

ena_runs <- if (!is.null(ena112)) {
  sort(
    unique(
      clean_chr(
        ena112$run_accession
      )
    )
  )
} else character()

srr_list_files <- find_recursive(
  file.path(DATA_ROOT, project),
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

reconcile112 <- bind_rows(
  tibble(
    Run_ID=local_runs,
    In_Local_FASTQ=TRUE
  ),
  tibble(
    Run_ID=ena_runs,
    In_ENA_Current=TRUE
  ),
  tibble(
    Run_ID=list_runs,
    In_SRR_Acc_List=TRUE
  )
) |>
  group_by(Run_ID) |>
  summarise(
    In_Local_FASTQ=any(
      coalesce(
        In_Local_FASTQ,
        FALSE
      )
    ),

    In_ENA_Current=any(
      coalesce(
        In_ENA_Current,
        FALSE
      )
    ),

    In_SRR_Acc_List=any(
      coalesce(
        In_SRR_Acc_List,
        FALSE
      )
    ),

    .groups="drop"
  ) |>
  arrange(Run_ID)

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
    "ENA_not_local",
    "SRR_list_not_current_ENA"
  ),

  Value=c(
    length(local_fastqs),
    length(local_runs),
    length(ena_runs),
    length(list_runs),

    sum(
      reconcile112$In_Local_FASTQ &
      !reconcile112$In_ENA_Current
    ),

    sum(
      reconcile112$In_ENA_Current &
      !reconcile112$In_Local_FASTQ
    ),

    sum(
      reconcile112$In_SRR_Acc_List &
      !reconcile112$In_ENA_Current
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

# ------------------------------------------------------------
# C. PRJEB68229 supplementary DOCX
# ------------------------------------------------------------
project <- "PRJEB68229"

supp682_url <- paste0(
  "https://media.springernature.com/original/",
  "springer-static/esm/art%3A10.1186%2Fs13613-024-01407-x/",
  "MediaObjects/13613_2024_1407_MOESM1_ESM.docx"
)

supp682_path <- file.path(
  DATA_ROOT,
  project,
  "00_metadata",
  "step76_public_supplement",
  "PRJEB68229_Additional_file_1.docx"
)

supp682_status <- safe_download(
  supp682_url,
  supp682_path
)

docx682 <- tibble()

if (
  file.exists(supp682_path) &&
  file.info(supp682_path)$size > 0
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

# ------------------------------------------------------------
# D. PRJNA851469 Nature supplementary files
# ------------------------------------------------------------
project <- "PRJNA851469"

nature_url <- "https://www.nature.com/articles/s41591-023-02243-5"

nature_links <- tibble()

try({

  page <- read_html(
    request(nature_url) |>
      req_user_agent(UA) |>
      req_timeout(120) |>
      req_perform() |>
      resp_body_string()
  )

  anchors <- html_elements(
    page,
    "a"
  )

  hrefs <- html_attr(
    anchors,
    "href"
  )

  texts <- html_text2(
    anchors
  )

  nature_links <- tibble(
    text=texts,
    href=hrefs
  ) |>
    filter(
      !is.na(href)
    ) |>
    filter(
      str_detect(
        href,
        regex(
          "\\.(xlsx|xls|csv|pdf)(\\?|$)|supplement",
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
  "step76_public_supplement"
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

    nm <- basename(
      sub(
        "\\?.*$",
        "",
        url
      )
    )

    if (nm=="") {
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

    status <- safe_download(
      url,
      dest
    )

    download851[[length(download851)+1]] <- tibble(
      text=nature_links$text[i],
      url=url,
      local_path=dest,
      status=status
    )
  }
}

download851_df <- if (length(download851)>0) {
  bind_rows(download851)
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

# Profile downloaded spreadsheets.
profile851 <- list()

for (p in list.files(
  supp851_dir,
  pattern="\\.(xlsx|xls)$",
  full.names=TRUE,
  ignore.case=TRUE
)) {

  sheets <- tryCatch(
    excel_sheets(p),
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

    profile851[[length(profile851)+1]] <- tibble(
      file=basename(p),
      sheet=sh,
      rows=nrow(x),
      columns=ncol(x),
      column_names=paste(
        names(x),
        collapse=";"
      ),
      patient_field=any(
        str_detect(
          names(x),
          regex(
            "patient|subject|participant",
            ignore_case=TRUE
          )
        )
      ),
      sample_field=any(
        str_detect(
          names(x),
          regex(
            "sample|biosample|specimen",
            ignore_case=TRUE
          )
        )
      ),
      time_field=any(
        str_detect(
          names(x),
          regex(
            "day|time|visit",
            ignore_case=TRUE
          )
        )
      ),
      infection_field=any(
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

profile851_df <- if (length(profile851)>0) {
  bind_rows(profile851)
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

# ------------------------------------------------------------
# E. CRA002354 public metadata + supplements
# ------------------------------------------------------------
project <- "CRA002354"

cra_dir <- file.path(
  DATA_ROOT,
  project,
  "00_metadata",
  "step76_public_metadata"
)

dir.create(
  cra_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

cra_page_url <- "https://ngdc.cncb.ac.cn/gsa/browse/CRA002354"

cra_links <- tibble()

try({

  html <- request(cra_page_url) |>
    req_user_agent(UA) |>
    req_timeout(120) |>
    req_perform() |>
    resp_body_string() |>
    read_html()

  a <- html_elements(
    html,
    "a"
  )

  cra_links <- tibble(
    text=html_text2(a),
    href=html_attr(a, "href")
  ) |>
    filter(
      !is.na(href)
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

# Try to locate official CRA002354.xlsx from page.
xlsx_links <- cra_links |>
  filter(
    str_detect(
      href,
      regex(
        "\\.xlsx(\\?|$)",
        ignore_case=TRUE
      )
    ) |
    str_detect(
      text,
      regex(
        "CRA002354\\.xlsx",
        ignore_case=TRUE
      )
    )
  )

cra_metadata_status <- "NOT_FOUND"

if (nrow(xlsx_links)>0) {

  for (i in seq_len(nrow(xlsx_links))) {

    url <- xlsx_links$href[i]

    if (str_starts(url, "/")) {
      url <- paste0(
        "https://ngdc.cncb.ac.cn",
        url
      )
    }

    dest <- file.path(
      cra_dir,
      "CRA002354_GSA_metadata.xlsx"
    )

    st <- safe_download(
      url,
      dest
    )

    if (
      file.exists(dest) &&
      file.info(dest)$size > 1000
    ) {
      cra_metadata_status <- st
      break
    }
  }
}

# Official supplementary tables from PMC.
supp_candidates <- list(

  Table_S1=c(
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC8377022/bin/mmc5.xlsx",
    "https://pmc.ncbi.nlm.nih.gov/articles/instance/8377022/bin/mmc5.xlsx"
  ),

  Table_S2=c(
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC8377022/bin/mmc6.xlsx",
    "https://pmc.ncbi.nlm.nih.gov/articles/instance/8377022/bin/mmc6.xlsx"
  )
)

cra_supp_log <- list()

for (nm in names(supp_candidates)) {

  urls <- supp_candidates[[nm]]

  dest <- file.path(
    cra_dir,
    paste0(
      "CRA002354_",
      nm,
      ".xlsx"
    )
  )

  status <- "FAILED"

  for (url in urls) {

    status <- safe_download(
      url,
      dest,
      force=TRUE
    )

    if (
      file.exists(dest) &&
      file.info(dest)$size > 1000
    ) break
  }

  cra_supp_log[[nm]] <- tibble(
    resource=nm,
    status=status,
    local_path=dest
  )
}

cra_supp_log_df <- bind_rows(
  cra_supp_log
)

write_excel_csv(
  cra_supp_log_df,
  file.path(
    OUT_ROOT,
    "CRA002354_supplement_download_log.csv"
  ),
  na=""
)

# Profile every recovered CRA xlsx.
cra_profiles <- list()

cra_xlsx <- list.files(
  cra_dir,
  pattern="\\.xlsx$",
  full.names=TRUE,
  ignore.case=TRUE
)

for (p in cra_xlsx) {

  sheets <- tryCatch(
    excel_sheets(p),
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

    cra_profiles[[length(cra_profiles)+1]] <- tibble(
      file=basename(p),
      sheet=sh,
      rows=nrow(x),
      columns=ncol(x),
      column_names=paste(
        names(x),
        collapse=";"
      ),

      patient_field=any(
        str_detect(
          names(x),
          regex(
            "patient|subject|id",
            ignore_case=TRUE
          )
        )
      ),

      sample_field=any(
        str_detect(
          names(x),
          regex(
            "sample|specimen|run|biosample",
            ignore_case=TRUE
          )
        )
      ),

      time_field=any(
        str_detect(
          names(x),
          regex(
            "day|time|date",
            ignore_case=TRUE
          )
        )
      ),

      mortality_field=any(
        str_detect(
          names(x),
          regex(
            "mortality|death|survival",
            ignore_case=TRUE
          )
        )
      ),

      sepsis_field=any(
        str_detect(
          names(x),
          regex(
            "sepsis|septic|shock",
            ignore_case=TRUE
          )
        )
      ),

      infection_field=any(
        str_detect(
          names(x),
          regex(
            "infection",
            ignore_case=TRUE
          )
        )
      ),

      antibiotic_field=any(
        str_detect(
          names(x),
          regex(
            "antibiotic|carbapenem",
            ignore_case=TRUE
          )
        )
      )
    )
  }
}

cra_profiles_df <- if (length(cra_profiles)>0) {
  bind_rows(cra_profiles)
} else {
  tibble()
}

write_excel_csv(
  cra_profiles_df,
  file.path(
    OUT_ROOT,
    "CRA002354_public_metadata_profile.csv"
  ),
  na=""
)

# ------------------------------------------------------------
# F. Summary
# ------------------------------------------------------------
ena_summary <- if (length(ena_rescue_summary)>0) {
  bind_rows(ena_rescue_summary)
} else {
  tibble()
}

write_excel_csv(
  ena_summary,
  file.path(
    OUT_ROOT,
    "V2_step76_ENA_explicit_mapping_summary.csv"
  ),
  na=""
)

decision <- tibble(
  project=c(
    "PRJNA1125274",
    "PRJEB68229",
    "PRJNA851469",
    "CRA002354"
  ),

  scientific_priority=c(
    "VERY_HIGH",
    "HIGH",
    "HIGH",
    "HIGH"
  ),

  current_action=c(
    "RECONCILE_289_VS_LOCAL_RUNS_AND_PARSE_XML",
    "PARSE_SAMPLE_XML_AND_SUPPLEMENT",
    "PARSE_SAMPLE_XML_AND_PUBLIC_NATURE_TABLES",
    "PARSE_GSA_METADATA_AND_SUPPLEMENT_S1_S2"
  ),

  automatic_master_inclusion_now=FALSE,

  reason=c(
    "Need exact patient/time mapping and run-set reconciliation first.",
    "Need public sample-to-patient linkage before longitudinal use.",
    "Use public exact identifiers only; restricted metadata must not be inferred.",
    "Use official patient/sample supplementary tables; do not infer from opaque sample-name prefixes."
  )
)

write_excel_csv(
  decision,
  file.path(
    OUT_ROOT,
    "V2_step76_rescue_decision.csv"
  ),
  na=""
)

cat("\n============================================================\n")
cat("SEPSIS V2 - STEP 76 COMPLETE\n")
cat("R version: ", R.version.string, "\n", sep="")
cat("============================================================\n\n")

cat("ENA EXPLICIT MAPPING SUMMARY:\n")
print(
  ena_summary,
  n=Inf,
  width=Inf
)

cat("\nPRJNA1125274 RUN RECONCILIATION:\n")
print(
  reconcile112_summary,
  n=Inf,
  width=Inf
)

cat("\nCRA002354 supplementary recovery:\n")
print(
  cra_supp_log_df,
  n=Inf,
  width=Inf
)

cat("\nOutput folder:\n", OUT_ROOT, "\n", sep="")
cat("============================================================\n")
