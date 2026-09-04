# ============================================================
# Sepsis V2 - Step 69
# Final public-metadata rescue for unresolved cohorts
# R 4.4.0 / Windows
#
# Targets:
#   PRJEB67798
#   PRJEB82425
#
# Goals:
# 1) Query BORIS Portal DSpace REST API for DOI 10.48620/372
#    and download public sample/metadata files, especially
#    rectal_resection_samples_used.xlsx.
# 2) Re-download and fully unpack PRJEB82425 supplementary DOCX,
#    export text/table structure and embedded media.
# 3) Examine ENA sample codes without converting them into patient
#    IDs unless an explicit public mapping is found.
# 4) Produce a final manual-metadata-needed decision file.
#
# No raw sequencing files are modified.
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "httr2","jsonlite","readr","readxl","dplyr",
  "purrr","stringr","tibble","tidyr","officer"
)

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) > 0) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(httr2)
  library(jsonlite)
  library(readr)
  library(readxl)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(officer)
})

DATA_ROOT <- "E:/sepsis_project/data"
OUT_ROOT <- "E:/sepsis_project/results/V2_09_final_metadata_rescue"
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

UA <- "SepsisV2-R44-Step69/1.0"

BORIS_API <- "https://boris-portal.unibe.ch/server/api"

SPRINGER_824_URL <- paste0(
  "https://media.springernature.com/original/",
  "springer-static/esm/art%3A10.1186%2Fs12879-025-10825-6/",
  "MediaObjects/12879_2025_10825_MOESM1_ESM.docx"
)


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

safe_get_json <- function(url) {

  tryCatch(
    {
      resp <- request(url) |>
        req_user_agent(UA) |>
        req_headers(Accept = "application/json") |>
        req_timeout(90) |>
        req_perform()

      if (resp_status(resp) >= 300) return(NULL)

      resp_body_json(
        resp,
        simplifyVector = FALSE
      )
    },
    error = function(e) NULL
  )
}


safe_download <- function(url, dest, force = FALSE) {

  dir.create(
    dirname(dest),
    recursive = TRUE,
    showWarnings = FALSE
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
        req_headers(Accept = "*/*") |>
        req_timeout(300) |>
        req_perform(path = dest)

      paste0("HTTP_", resp_status(resp))
    },
    error = function(e) {
      paste0(
        "FAILED:",
        conditionMessage(e)
      )
    }
  )
}


first_href <- function(obj, rel) {

  if (
    is.null(obj) ||
    is.null(obj$`_links`) ||
    is.null(obj$`_links`[[rel]]) ||
    is.null(obj$`_links`[[rel]]$href)
  ) {
    return(NA_character_)
  }

  obj$`_links`[[rel]]$href
}


flatten_search_objects <- function(obj) {

  if (
    is.null(obj) ||
    is.null(obj$`_embedded`) ||
    is.null(obj$`_embedded`$searchResult)
  ) {
    return(list())
  }

  sr <- obj$`_embedded`$searchResult

  if (
    is.null(sr$`_embedded`) ||
    is.null(sr$`_embedded`$objects)
  ) {
    return(list())
  }

  sr$`_embedded`$objects
}


extract_dspace_object <- function(x) {

  if (!is.null(x$`_embedded`$indexableObject)) {
    return(x$`_embedded`$indexableObject)
  }

  if (!is.null(x$indexableObject)) {
    return(x$indexableObject)
  }

  NULL
}


# ------------------------------------------------------------
# A. PRJEB67798 - BORIS Portal API search
# ------------------------------------------------------------
project <- "PRJEB67798"

boris_dir <- file.path(
  DATA_ROOT,
  project,
  "00_metadata",
  "auto_fetched",
  "10_BORIS_public_data"
)

dir.create(
  boris_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


search_terms <- c(
  "10.48620/372",
  "Intestinal dysbiosis intraoperative predictor septic complications",
  "rectal_resection_samples_used"
)

boris_objects <- list()
boris_search_log <- list()


for (term in search_terms) {

  url <- paste0(
    BORIS_API,
    "/discover/search/objects?query=",
    URLencode(term, reserved = TRUE),
    "&size=100"
  )

  obj <- safe_get_json(url)

  objects <- flatten_search_objects(obj)

  boris_search_log[[length(boris_search_log) + 1]] <- tibble(
    Query = term,
    URL = url,
    Objects_Returned = length(objects)
  )

  if (length(objects) > 0) {
    boris_objects <- c(
      boris_objects,
      objects
    )
  }
}


boris_search_log_df <- bind_rows(
  boris_search_log
)

write_excel_csv(
  boris_search_log_df,
  file.path(
    OUT_ROOT,
    "PRJEB67798_BORIS_search_log.csv"
  ),
  na = ""
)


# Collect unique item objects / UUIDs
item_rows <- list()

for (x in boris_objects) {

  ob <- extract_dspace_object(x)

  if (is.null(ob)) next

  uuid <- ob$uuid %||% NA_character_
  name <- ob$name %||% NA_character_
  type <- ob$type %||% NA_character_

  item_rows[[length(item_rows) + 1]] <- tibble(
    UUID = as.character(uuid),
    Name = as.character(name),
    Type = as.character(type)
  )
}


boris_items <- if (length(item_rows) > 0) {
  bind_rows(item_rows) |>
    distinct()
} else {
  tibble(
    UUID = character(),
    Name = character(),
    Type = character()
  )
}


write_excel_csv(
  boris_items,
  file.path(
    OUT_ROOT,
    "PRJEB67798_BORIS_search_objects.csv"
  ),
  na = ""
)


# ------------------------------------------------------------
# Traverse item -> bundles -> bitstreams
# ------------------------------------------------------------
bitstream_rows <- list()


for (i in seq_len(nrow(boris_items))) {

  uuid <- boris_items$UUID[i]

  if (
    is.na(uuid) ||
    uuid == ""
  ) next


  item_url <- paste0(
    BORIS_API,
    "/core/items/",
    uuid
  )

  item_obj <- safe_get_json(
    item_url
  )

  if (is.null(item_obj)) next


  bundles_url <- first_href(
    item_obj,
    "bundles"
  )

  if (is.na(bundles_url)) {
    bundles_url <- paste0(
      item_url,
      "/bundles"
    )
  }


  bundles_obj <- safe_get_json(
    bundles_url
  )

  if (
    is.null(bundles_obj) ||
    is.null(bundles_obj$`_embedded`)
  ) next


  bundle_list <- bundles_obj$`_embedded`$bundles

  if (is.null(bundle_list)) next


  for (bundle in bundle_list) {

    bundle_name <- bundle$name %||% NA_character_
    bundle_uuid <- bundle$uuid %||% NA_character_

    bitstreams_url <- first_href(
      bundle,
      "bitstreams"
    )

    if (
      is.na(bitstreams_url) &&
      !is.na(bundle_uuid)
    ) {
      bitstreams_url <- paste0(
        BORIS_API,
        "/core/bundles/",
        bundle_uuid,
        "/bitstreams"
      )
    }

    if (is.na(bitstreams_url)) next


    bits_obj <- safe_get_json(
      bitstreams_url
    )

    if (
      is.null(bits_obj) ||
      is.null(bits_obj$`_embedded`) ||
      is.null(bits_obj$`_embedded`$bitstreams)
    ) next


    for (bit in bits_obj$`_embedded`$bitstreams) {

      bit_uuid <- bit$uuid %||% NA_character_
      bit_name <- bit$name %||% NA_character_

      content_url <- first_href(
        bit,
        "content"
      )

      if (
        is.na(content_url) &&
        !is.na(bit_uuid)
      ) {
        content_url <- paste0(
          BORIS_API,
          "/core/bitstreams/",
          bit_uuid,
          "/content"
        )
      }


      bitstream_rows[[length(bitstream_rows) + 1]] <- tibble(
        Item_UUID = uuid,
        Item_Name = boris_items$Name[i],
        Bundle = as.character(bundle_name),
        Bitstream_UUID = as.character(bit_uuid),
        FileName = as.character(bit_name),
        Content_URL = as.character(content_url)
      )
    }
  }
}


boris_bitstreams <- if (
  length(bitstream_rows) > 0
) {

  bind_rows(bitstream_rows) |>
    distinct()

} else {

  tibble(
    Item_UUID = character(),
    Item_Name = character(),
    Bundle = character(),
    Bitstream_UUID = character(),
    FileName = character(),
    Content_URL = character()
  )
}


write_excel_csv(
  boris_bitstreams,
  file.path(
    OUT_ROOT,
    "PRJEB67798_BORIS_bitstream_inventory.csv"
  ),
  na = ""
)


# ------------------------------------------------------------
# Download likely useful public files
# ------------------------------------------------------------
download_log <- list()

if (nrow(boris_bitstreams) > 0) {

  useful <- boris_bitstreams |>
    filter(
      str_detect(
        FileName,
        regex(
          paste(
            c(
              "rectal_resection",
              "sample",
              "metadata",
              "meta_data",
              "patient",
              "clinical",
              "subject",
              "\\.xlsx$",
              "\\.xls$",
              "\\.csv$",
              "\\.tsv$"
            ),
            collapse = "|"
          ),
          ignore_case = TRUE
        )
      )
    )


  for (i in seq_len(nrow(useful))) {

    if (
      is.na(useful$Content_URL[i]) ||
      useful$Content_URL[i] == ""
    ) next

    safe_name <- gsub(
      '[<>:"/\\\\|?*]+',
      "_",
      useful$FileName[i]
    )

    dest <- file.path(
      boris_dir,
      safe_name
    )

    status <- safe_download(
      useful$Content_URL[i],
      dest
    )

    download_log[[length(download_log) + 1]] <- tibble(
      FileName = useful$FileName[i],
      Status = status,
      LocalPath = dest,
      Source = useful$Content_URL[i]
    )
  }
}


boris_download_log <- if (
  length(download_log) > 0
) {

  bind_rows(download_log)

} else {

  tibble(
    FileName = character(),
    Status = character(),
    LocalPath = character(),
    Source = character()
  )
}


write_excel_csv(
  boris_download_log,
  file.path(
    OUT_ROOT,
    "PRJEB67798_BORIS_download_log.csv"
  ),
  na = ""
)


# ------------------------------------------------------------
# Profile downloaded BORIS tables
# ------------------------------------------------------------
profile_table <- function(path) {

  ext <- tolower(
    tools::file_ext(path)
  )

  out <- list()

  if (ext %in% c("xlsx","xls")) {

    sheets <- tryCatch(
      excel_sheets(path),
      error = function(e) character()
    )

    for (sh in sheets) {

      x <- tryCatch(
        read_excel(
          path,
          sheet = sh,
          n_max = 5000,
          .name_repair = "unique"
        ),
        error = function(e) NULL
      )

      if (is.null(x)) next

      out[[length(out) + 1]] <- tibble(
        File = basename(path),
        Sheet = sh,
        Rows = nrow(x),
        Columns = ncol(x),
        ColumnNames = paste(
          names(x),
          collapse = ";"
        ),
        Has_Patient_Field = any(
          str_detect(
            names(x),
            regex(
              "patient|subject|participant",
              ignore_case = TRUE
            )
          )
        ),
        Has_Sample_Field = any(
          str_detect(
            names(x),
            regex(
              "sample|biosample|specimen",
              ignore_case = TRUE
            )
          )
        ),
        Has_Time_Field = any(
          str_detect(
            names(x),
            regex(
              "time|day|visit|date|T1|T2|T3",
              ignore_case = TRUE
            )
          )
        ),
        Has_Outcome_Field = any(
          str_detect(
            names(x),
            regex(
              "outcome|mortality|SSI|complication",
              ignore_case = TRUE
            )
          )
        )
      )
    }

  } else if (ext %in% c("csv","tsv","txt")) {

    x <- tryCatch(
      {
        if (ext == "csv") {
          read_csv(
            path,
            n_max = 5000,
            show_col_types = FALSE
          )
        } else {
          read_tsv(
            path,
            n_max = 5000,
            show_col_types = FALSE
          )
        }
      },
      error = function(e) NULL
    )

    if (!is.null(x)) {

      out[[1]] <- tibble(
        File = basename(path),
        Sheet = "",
        Rows = nrow(x),
        Columns = ncol(x),
        ColumnNames = paste(
          names(x),
          collapse = ";"
        ),
        Has_Patient_Field = any(
          str_detect(
            names(x),
            regex(
              "patient|subject|participant",
              ignore_case = TRUE
            )
          )
        ),
        Has_Sample_Field = any(
          str_detect(
            names(x),
            regex(
              "sample|biosample|specimen",
              ignore_case = TRUE
            )
          )
        ),
        Has_Time_Field = any(
          str_detect(
            names(x),
            regex(
              "time|day|visit|date|T1|T2|T3",
              ignore_case = TRUE
            )
          )
        ),
        Has_Outcome_Field = any(
          str_detect(
            names(x),
            regex(
              "outcome|mortality|SSI|complication",
              ignore_case = TRUE
            )
          )
        )
      )
    }
  }

  if (length(out) > 0) {
    bind_rows(out)
  } else {
    tibble()
  }
}


boris_profiles <- list()

if (dir.exists(boris_dir)) {

  local_tables <- list.files(
    boris_dir,
    pattern = "\\.(xlsx|xls|csv|tsv|txt)$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  for (p in local_tables) {
    z <- profile_table(p)
    if (nrow(z) > 0) {
      boris_profiles[[length(boris_profiles) + 1]] <- z
    }
  }
}


boris_profile_df <- if (
  length(boris_profiles) > 0
) {

  bind_rows(boris_profiles)

} else {

  tibble()
}


write_excel_csv(
  boris_profile_df,
  file.path(
    OUT_ROOT,
    "PRJEB67798_BORIS_table_profile.csv"
  ),
  na = ""
)


public_677_patient_mapping <- (
  nrow(boris_profile_df) > 0 &&
  any(
    boris_profile_df$Has_Patient_Field &
    boris_profile_df$Has_Sample_Field
  )
)


# ============================================================
# B. PRJEB82425 - supplement and ENA diagnostic
# ============================================================
project <- "PRJEB82425"

supp_dir <- file.path(
  DATA_ROOT,
  project,
  "00_metadata",
  "auto_fetched",
  "07_known_public_supplement"
)

dir.create(
  supp_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


supp_docx <- file.path(
  supp_dir,
  "PRJEB82425_Supplementary_Material_1_STEP69.docx"
)


supp_status <- safe_download(
  SPRINGER_824_URL,
  supp_docx,
  force = TRUE
)


# ------------------------------------------------------------
# Full DOCX structure
# ------------------------------------------------------------
sum824 <- tibble()

if (
  file.exists(supp_docx) &&
  file.info(supp_docx)$size > 0
) {

  sum824 <- tryCatch(
    officer::docx_summary(
      officer::read_docx(
        supp_docx
      )
    ) |>
      as_tibble(),
    error = function(e) tibble()
  )
}


if (nrow(sum824) > 0) {

  write_excel_csv(
    sum824,
    file.path(
      OUT_ROOT,
      "PRJEB82425_DOCX_full_summary_STEP69.csv"
    ),
    na = ""
  )
}


table_cells824 <- tibble()

if (
  nrow(sum824) > 0 &&
  all(
    c(
      "table_index",
      "row_id",
      "cell_id",
      "text"
    ) %in% names(sum824)
  )
) {

  table_cells824 <- sum824 |>
    filter(
      !is.na(table_index),
      !is.na(row_id),
      !is.na(cell_id)
    ) |>
    select(
      table_index,
      row_id,
      cell_id,
      text
    ) |>
    arrange(
      table_index,
      row_id,
      cell_id
    )


  write_excel_csv(
    table_cells824,
    file.path(
      OUT_ROOT,
      "PRJEB82425_DOCX_table_cells_STEP69.csv"
    ),
    na = ""
  )
}


# ------------------------------------------------------------
# Extract embedded DOCX media for later visual review
# ------------------------------------------------------------
media_dir <- file.path(
  OUT_ROOT,
  "PRJEB82425_supplement_media"
)

dir.create(
  media_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


tmp_docx <- tempfile(
  "prjeb82425_docx_"
)

dir.create(
  tmp_docx
)


unzip_ok <- tryCatch(
  {
    unzip(
      supp_docx,
      exdir = tmp_docx
    )
    TRUE
  },
  error = function(e) FALSE
)


media_files <- character()

if (unzip_ok) {

  source_media <- file.path(
    tmp_docx,
    "word",
    "media"
  )

  if (dir.exists(source_media)) {

    media_files <- list.files(
      source_media,
      full.names = TRUE
    )

    if (length(media_files) > 0) {

      file.copy(
        media_files,
        media_dir,
        overwrite = TRUE
      )
    }
  }
}


unlink(
  tmp_docx,
  recursive = TRUE,
  force = TRUE
)


# ------------------------------------------------------------
# ENA metadata profile
# ------------------------------------------------------------
ena_files <- list.files(
  file.path(
    DATA_ROOT,
    project,
    "00_metadata"
  ),
  pattern = "ENA.*(tsv|csv)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)


ena824 <- NULL

if (length(ena_files) > 0) {

  p <- ena_files[1]

  ena824 <- tryCatch(
    {
      if (
        tolower(
          tools::file_ext(p)
        ) == "csv"
      ) {
        read_csv(
          p,
          show_col_types = FALSE
        )
      } else {
        read_tsv(
          p,
          show_col_types = FALSE
        )
      }
    },
    error = function(e) NULL
  )
}


ena_profile824 <- tibble()

if (!is.null(ena824)) {

  ena_profile824 <- map_dfr(
    names(ena824),
    function(nm) {

      x <- clean_chr(
        ena824[[nm]]
      )

      tibble(
        Column = nm,
        NonMissing = sum(
          !is.na(x)
        ),
        Unique_N = n_distinct(
          x,
          na.rm = TRUE
        ),
        Examples = paste(
          head(
            unique(
              x[
                !is.na(x)
              ]
            ),
            10
          ),
          collapse = " | "
        )
      )
    }
  )


  write_excel_csv(
    ena_profile824,
    file.path(
      OUT_ROOT,
      "PRJEB82425_ENA_column_profile.csv"
    ),
    na = ""
  )
}


# ------------------------------------------------------------
# Save sample-code structure ONLY as a diagnostic.
# It must not be treated as patient mapping.
# ------------------------------------------------------------
sample_code_diag <- tibble()

if (!is.null(ena824)) {

  sample_col <- names(ena824)[
    str_detect(
      names(ena824),
      regex(
        "^sample_alias$|sample.?name",
        ignore_case = TRUE
      )
    )
  ]

  run_col <- names(ena824)[
    str_detect(
      names(ena824),
      regex(
        "^run_accession$",
        ignore_case = TRUE
      )
    )
  ]

  date_col <- names(ena824)[
    str_detect(
      names(ena824),
      regex(
        "^collection_date$",
        ignore_case = TRUE
      )
    )
  ]


  if (
    length(sample_col) > 0 &&
    length(run_col) > 0
  ) {

    sid <- clean_chr(
      ena824[[sample_col[1]]]
    )

    sample_code_diag <- tibble(
      Sample_Alias = sid,
      Run_ID = clean_chr(
        ena824[[run_col[1]]]
      ),
      Collection_Date = if (
        length(date_col) > 0
      ) {
        clean_chr(
          ena824[[date_col[1]]]
        )
      } else {
        NA_character_
      },
      Numeric_Sample_Code = suppressWarnings(
        as.numeric(
          sid
        )
      ),
      Note = paste(
        "DIAGNOSTIC ONLY.",
        "Do not derive patient ID from this numeric code",
        "without an explicit public mapping."
      )
    )


    write_excel_csv(
      sample_code_diag,
      file.path(
        OUT_ROOT,
        "PRJEB82425_sample_code_DIAGNOSTIC_ONLY.csv"
      ),
      na = ""
    )
  }
}


# Does the supplement expose explicit sample/patient IDs?
supp_text <- if (
  nrow(sum824) > 0 &&
  "text" %in% names(sum824)
) {

  paste(
    clean_chr(
      sum824$text
    ),
    collapse = " "
  )

} else {

  ""
}


supp_has_patient_identifier <- str_detect(
  supp_text,
  regex(
    "patient.?id|subject.?id|participant.?id",
    ignore_case = TRUE
  )
)


supp_has_sample_identifier <- str_detect(
  supp_text,
  regex(
    "sample.?id|sample.?number|specimen.?id",
    ignore_case = TRUE
  )
)


public_824_patient_mapping <- (
  supp_has_patient_identifier &&
  supp_has_sample_identifier
)


# ============================================================
# C. Final manual-collection decision
# ============================================================
manual_decision <- tibble(

  Project = c(
    "PRJNA516701",
    "PRJNA578267",
    "PRJNA595346",
    "PRJNA884103",
    "PRJEB37289",
    "PRJEB67798",
    "PRJEB82425"
  ),

  Current_Metadata_Status = c(
    "READY",
    "READY",
    "READY_EXPLICIT_NCBI_BIOSAMPLE_MAPPING",
    "READY_ALL_RUNS_MAPPED",
    "STATIC_BACKGROUND_ONLY",
    if (
      public_677_patient_mapping
    ) {
      "PUBLIC_PATIENT_SAMPLE_MAPPING_RECOVERED"
    } else {
      "PUBLIC_PATIENT_MAPPING_NOT_RECOVERED"
    },
    if (
      public_824_patient_mapping
    ) {
      "SUPPLEMENT_EXPLICIT_MAPPING_RECOVERED"
    } else {
      "SUPPLEMENT_NO_EXPLICIT_PATIENT_SAMPLE_MAPPING"
    }
  ),

  Human_Manual_Work_Now = c(
    "NONE",
    "NONE",
    "NONE",
    "NONE",
    "NONE",
    if (
      public_677_patient_mapping
    ) {
      "NONE"
    } else {
      paste(
        "AUTHOR_CONTACT LIKELY only if this cohort is required",
        "for patient-level longitudinal analysis."
      )
    },
    if (
      public_824_patient_mapping
    ) {
      "NONE"
    } else {
      paste(
        "REVIEW EXTRACTED SUPPLEMENT MEDIA FIRST.",
        "If no explicit patient-sample map is present,",
        "author contact may be needed."
      )
    }
  ),

  Mandatory_For_Main_V2 = c(
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO"
  ),

  Recommendation = c(
    "Keep in ready metadata pool",
    "Keep as longitudinal non-sepsis control",
    "Keep gut subset only (rectal/stool)",
    "Keep as neonatal BSI external validation; shotgun processing separate",
    "Use as static ICU/sepsis background validation",
    if (
      public_677_patient_mapping
    ) {
      "Use recovered public mapping"
    } else {
      paste(
        "Do not delay V2.",
        "Exclude from patient-level analysis unless author data are obtained."
      )
    },
    if (
      public_824_patient_mapping
    ) {
      "Use recovered mapping"
    } else {
      paste(
        "Do not infer patient IDs from numeric sample aliases.",
        "Keep pending while reviewing supplement images."
      )
    }
  )
)


write_excel_csv(
  manual_decision,
  file.path(
    OUT_ROOT,
    "V2_FINAL_manual_metadata_needed_after69.csv"
  ),
  na = ""
)


# ============================================================
# Console
# ============================================================
cat(
  "\n============================================================\n"
)

cat(
  "SEPSIS V2 - STEP 69 COMPLETE\n"
)

cat(
  "R version: ",
  R.version.string,
  "\n",
  sep = ""
)

cat(
  "============================================================\n\n"
)

cat(
  "PRJEB67798 BORIS search objects: ",
  nrow(boris_items),
  "\n",
  sep = ""
)

cat(
  "PRJEB67798 BORIS bitstreams: ",
  nrow(boris_bitstreams),
  "\n",
  sep = ""
)

cat(
  "PRJEB67798 public patient mapping recovered: ",
  public_677_patient_mapping,
  "\n",
  sep = ""
)

cat(
  "PRJEB82425 supplement download: ",
  supp_status,
  "\n",
  sep = ""
)

cat(
  "PRJEB82425 extracted media files: ",
  length(media_files),
  "\n",
  sep = ""
)

cat(
  "PRJEB82425 explicit supplement mapping recovered: ",
  public_824_patient_mapping,
  "\n\n",
  sep = ""
)


print(
  manual_decision,
  n = Inf,
  width = Inf
)


cat(
  "\nOutputs:\n",
  OUT_ROOT,
  "\n"
)

cat(
  "============================================================\n"
)
