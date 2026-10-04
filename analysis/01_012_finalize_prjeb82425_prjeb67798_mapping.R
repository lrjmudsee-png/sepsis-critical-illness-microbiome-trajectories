# ============================================================
# Sepsis V2 - Step 70
# Final exact mapping rescue for PRJEB82425 + PRJEB67798
# R 4.4.0 / Windows
#
# PRINCIPLES
# - Use only explicit public identifiers.
# - Never infer patient IDs from numeric sample_alias.
# - PRJEB82425 patient/time come explicitly from sample_title.
# - Body site is classified only when ENA experiment XML
#   explicitly identifies V1-V2 vs V3-V4 or related protocol text.
# - PRJEB67798 uses exact SampleID -> ENA identifier matches.
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr","readxl","dplyr","tidyr","stringr",
  "purrr","tibble","httr2","xml2"
)

missing <- pkgs[
  !vapply(
    pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing) > 0) {

  message(
    "Installing missing packages: ",
    paste(missing, collapse = ", ")
  )

  install.packages(
    missing,
    repos = "https://cloud.r-project.org"
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
})


# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
DATA_ROOT <- "E:/sepsis_project/data"

STEP68_ROOT <- paste0(
  "E:/sepsis_project/results/",
  "V2_08_remaining_metadata_resolution"
)

STEP69_ROOT <- paste0(
  "E:/sepsis_project/results/",
  "V2_09_final_metadata_rescue"
)

OUT_ROOT <- paste0(
  "E:/sepsis_project/results/",
  "V2_10_exact_metadata_mapping"
)

dir.create(
  OUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)

UA <- "SepsisV2-R44-Step70/1.0"


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
      "na",
      "nan",
      "n/a",
      "null",
      "none"
    )
  ] <- NA_character_

  x
}


norm_id <- function(x) {

  x <- clean_chr(x)

  ifelse(
    is.na(x),
    NA_character_,
    toupper(
      gsub(
        "[^A-Za-z0-9._-]",
        "",
        x
      )
    )
  )
}


find_latest <- function(
  folder,
  pattern
) {

  if (!dir.exists(folder)) {
    return(NA_character_)
  }

  x <- list.files(
    folder,
    pattern = pattern,
    full.names = TRUE
  )

  if (length(x) == 0) {
    return(NA_character_)
  }

  x[
    which.max(
      file.info(x)$mtime
    )
  ]
}


find_recursive <- function(
  folder,
  pattern
) {

  if (!dir.exists(folder)) {
    return(character())
  }

  x <- list.files(
    folder,
    recursive = TRUE,
    full.names = TRUE
  )

  x[
    str_detect(
      basename(x),
      regex(
        pattern,
        ignore_case = TRUE
      )
    )
  ]
}


safe_csv <- function(path) {

  if (
    length(path) == 0 ||
    is.na(path) ||
    !file.exists(path)
  ) {
    return(NULL)
  }

  tryCatch(
    suppressMessages(
      read_csv(
        path,
        show_col_types = FALSE,
        progress = FALSE,
        name_repair = "unique"
      )
    ),
    error = function(e) NULL
  )
}


safe_tsv <- function(path) {

  if (
    length(path) == 0 ||
    is.na(path) ||
    !file.exists(path)
  ) {
    return(NULL)
  }

  tryCatch(
    suppressMessages(
      read_tsv(
        path,
        show_col_types = FALSE,
        progress = FALSE,
        name_repair = "unique"
      )
    ),
    error = function(e) NULL
  )
}


read_ena <- function(project) {

  base <- file.path(
    DATA_ROOT,
    project,
    "00_metadata"
  )

  p <- find_recursive(
    base,
    "^01_ENA_read_run_metadata\\.tsv$"
  )

  if (length(p) == 0) {

    p <- find_recursive(
      base,
      paste0(
        project,
        ".*manifest.*\\.tsv$"
      )
    )
  }

  if (length(p) == 0) {
    return(NULL)
  }

  safe_tsv(
    p[1]
  )
}


download_xml <- function(
  accession,
  dest
) {

  if (
    file.exists(dest) &&
    file.info(dest)$size > 100
  ) {
    return("EXISTS")
  }

  url <- paste0(
    "https://www.ebi.ac.uk/ena/browser/api/xml/",
    accession
  )

  tryCatch(
    {

      resp <- request(url) |>
        req_user_agent(UA) |>
        req_headers(
          Accept = "application/xml"
        ) |>
        req_timeout(90) |>
        req_perform(
          path = dest
        )

      paste0(
        "HTTP_",
        resp_status(resp)
      )
    },
    error = function(e) {

      paste0(
        "FAILED:",
        conditionMessage(e)
      )
    }
  )
}


xml_protocol_text <- function(path) {

  if (
    !file.exists(path) ||
    file.info(path)$size < 100
  ) {
    return(NA_character_)
  }

  tryCatch(
    {

      doc <- read_xml(path)

      nodes <- xml_find_all(
        doc,
        paste0(
          "//DESIGN_DESCRIPTION | ",
          "//LIBRARY_NAME | ",
          "//LIBRARY_CONSTRUCTION_PROTOCOL | ",
          "//EXPERIMENT_ATTRIBUTE/TAG | ",
          "//EXPERIMENT_ATTRIBUTE/VALUE"
        )
      )

      txt <- paste(
        xml_text(nodes),
        collapse = " | "
      )

      clean_chr(txt)

    },
    error = function(e) {
      NA_character_
    }
  )
}


classify_amplicon_site <- function(txt) {

  txt2 <- tolower(
    ifelse(
      is.na(txt),
      "",
      txt
    )
  )

  # Rectal swabs in the publication were amplified at V3-V4.
  rectal_hit <- str_detect(
    txt2,
    regex(
      paste(
        c(
          "v3.?v4",
          "v3-v4",
          "341f",
          "805r",
          "806r",
          "cctacggg",
          "ggactachv"
        ),
        collapse = "|"
      ),
      ignore_case = TRUE
    )
  )

  # ETA/tracheal aspirates in the publication were amplified at V1-V2.
  tracheal_hit <- str_detect(
    txt2,
    regex(
      paste(
        c(
          "v1.?v2",
          "v1-v2",
          "27f",
          "338r",
          "agagtttgat",
          "tgctgcctcc"
        ),
        collapse = "|"
      ),
      ignore_case = TRUE
    )
  )

  case_when(
    rectal_hit & !tracheal_hit ~ "rectal",
    tracheal_hit & !rectal_hit ~ "tracheal_ETA",
    TRUE ~ NA_character_
  )
}


count_longitudinal <- function(
  df,
  patient_col = "Patient_ID",
  time_col = "Time_Raw"
) {

  if (
    is.null(df) ||
    nrow(df) == 0
  ) {
    return(
      tibble(
        Patients = 0,
        Patients_GE2_Timepoints = 0,
        Patients_GE3_Timepoints = 0
      )
    )
  }

  tmp <- df |>
    filter(
      !is.na(.data[[patient_col]]),
      !is.na(.data[[time_col]])
    ) |>
    distinct(
      .data[[patient_col]],
      .data[[time_col]]
    ) |>
    count(
      .data[[patient_col]],
      name = "n_time"
    )

  tibble(
    Patients = nrow(tmp),
    Patients_GE2_Timepoints = sum(
      tmp$n_time >= 2
    ),
    Patients_GE3_Timepoints = sum(
      tmp$n_time >= 3
    )
  )
}


# ============================================================
# A. PRJEB82425
# ============================================================
project824 <- "PRJEB82425"

ena824 <- read_ena(
  project824
)

map824 <- tibble()
xml_log824 <- tibble()


if (!is.null(ena824)) {

  required_cols <- c(
    "run_accession",
    "sample_accession",
    "sample_alias",
    "sample_title",
    "sample_description",
    "experiment_accession"
  )

  missing_cols <- setdiff(
    required_cols,
    names(ena824)
  )

  if (length(missing_cols) > 0) {

    stop(
      paste0(
        "PRJEB82425 ENA metadata missing columns: ",
        paste(
          missing_cols,
          collapse = ", "
        )
      )
    )
  }


  # ----------------------------------------------------------
  # Parse explicit title:
  # P03_Infection_D5
  # P44_Inclusion
  # P15_Discharge
  # ----------------------------------------------------------
  parsed824 <- ena824 |>
    mutate(
      Sample_Title = clean_chr(
        sample_title
      ),

      Patient_ID = str_match(
        Sample_Title,
        "^(P[0-9]+)_"
      )[,2],

      Time_Raw = str_replace(
        Sample_Title,
        "^P[0-9]+_",
        ""
      ),

      Infection_Group = clean_chr(
        sample_description
      ),

      Time_Class = case_when(

        str_detect(
          Time_Raw,
          regex(
            "^Inclusion$",
            ignore_case = TRUE
          )
        ) ~ "Inclusion",

        str_detect(
          Time_Raw,
          regex(
            "^Infection_D1$",
            ignore_case = TRUE
          )
        ) ~ "Infection_D1",

        str_detect(
          Time_Raw,
          regex(
            "^Infection_D5$",
            ignore_case = TRUE
          )
        ) ~ "Infection_D5",

        str_detect(
          Time_Raw,
          regex(
            "^Discharge$",
            ignore_case = TRUE
          )
        ) ~ "Discharge",

        TRUE ~ Time_Raw
      )
    )


  # ----------------------------------------------------------
  # Fetch experiment XML to identify V1-V2 vs V3-V4.
  # ----------------------------------------------------------
  xml_dir824 <- file.path(
    DATA_ROOT,
    project824,
    "00_metadata",
    "auto_fetched",
    "11_ENA_experiment_XML"
  )

  dir.create(
    xml_dir824,
    recursive = TRUE,
    showWarnings = FALSE
  )


  experiments <- sort(
    unique(
      clean_chr(
        parsed824$experiment_accession
      )
    )
  )

  experiments <- experiments[
    !is.na(experiments)
  ]


  xml_rows <- vector(
    "list",
    length(experiments)
  )


  for (i in seq_along(experiments)) {

    acc <- experiments[i]

    dest <- file.path(
      xml_dir824,
      paste0(
        acc,
        ".xml"
      )
    )

    status <- download_xml(
      acc,
      dest
    )

    protocol <- xml_protocol_text(
      dest
    )

    site <- classify_amplicon_site(
      protocol
    )

    xml_rows[[i]] <- tibble(
      experiment_accession = acc,
      XML_Status = status,
      Body_Site = site,
      Protocol_Text = protocol
    )

    if (
      i %% 20 == 0 ||
      i == length(experiments)
    ) {

      message(
        "PRJEB82425 experiment XML: ",
        i,
        "/",
        length(experiments)
      )
    }
  }


  xml_log824 <- bind_rows(
    xml_rows
  )


  write_excel_csv(
    xml_log824,
    file.path(
      OUT_ROOT,
      "PRJEB82425_experiment_XML_classification.csv"
    ),
    na = ""
  )


  # ----------------------------------------------------------
  # Merge explicit patient/time + body site
  # ----------------------------------------------------------
  map824 <- parsed824 |>
    left_join(
      xml_log824 |>
        select(
          experiment_accession,
          Body_Site,
          XML_Status
        ),
      by = "experiment_accession"
    ) |>
    transmute(
      Project = project824,

      Patient_ID = clean_chr(
        Patient_ID
      ),

      Sample_ID = clean_chr(
        sample_alias
      ),

      BioSample = clean_chr(
        sample_accession
      ),

      Run_ID = clean_chr(
        run_accession
      ),

      Experiment_ID = clean_chr(
        experiment_accession
      ),

      Sample_Title,

      Time_Raw = clean_chr(
        Time_Raw
      ),

      Time_Class = clean_chr(
        Time_Class
      ),

      Infection_Group = clean_chr(
        Infection_Group
      ),

      Body_Site = clean_chr(
        Body_Site
      ),

      Collection_Date_Shifted = if (
        "collection_date" %in% names(parsed824)
      ) {
        clean_chr(
          parsed824$collection_date
        )
      } else {
        NA_character_
      },

      Patient_Mapping_Method =
        "Explicit ENA sample_title prefix Pxx",

      Time_Mapping_Method =
        "Explicit ENA sample_title event label",

      Body_Site_Mapping_Method = ifelse(
        !is.na(Body_Site),
        "ENA experiment XML amplicon region",
        NA_character_
      )
    )


  write_excel_csv(
    map824,
    file.path(
      OUT_ROOT,
      "PRJEB82425_explicit_sample_patient_time_map.csv"
    ),
    na = ""
  )


  # Gut subset only when body site is explicitly resolved.
  gut824 <- map824 |>
    filter(
      Body_Site == "rectal"
    )


  write_excel_csv(
    gut824,
    file.path(
      OUT_ROOT,
      "PRJEB82425_GUT_ONLY_explicit_map.csv"
    ),
    na = ""
  )
}


# ============================================================
# B. PRJEB67798
# ============================================================
project677 <- "PRJEB67798"

ena677 <- read_ena(
  project677
)

boris_dir677 <- file.path(
  DATA_ROOT,
  project677,
  "00_metadata",
  "auto_fetched",
  "10_BORIS_public_data"
)


map677 <- tibble()
candidate677 <- tibble()


if (
  !is.null(ena677) &&
  dir.exists(boris_dir677)
) {

  # ----------------------------------------------------------
  # Primary public table:
  # #SampleID + condition + patients
  # ----------------------------------------------------------
  chip_path <- file.path(
    boris_dir677,
    "figure1_Chip2_3.tsv"
  )

  chip <- safe_tsv(
    chip_path
  )


  crosswalks <- list()


  if (!is.null(chip)) {

    sample_col <- names(chip)[
      str_detect(
        names(chip),
        regex(
          "sample",
          ignore_case = TRUE
        )
      )
    ]

    patient_col <- names(chip)[
      str_detect(
        names(chip),
        regex(
          "^patients?$|patient.?id|subject",
          ignore_case = TRUE
        )
      )
    ]

    condition_col <- names(chip)[
      str_detect(
        names(chip),
        regex(
          "^condition$|time|tp",
          ignore_case = TRUE
        )
      )
    ]


    if (
      length(sample_col) > 0 &&
      length(patient_col) > 0
    ) {

      cw <- tibble(
        Public_Sample_ID = clean_chr(
          chip[[sample_col[1]]]
        ),

        Public_Patient_ID = clean_chr(
          chip[[patient_col[1]]]
        ),

        Public_Time_Condition = if (
          length(condition_col) > 0
        ) {
          clean_chr(
            chip[[condition_col[1]]]
          )
        } else {
          NA_character_
        },

        Crosswalk_Source =
          "BORIS figure1_Chip2_3.tsv"
      ) |>
        filter(
          !is.na(Public_Sample_ID),
          !is.na(Public_Patient_ID)
        )


      crosswalks[[
        length(crosswalks) + 1
      ]] <- cw
    }
  }


  # ----------------------------------------------------------
  # Also inspect alpha-diversity sheets because they contain
  # Sample + patients explicitly.
  # ----------------------------------------------------------
  rectal_files <- list.files(
    boris_dir677,
    pattern = "^figure1_rectal_resection_.*\\.xlsx$",
    full.names = TRUE,
    ignore.case = TRUE
  )


  for (p in rectal_files) {

    sheets <- tryCatch(
      excel_sheets(p),
      error = function(e) character()
    )

    alpha_sheets <- sheets[
      str_detect(
        sheets,
        regex(
          "alpha",
          ignore_case = TRUE
        )
      )
    ]


    for (sh in alpha_sheets) {

      x <- tryCatch(
        read_excel(
          p,
          sheet = sh,
          .name_repair = "unique"
        ),
        error = function(e) NULL
      )

      if (is.null(x)) next


      s_col <- names(x)[
        str_detect(
          names(x),
          regex(
            "^sample$|sample.?id",
            ignore_case = TRUE
          )
        )
      ]

      p_col <- names(x)[
        str_detect(
          names(x),
          regex(
            "^patients?$|patient.?id",
            ignore_case = TRUE
          )
        )
      ]


      if (
        length(s_col) > 0 &&
        length(p_col) > 0
      ) {

        cw <- tibble(
          Public_Sample_ID = clean_chr(
            x[[s_col[1]]]
          ),

          Public_Patient_ID = clean_chr(
            x[[p_col[1]]]
          ),

          Public_Time_Condition = NA_character_,

          Crosswalk_Source = paste0(
            "BORIS ",
            basename(p),
            "::",
            sh
          )
        ) |>
          filter(
            !is.na(Public_Sample_ID),
            !is.na(Public_Patient_ID)
          )


        crosswalks[[
          length(crosswalks) + 1
        ]] <- cw
      }
    }
  }


  if (length(crosswalks) > 0) {

    public_crosswalk <- bind_rows(
      crosswalks
    ) |>
      distinct()


    write_excel_csv(
      public_crosswalk,
      file.path(
        OUT_ROOT,
        "PRJEB67798_public_sample_patient_crosswalk.csv"
      ),
      na = ""
    )


    # --------------------------------------------------------
    # Build long ENA identifier lookup.
    # Exact normalized matches only.
    # --------------------------------------------------------
    possible_id_cols <- intersect(
      c(
        "sample_alias",
        "sample_title",
        "sample_accession",
        "secondary_sample_accession",
        "run_alias",
        "experiment_alias",
        "experiment_title"
      ),
      names(ena677)
    )


    ena_long <- ena677 |>
      mutate(
        .ENA_Row = row_number()
      ) |>
      select(
        .ENA_Row,
        any_of(
          c(
            "run_accession",
            "sample_accession",
            "secondary_sample_accession",
            "sample_alias",
            "sample_title",
            "sample_description",
            "collection_date"
          )
        ),
        all_of(
          setdiff(
            possible_id_cols,
            c(
              "sample_accession",
              "secondary_sample_accession",
              "sample_alias",
              "sample_title"
            )
          )
        )
      ) |>
      pivot_longer(
        cols = all_of(
          possible_id_cols
        ),
        names_to = "ENA_ID_Field",
        values_to = "ENA_ID_Value"
      ) |>
      mutate(
        Match_ID = norm_id(
          ENA_ID_Value
        )
      ) |>
      filter(
        !is.na(Match_ID)
      )


    candidate677 <- public_crosswalk |>
      mutate(
        Match_ID = norm_id(
          Public_Sample_ID
        )
      ) |>
      inner_join(
        ena_long,
        by = "Match_ID"
      ) |>
      distinct()


    write_excel_csv(
      candidate677,
      file.path(
        OUT_ROOT,
        "PRJEB67798_exact_ID_match_candidates.csv"
      ),
      na = ""
    )


    # --------------------------------------------------------
    # Keep exact mappings only.
    # If one Public Sample ID maps to conflicting runs/patients,
    # flag it and do not silently choose one.
    # --------------------------------------------------------
    ambiguity677 <- candidate677 |>
      group_by(
        Public_Sample_ID,
        Public_Patient_ID
      ) |>
      summarise(
        N_ENA_Rows = n_distinct(
          .ENA_Row
        ),
        N_Runs = n_distinct(
          run_accession[
            !is.na(
              run_accession
            )
          ]
        ),
        .groups = "drop"
      )


    write_excel_csv(
      ambiguity677,
      file.path(
        OUT_ROOT,
        "PRJEB67798_exact_match_ambiguity_QC.csv"
      ),
      na = ""
    )


    good_ids <- ambiguity677 |>
      filter(
        N_ENA_Rows == 1,
        N_Runs <= 1
      ) |>
      select(
        Public_Sample_ID,
        Public_Patient_ID
      )


    map677 <- candidate677 |>
      inner_join(
        good_ids,
        by = c(
          "Public_Sample_ID",
          "Public_Patient_ID"
        )
      ) |>
      transmute(
        Project = project677,

        Patient_ID = clean_chr(
          Public_Patient_ID
        ),

        Sample_ID = clean_chr(
          Public_Sample_ID
        ),

        Run_ID = clean_chr(
          run_accession
        ),

        BioSample = clean_chr(
          sample_accession
        ),

        Time_Raw = clean_chr(
          Public_Time_Condition
        ),

        Sample_Title = if (
          "sample_title" %in% names(candidate677)
        ) {
          clean_chr(
            sample_title
          )
        } else {
          NA_character_
        },

        Sample_Description = if (
          "sample_description" %in% names(candidate677)
        ) {
          clean_chr(
            sample_description
          )
        } else {
          NA_character_
        },

        ENA_ID_Field,

        Crosswalk_Source,

        Mapping_Method =
          "Exact normalized public SampleID -> ENA identifier match"
      ) |>
      distinct()


    write_excel_csv(
      map677,
      file.path(
        OUT_ROOT,
        "PRJEB67798_exact_sample_patient_run_map.csv"
      ),
      na = ""
    )
  }
}


# ============================================================
# C. QC
# ============================================================
qc_rows <- list()


if (nrow(map824) > 0) {

  long824 <- count_longitudinal(
    map824
  )

  qc_rows[[length(qc_rows) + 1]] <- tibble(
    Project = "PRJEB82425",
    Rows = nrow(map824),
    Unique_Patients = n_distinct(
      map824$Patient_ID,
      na.rm = TRUE
    ),
    Unique_Samples = n_distinct(
      map824$Sample_ID,
      na.rm = TRUE
    ),
    Unique_Runs = n_distinct(
      map824$Run_ID,
      na.rm = TRUE
    ),
    Missing_Patient = sum(
      is.na(
        map824$Patient_ID
      )
    ),
    Missing_Time = sum(
      is.na(
        map824$Time_Raw
      )
    ),
    Body_Site_Resolved = sum(
      !is.na(
        map824$Body_Site
      )
    ),
    Gut_Rows = sum(
      map824$Body_Site == "rectal",
      na.rm = TRUE
    ),
    Patients_GE2_Timepoints =
      long824$Patients_GE2_Timepoints,
    Patients_GE3_Timepoints =
      long824$Patients_GE3_Timepoints
  )
}


if (nrow(map677) > 0) {

  long677 <- count_longitudinal(
    map677
  )

  qc_rows[[length(qc_rows) + 1]] <- tibble(
    Project = "PRJEB67798",
    Rows = nrow(map677),
    Unique_Patients = n_distinct(
      map677$Patient_ID,
      na.rm = TRUE
    ),
    Unique_Samples = n_distinct(
      map677$Sample_ID,
      na.rm = TRUE
    ),
    Unique_Runs = n_distinct(
      map677$Run_ID,
      na.rm = TRUE
    ),
    Missing_Patient = sum(
      is.na(
        map677$Patient_ID
      )
    ),
    Missing_Time = sum(
      is.na(
        map677$Time_Raw
      )
    ),
    Body_Site_Resolved = NA_integer_,
    Gut_Rows = NA_integer_,
    Patients_GE2_Timepoints =
      long677$Patients_GE2_Timepoints,
    Patients_GE3_Timepoints =
      long677$Patients_GE3_Timepoints
  )
}


qc <- if (
  length(qc_rows) > 0
) {

  bind_rows(
    qc_rows
  )

} else {

  tibble()
}


write_excel_csv(
  qc,
  file.path(
    OUT_ROOT,
    "V2_step70_mapping_QC.csv"
  ),
  na = ""
)


# ============================================================
# D. Decision table
# ============================================================
decision <- tibble(
  Project = c(
    "PRJEB82425",
    "PRJEB67798"
  ),

  Patient_Mapping = c(
    if (
      nrow(map824) > 0 &&
      sum(is.na(map824$Patient_ID)) == 0
    ) {
      "RECOVERED_EXPLICITLY_FROM_ENA_SAMPLE_TITLE"
    } else {
      "INCOMPLETE"
    },

    if (
      nrow(map677) > 0
    ) {
      "RECOVERED_FOR_EXACTLY_MATCHED_PUBLIC_SAMPLES"
    } else {
      "NO_EXACT_ENA_LINK_RECOVERED"
    }
  ),

  Time_Mapping = c(
    if (
      nrow(map824) > 0 &&
      sum(is.na(map824$Time_Raw)) == 0
    ) {
      "RECOVERED_EXPLICITLY_FROM_ENA_SAMPLE_TITLE"
    } else {
      "INCOMPLETE"
    },

    if (
      nrow(map677) > 0 &&
      any(
        !is.na(
          map677$Time_Raw
        )
      )
    ) {
      "PARTIALLY_RECOVERED_FROM_PUBLIC_CONDITION"
    } else {
      "NOT_YET_RESOLVED"
    }
  ),

  Body_Site = c(
    if (
      nrow(map824) > 0 &&
      sum(!is.na(map824$Body_Site)) > 0
    ) {
      paste0(
        "EXPLICIT_EXPERIMENT_PROTOCOL_RESOLVED_",
        sum(
          !is.na(
            map824$Body_Site
          )
        ),
        "_OF_",
        nrow(map824)
      )
    } else {
      "NOT_RESOLVED_FROM_EXPERIMENT_XML"
    },

    "RECTAL_RESECTION_PUBLIC_DATA"
  ),

  Manual_Action_Now = c(
    if (
      nrow(map824) > 0 &&
      sum(!is.na(map824$Body_Site)) == nrow(map824)
    ) {
      "NONE"
    } else {
      "NONE YET - inspect unresolved experiment protocol rows only"
    },

    if (
      nrow(map677) > 0
    ) {
      "NONE FOR RECOVERED EXACTLY MATCHED SUBSET"
    } else {
      "AUTHOR CONTACT ONLY IF THIS COHORT IS ESSENTIAL"
    }
  )
)


write_excel_csv(
  decision,
  file.path(
    OUT_ROOT,
    "V2_metadata_decision_after70.csv"
  ),
  na = ""
)


# ============================================================
# E. Console
# ============================================================
cat(
  "\n============================================================\n"
)

cat(
  "SEPSIS V2 - STEP 70 COMPLETE\n"
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
  "QC:\n"
)

print(
  qc,
  n = Inf,
  width = Inf
)


cat(
  "\nDecision:\n"
)

print(
  decision,
  n = Inf,
  width = Inf
)


cat(
  "\nOutput folder:\n",
  OUT_ROOT,
  "\n",
  sep = ""
)

cat(
  "============================================================\n"
)
