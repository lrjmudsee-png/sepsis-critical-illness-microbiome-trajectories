# ============================================================
# Sepsis V2 - Step92A
# Ecological Instability Index (EII) feasibility audit
#
# PURPOSE:
# Before building a predictive score, audit whether frozen
# longitudinal cohorts contain sufficient clinical variables.
#
# This step:
# 1. Audits clinical metadata availability
# 2. Searches for outcome/severity variables
# 3. Creates an EII-ready registry
# 4. Does NOT fit prediction models yet
#
# Target clinical variables:
# mortality
# SOFA
# APACHE
# lactate
# vasopressor
# ventilation
# ICU length of stay
# antibiotic exposure
# inflammatory markers
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c(
  "readr",
  "dplyr",
  "stringr",
  "tibble",
  "purrr"
)

missing <- pkgs[
  !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if(length(missing)){
  install.packages(
    missing,
    repos="https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
})

ROOT <- "E:/sepsis_project"

STEP87B <- file.path(
  ROOT,
  "results",
  "V2_27B_ANALYSIS_OBJECTS_AND_REPLICATE_FREEZE"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_32A_STEP92A_EII_FEASIBILITY_AUDIT"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

registry_path <- file.path(
  STEP87B,
  "V2_STEP87B_analysis_object_registry.csv"
)

if(!file.exists(registry_path)){
  stop("Step87B registry missing.")
}

registry <- read_csv(
  registry_path,
  show_col_types=FALSE
)

keyword_groups <- list(

  mortality=c(
    "death",
    "mortality",
    "survival",
    "dead",
    "28day"
  ),

  severity=c(
    "sofa",
    "apache",
    "score",
    "severity",
    "organ"
  ),

  shock=c(
    "vasopressor",
    "norepinephrine",
    "noradrenaline",
    "pressor",
    "shock"
  ),

  respiratory=c(
    "ventilation",
    "mechanical",
    "intub",
    "oxygen"
  ),

  inflammation=c(
    "crp",
    "pcr",
    "pct",
    "procalcitonin",
    "il6",
    "interleukin"
  ),

  antibiotic=c(
    "antibiotic",
    "abx",
    "drug",
    "treatment"
  ),

  icu_course=c(
    "icu",
    "length",
    "stay",
    "los",
    "hospital"
  ),

  laboratory=c(
    "lactate",
    "wbc",
    "white",
    "platelet",
    "creatinine",
    "bilirubin"
  )
)

score_column <- function(x){

  low <- tolower(x)

  scores <- map_int(
    keyword_groups,
    function(keys){
      sum(
        str_detect(
          low,
          paste(keys, collapse="|")
        )
      )
    }
  )

  max(scores)
}

category_column <- function(x){

  low <- tolower(x)

  hit <- names(keyword_groups)[
    map_lgl(
      keyword_groups,
      function(keys){
        str_detect(
          low,
          paste(keys, collapse="|")
        )
      }
    )
  ]

  if(length(hit)==0){
    return("")
  }

  paste(hit, collapse=";")
}

all_rows <- list()

for(proj in registry$project){

  row <- registry %>%
    filter(project==proj)

  obj_path <- row$analysis_object_path[1]

  if(!file.exists(obj_path)){
    next
  }

  obj <- readRDS(obj_path)

  md <- as_tibble(obj$metadata)

  if(!("patient_id" %in% names(md))){
    next
  }

  for(col in names(md)){

    if(col=="patient_id"){
      next
    }

    all_rows[[length(all_rows)+1]] <-
      tibble(
        project=proj,
        column=col,
        n_samples=nrow(md),
        n_nonmissing=sum(
          !is.na(md[[col]])
        ),
        n_unique=dplyr::n_distinct(
          md[[col]],
          na.rm=TRUE
        ),
        category=category_column(col),
        keyword_score=score_column(col)
      )
  }
}

audit <- bind_rows(all_rows) %>%
  arrange(
    desc(keyword_score)
  )

write_csv(
  audit,
  file.path(
    OUT,
    "V2_STEP92A_clinical_variable_inventory.csv"
  )
)

candidate <- audit %>%
  filter(
    keyword_score>0
  )

write_csv(
  candidate,
  file.path(
    OUT,
    "V2_STEP92A_EII_candidate_clinical_variables.csv"
  )
)

summary <- candidate %>%
  count(
    project,
    category,
    name="n_candidates"
  )

write_csv(
  summary,
  file.path(
    OUT,
    "V2_STEP92A_project_clinical_availability_summary.csv"
  )
)

writeLines(
  c(
    paste0(
      "Completed: ",
      Sys.time()
    ),
    "STEP92A COMPLETE",
    "Purpose: clinical feasibility audit before EII construction."
  ),
  file.path(
    OUT,
    "_STEP92A_COMPLETE.ok"
  )
)

cat(
  "STEP92A COMPLETE\n"
)
cat(
  "Output:",
  OUT,
  "\n"
)
