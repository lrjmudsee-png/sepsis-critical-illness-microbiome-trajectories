# ============================================================
# Sepsis V2 - Step 85B
# Taxonomy QC + non-target feature audit + library-depth audit
#
# PURPOSE
#   Review the 8 taxonomy-complete cohorts BEFORE freezing any
#   feature/sample exclusion rule.
#
# THIS STEP DOES NOT MODIFY 02_asv_selected.
# THIS STEP DOES NOT EXCLUDE LOW-DEPTH SAMPLES.
#
# It will:
#   1) match each final ASV table to SILVA 138.2 taxonomy;
#   2) quantify Bacteria / Archaea / unassigned / other kingdoms;
#   3) flag chloroplast / mitochondria-like annotations;
#   4) quantify BOTH ASV-count and read-abundance contribution;
#   5) create a PREVIEW table after retaining prokaryotic targets
#      and removing obvious organelle annotations;
#   6) audit sample library sizes before/after that preview;
#   7) report candidate depth thresholds but NOT apply one.
#
# R 4.4.0 / Windows
# ============================================================

options(stringsAsFactors = FALSE)

cran <- c("readr","dplyr","stringr","tibble","purrr")
missing <- cran[
  !vapply(cran, requireNamespace, logical(1), quietly=TRUE)
]

if (length(missing)) {
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

CANON <- file.path(
  ROOT,
  "data",
  "_V2_ANALYSIS_READY"
)

PROJECT_ROOT <- file.path(
  CANON,
  "01_PROJECTS"
)

OUT <- file.path(
  ROOT,
  "results",
  "V2_25B2_TAXONOMY_QC_AND_DEPTH_AUDIT_FIXED"
)

dir.create(
  OUT,
  recursive=TRUE,
  showWarnings=FALSE
)

LOG <- file.path(
  OUT,
  "_STEP85B2_runtime_checkpoints.txt"
)

ERR <- file.path(
  OUT,
  "_STEP85B2_FATAL_ERROR.txt"
)

if (file.exists(ERR)) unlink(ERR)

ck <- function(x) {
  cat(
    paste0(x,": ",Sys.time(),"\n"),
    file=LOG,
    append=TRUE
  )
}

safe_csv <- function(p) {
  suppressMessages(
    read_csv(
      p,
      show_col_types=FALSE,
      progress=FALSE,
      name_repair="unique"
    )
  )
}

valid_rds_matrix <- function(p) {
  if (!file.exists(p)) return(FALSE)

  tryCatch({
    x <- readRDS(p)
    is.matrix(x) && nrow(x)>0 && ncol(x)>0
  }, error=function(e) FALSE)
}

md5_one <- function(p) {
  if (!file.exists(p)) return(NA_character_)
  unname(tools::md5sum(p))
}

registry <- tribble(
  ~project, ~expected_samples, ~preferred_pattern,

  "PRJNA430161",
  26L,
  "seqtab_final",

  "PRJNA691455",
  49L,
  "seqtab_final",

  "PRJNA978257",
  13L,
  "seqtab_final",

  "PRJEB82425",
  92L,
  "seqtab_final_step84B",

  "PRJNA516701",
  170L,
  "seqtab_final_step84C_QC_REVIEWED",

  "PRJNA578267",
  207L,
  "seqtab_final_step84B4B_minParent8",

  "PRJNA851469",
  119L,
  "seqtab_final_step84B",

  "PRJNA1166732",
  120L,
  "seqtab_final_step84B"
)

choose_selected_asv <- function(project, preferred_pattern) {

  d <- file.path(
    PROJECT_ROOT,
    project,
    "02_asv_selected"
  )

  files <- list.files(
    d,
    pattern="\\.rds$",
    full.names=TRUE,
    recursive=FALSE
  )

  files <- files[
    vapply(
      files,
      valid_rds_matrix,
      logical(1)
    )
  ]

  if (!length(files)) {
    stop(
      project,
      ": no valid selected ASV matrix RDS."
    )
  }

  p <- files[
    str_detect(
      basename(files),
      regex(
        preferred_pattern,
        ignore_case=TRUE
      )
    )
  ]

  candidates <- if (length(p)) p else files

  if (length(candidates)>1) {

    hashes <- vapply(
      candidates,
      md5_one,
      character(1)
    )

    if (
      length(
        unique(hashes)
      )>1
    ) {
      stop(
        project,
        ": multiple non-identical selected ASV candidates."
      )
    }
  }

  candidates[1]
}

taxonomy_path <- function(project) {

  p <- file.path(
    PROJECT_ROOT,
    project,
    "04_taxonomy_qc",
    paste0(
      project,
      "_taxonomy_SILVA1382_minBoot80.csv"
    )
  )

  if (!file.exists(p)) {
    stop(
      project,
      ": taxonomy CSV missing: ",
      p
    )
  }

  p
}

normalize_taxonomy <- function(tax) {

  rank_cols <- intersect(
    c(
      "Kingdom",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus"
    ),
    names(tax)
  )

  if (!length(rank_cols)) {
    stop("No taxonomy rank columns detected.")
  }

  # Do NOT use select(.) inside mutate().
  # In the user's R/dplyr environment the magrittr "." pronoun is
  # not available in that nested apply() context, which caused the
  # original Step85B failure.
  rank_df <- tax[
    ,
    rank_cols,
    drop=FALSE
  ]

  taxonomy_concat <- apply(
    rank_df,
    1,
    function(x) {
      paste(
        x[
          !is.na(x) &
          x!=""
        ],
        collapse=";"
      )
    }
  )

  tax$taxonomy_concat <- taxonomy_concat

  # Guard in case a rank column is unexpectedly absent.
  if (!"Kingdom" %in% names(tax)) {
    tax$Kingdom <- NA_character_
  }

  tax |>
    mutate(
      kingdom_norm=case_when(
        is.na(Kingdom) |
        Kingdom=="" ~
          "UNASSIGNED",

        str_to_lower(Kingdom)=="bacteria" ~
          "BACTERIA",

        str_to_lower(Kingdom)=="archaea" ~
          "ARCHAEA",

        TRUE ~
          "OTHER_KINGDOM"
      ),

      flag_chloroplast=str_detect(
        str_to_lower(
          coalesce(
            taxonomy_concat,
            ""
          )
        ),
        "chloroplast|chloroplastida"
      ),

      flag_mitochondria=str_detect(
        str_to_lower(
          coalesce(
            taxonomy_concat,
            ""
          )
        ),
        "mitochond|mitochondria"
      ),

      flag_organelle=
        flag_chloroplast |
        flag_mitochondria,

      target_prokaryote=
        kingdom_norm %in%
        c(
          "BACTERIA",
          "ARCHAEA"
        ),

      preview_keep=
        target_prokaryote &
        !flag_organelle
    )
}

depth_thresholds <- c(
  100,
  500,
  1000,
  2000,
  5000,
  10000
)

main <- function() {

  ck("STEP85B2 STARTED")

  project_summaries <- list()
  depth_all <- list()
  depth_threshold_all <- list()
  feature_class_all <- list()

  for (i in seq_len(nrow(registry))) {

    project <- registry$project[i]

    ck(
      paste0(
        project,
        " STARTED"
      )
    )

    seq_path <- choose_selected_asv(
      project,
      registry$preferred_pattern[i]
    )

    tax_path <- taxonomy_path(
      project
    )

    seqtab <- readRDS(
      seq_path
    )

    tax <- safe_csv(
      tax_path
    )

    if (
      nrow(seqtab)!=
      registry$expected_samples[i]
    ) {
      stop(
        project,
        ": sample-count guard failed."
      )
    }

    if (
      !"ASV_sequence" %in%
      names(tax)
    ) {
      stop(
        project,
        ": ASV_sequence missing from taxonomy CSV."
      )
    }

    if (
      anyDuplicated(
        tax$ASV_sequence
      )
    ) {
      stop(
        project,
        ": duplicated taxonomy ASV_sequence."
      )
    }

    seqs <- colnames(seqtab)

    if (
      !setequal(
        seqs,
        tax$ASV_sequence
      )
    ) {
      stop(
        project,
        ": taxonomy sequences do not exactly match selected ASV table."
      )
    }

    tax <- tax[
      match(
        seqs,
        tax$ASV_sequence
      ),
      ,
      drop=FALSE
    ]

    if (
      !identical(
        tax$ASV_sequence,
        seqs
      )
    ) {
      stop(
        project,
        ": taxonomy reorder guard failed."
      )
    }

    tax <- normalize_taxonomy(
      tax
    )

    abund <- colSums(
      seqtab
    )

    tax <- tax |>
      mutate(
        total_abundance_from_seqtab=
          as.numeric(abund),

        abundance_fraction=
          total_abundance_from_seqtab /
          sum(
            total_abundance_from_seqtab
          )
      )

    feature_class <- tax |>
      mutate(
        project=project,
        feature_status=case_when(
          preview_keep ~
            "KEEP_TARGET_PROKARYOTE",

          flag_organelle ~
            "REMOVE_ORGANELLE",

          kingdom_norm=="UNASSIGNED" ~
            "REVIEW_UNASSIGNED_KINGDOM",

          kingdom_norm=="OTHER_KINGDOM" ~
            "REMOVE_NON_TARGET_KINGDOM",

          TRUE ~
            "REVIEW_OTHER"
        )
      ) |>
      count(
        project,
        feature_status,
        name="ASVs"
      )

    feature_reads <- tax |>
      mutate(
        project=project,
        feature_status=case_when(
          preview_keep ~
            "KEEP_TARGET_PROKARYOTE",

          flag_organelle ~
            "REMOVE_ORGANELLE",

          kingdom_norm=="UNASSIGNED" ~
            "REVIEW_UNASSIGNED_KINGDOM",

          kingdom_norm=="OTHER_KINGDOM" ~
            "REMOVE_NON_TARGET_KINGDOM",

          TRUE ~
            "REVIEW_OTHER"
        )
      ) |>
      group_by(
        project,
        feature_status
      ) |>
      summarise(
        reads=sum(
          total_abundance_from_seqtab
        ),
        .groups="drop"
      )

    feature_class <- full_join(
      feature_class,
      feature_reads,
      by=c(
        "project",
        "feature_status"
      )
    ) |>
      mutate(
        ASV_pct=round(
          100*ASVs/ncol(seqtab),
          3
        ),

        read_pct=round(
          100*reads/sum(seqtab),
          3
        )
      )

    feature_class_all[[
      length(feature_class_all)+1
    ]] <- feature_class

    write_excel_csv(
      tax,
      file.path(
        OUT,
        paste0(
          project,
          "_taxonomy_feature_flags.csv"
        )
      ),
      na=""
    )

    keep <- tax$preview_keep

    if (!any(keep)) {
      stop(
        project,
        ": preview filter retained zero ASVs."
      )
    }

    seq_preview <- seqtab[
      ,
      keep,
      drop=FALSE
    ]

    before_reads <- rowSums(
      seqtab
    )

    after_reads <- rowSums(
      seq_preview
    )

    depth <- tibble(
      project=project,
      sample_id=rownames(seqtab),
      reads_before_taxonomy_preview=
        before_reads,
      reads_after_taxonomy_preview=
        after_reads,
      reads_removed=
        before_reads-after_reads,
      removed_pct=round(
        ifelse(
          before_reads>0,
          100*(before_reads-after_reads)/before_reads,
          NA_real_
        ),
        3
      )
    )

    depth_all[[
      length(depth_all)+1
    ]] <- depth

    write_excel_csv(
      depth,
      file.path(
        OUT,
        paste0(
          project,
          "_sample_depth_before_after_taxonomy_preview.csv"
        )
      ),
      na=""
    )

    threshold_df <- map_dfr(
      depth_thresholds,
      function(thr) {

        tibble(
          project=project,
          threshold=thr,

          samples_below_before=
            sum(
              before_reads<thr
            ),

          samples_below_after=
            sum(
              after_reads<thr
            ),

          total_samples=
            length(after_reads),

          pct_below_after=round(
            100*mean(
              after_reads<thr
            ),
            2
          )
        )
      }
    )

    depth_threshold_all[[
      length(depth_threshold_all)+1
    ]] <- threshold_df

    project_summaries[[
      length(project_summaries)+1
    ]] <- tibble(
      project=project,
      samples=nrow(seqtab),

      ASVs_before=ncol(seqtab),
      reads_before=sum(seqtab),

      Bacteria_ASVs=sum(
        tax$kingdom_norm=="BACTERIA"
      ),

      Archaea_ASVs=sum(
        tax$kingdom_norm=="ARCHAEA"
      ),

      unassigned_kingdom_ASVs=sum(
        tax$kingdom_norm=="UNASSIGNED"
      ),

      other_kingdom_ASVs=sum(
        tax$kingdom_norm=="OTHER_KINGDOM"
      ),

      organelle_ASVs=sum(
        tax$flag_organelle
      ),

      chloroplast_ASVs=sum(
        tax$flag_chloroplast
      ),

      mitochondria_ASVs=sum(
        tax$flag_mitochondria
      ),

      preview_keep_ASVs=sum(
        keep
      ),

      preview_keep_ASV_pct=round(
        100*mean(
          keep
        ),
        2
      ),

      preview_keep_reads=sum(
        seq_preview
      ),

      preview_keep_read_pct=round(
        100*sum(
          seq_preview
        )/
        sum(
          seqtab
        ),
        2
      ),

      samples_zero_after_preview=sum(
        after_reads==0
      ),

      minimum_reads_after_preview=min(
        after_reads
      ),

      median_reads_after_preview=median(
        after_reads
      ),

      selected_asv_path=seq_path,
      taxonomy_path=tax_path
    )

    # Save preview only to results, NEVER promote.
    saveRDS(
      seq_preview,
      file.path(
        OUT,
        paste0(
          project,
          "_PREVIEW_prokaryote_nonorganelle_seqtab.rds"
        )
      )
    )

    ck(
      paste0(
        project,
        " COMPLETE"
      )
    )

    rm(
      seqtab,
      seq_preview,
      tax
    )

    gc(
      verbose=FALSE
    )
  }

  project_summary <- bind_rows(
    project_summaries
  )

  depth <- bind_rows(
    depth_all
  )

  thresholds <- bind_rows(
    depth_threshold_all
  )

  feature_class <- bind_rows(
    feature_class_all
  )

  write_excel_csv(
    project_summary,
    file.path(
      OUT,
      "V2_STEP85B_project_taxonomy_QC_summary.csv"
    ),
    na=""
  )

  write_excel_csv(
    feature_class,
    file.path(
      OUT,
      "V2_STEP85B_feature_class_read_abundance_summary.csv"
    ),
    na=""
  )

  write_excel_csv(
    depth,
    file.path(
      OUT,
      "V2_STEP85B_all_sample_depth_audit.csv"
    ),
    na=""
  )

  write_excel_csv(
    thresholds,
    file.path(
      OUT,
      "V2_STEP85B_depth_threshold_candidates.csv"
    ),
    na=""
  )

  global <- tibble(
    metric=c(
      "projects",
      "samples",
      "ASVs_before",
      "reads_before",
      "preview_keep_ASVs",
      "preview_keep_reads",
      "samples_zero_after_preview"
    ),
    value=c(
      nrow(project_summary),
      sum(project_summary$samples),
      sum(project_summary$ASVs_before),
      sum(project_summary$reads_before),
      sum(project_summary$preview_keep_ASVs),
      sum(project_summary$preview_keep_reads),
      sum(project_summary$samples_zero_after_preview)
    )
  )

  write_excel_csv(
    global,
    file.path(
      OUT,
      "V2_STEP85B_global_summary.csv"
    ),
    na=""
  )

  readme <- c(
    "SEPSIS V2 - STEP85B TAXONOMY QC + DEPTH AUDIT",
    paste0(
      "Created: ",
      Sys.time()
    ),
    "",
    "THIS IS AN AUDIT STEP.",
    "No selected ASV table is modified.",
    "No low-depth sample threshold is applied.",
    "",
    "PREVIEW FEATURE RULE",
    "Keep SILVA Kingdom Bacteria or Archaea;",
    "flag/remove obvious chloroplast or mitochondria-like annotations;",
    "do not automatically rescue Kingdom-unassigned features.",
    "",
    "WHY",
    "Before longitudinal alpha/beta diversity analysis, one common and documented target-feature rule and one sample-depth rule should be frozen across cohorts.",
    "",
    "NEXT",
    "Review project_taxonomy_QC_summary, feature_class_read_abundance_summary, and depth_threshold_candidates.",
    "Then freeze Step85C/86 filtering rules."
  )

  writeLines(
    readme,
    file.path(
      OUT,
      "README_STEP85B.txt"
    ),
    useBytes=TRUE
  )

  writeLines(
    c(
      paste0(
        "Completed: ",
        Sys.time()
      ),
      "Status: STEP85B2 COMPLETE",
      "No final ASV table modified.",
      "No sample-depth cutoff applied."
    ),
    file.path(
      OUT,
      "_STEP85B_COMPLETE.ok"
    )
  )

  ck("STEP85B2 COMPLETE")

  cat("\n============================================================\n")
  cat("SEPSIS V2 - STEP85B2 COMPLETE\n")
  cat("============================================================\n\n")

  print(
    project_summary,
    n=Inf,
    width=Inf
  )

  cat("\nDEPTH THRESHOLD AUDIT:\n")

  print(
    thresholds,
    n=Inf,
    width=Inf
  )

  cat("\nOutput:\n")
  cat(
    OUT,
    "\n"
  )
}

tryCatch(
  main(),
  error=function(e) {

    msg <- c(
      paste0(
        "STEP85B2 FATAL ERROR: ",
        Sys.time()
      ),

      paste0(
        "Message: ",
        conditionMessage(e)
      ),

      paste0(
        "Call: ",
        paste(
          deparse(
            conditionCall(e)
          ),
          collapse=" "
        )
      )
    )

    writeLines(
      msg,
      ERR,
      useBytes=TRUE
    )

    ck(
      "STEP85B2 FAILED"
    )

    message(
      paste(
        msg,
        collapse="\n"
      )
    )

    quit(
      save="no",
      status=1,
      runLast=FALSE
    )
  }
)
