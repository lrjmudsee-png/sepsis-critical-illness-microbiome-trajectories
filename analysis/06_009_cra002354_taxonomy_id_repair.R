# ============================================================
# Sepsis V2 - Step94B1E
# CRA002354 SILVA TAXONOMY OTU-ID REPAIR
#
# B1D showed 0% abundance-weighted taxonomy coverage, despite B1B
# reporting >10k OTUs with genus assignment.
#
# Likely cause:
# dada2::assignTaxonomy() output row names are sequence strings rather
# than OTU_1 / OTU_2 identifiers. The raw taxonomy classifications are
# retained, but downstream joins by otu_id fail.
#
# This script repairs taxonomy IDs by exact centroid-sequence matching,
# audits the mapping, recalculates abundance-weighted taxonomy coverage,
# and rebuilds corrected genus abundance.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Bioconductor package 'Biostrings' is required.")
}

ROOT <- "E:/sepsis_project"

DATA_DIR <- file.path(
  ROOT, "data", "CRA002354", "03_vsearch97_silva1382"
)

OUT <- file.path(
  ROOT, "results",
  "V2_34B1E_CRA002354_TAXONOMY_ID_REPAIR"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

CENTROIDS <- file.path(
  DATA_DIR,
  "CRA002354_OTU97_centroids.fa"
)

TAX_RAW <- file.path(
  DATA_DIR,
  "CRA002354_OTU97_SILVA1382_taxonomy_minBoot80.csv"
)

COUNTS <- file.path(
  DATA_DIR,
  "CRA002354_OTU97_counts_NONCHIMERIC_CORRECTED.csv"
)

needed <- c(CENTROIDS, TAX_RAW, COUNTS)
missing <- needed[!file.exists(needed)]

if (length(missing) > 0) {
  stop(paste("Missing required file(s):", paste(missing, collapse = "; ")))
}

dna <- Biostrings::readDNAStringSet(CENTROIDS)

centroid <- tibble(
  otu_id_fixed = sub(";.*$", "", names(dna)),
  centroid_sequence = toupper(as.character(dna))
)

if (anyDuplicated(centroid$otu_id_fixed)) {
  stop("Duplicated centroid OTU IDs detected.")
}

tax <- read_csv(
  TAX_RAW,
  show_col_types = FALSE,
  name_repair = "unique"
)

if (!"otu_id" %in% names(tax)) {
  stop("Raw taxonomy file has no otu_id column.")
}

tax_key <- toupper(as.character(tax$otu_id))

id_overlap <- mean(tax_key %in% toupper(centroid$otu_id_fixed))
seq_overlap <- mean(tax_key %in% centroid$centroid_sequence)

mapping_mode <- NA_character_

if (id_overlap >= 0.95) {

  mapping_mode <- "OTU_ID_ALREADY_VALID"

  fixed_ids <- centroid$otu_id_fixed[
    match(
      toupper(tax$otu_id),
      toupper(centroid$otu_id_fixed)
    )
  ]

} else if (seq_overlap >= 0.95) {

  mapping_mode <- "REPAIRED_BY_EXACT_CENTROID_SEQUENCE"

  fixed_ids <- centroid$otu_id_fixed[
    match(
      toupper(tax$otu_id),
      centroid$centroid_sequence
    )
  ]

} else {

  audit <- tibble(
    taxonomy_rows = nrow(tax),
    centroid_rows = nrow(centroid),
    direct_OTU_ID_overlap = id_overlap,
    exact_sequence_overlap = seq_overlap,
    status = "UNRESOLVED_DO_NOT_USE_ROW_ORDER_FALLBACK"
  )

  write_csv(
    audit,
    file.path(OUT, "01_TAXONOMY_ID_MAPPING_AUDIT.csv")
  )

  stop(
    paste0(
      "Taxonomy ID mapping unresolved. direct ID overlap=",
      signif(id_overlap,4),
      "; sequence overlap=",
      signif(seq_overlap,4),
      ". No unsafe row-order fallback was used."
    )
  )
}

tax_fixed <- tax
tax_fixed$otu_id_raw <- tax_fixed$otu_id
tax_fixed$otu_id <- fixed_ids

mapped_fraction <- mean(!is.na(tax_fixed$otu_id))

audit <- tibble(
  taxonomy_rows = nrow(tax),
  centroid_rows = nrow(centroid),
  direct_OTU_ID_overlap_before = id_overlap,
  exact_sequence_overlap_before = seq_overlap,
  mapping_mode = mapping_mode,
  mapped_fraction_after = mapped_fraction,
  duplicated_fixed_OTU_IDs =
    sum(duplicated(tax_fixed$otu_id[!is.na(tax_fixed$otu_id)]))
)

write_csv(
  audit,
  file.path(OUT, "01_TAXONOMY_ID_MAPPING_AUDIT.csv")
)

if (mapped_fraction < 0.95) {
  stop("Less than 95% taxonomy rows mapped to OTU IDs after repair.")
}

write_csv(
  tax_fixed,
  file.path(
    DATA_DIR,
    "CRA002354_OTU97_SILVA1382_taxonomy_minBoot80_ID_FIXED.csv"
  )
)

# ------------------------------------------------------------
# Corrected abundance-weighted coverage
# ------------------------------------------------------------

cnt <- read_csv(
  COUNTS,
  show_col_types = FALSE,
  name_repair = "minimal"
)

otu_cols <- setdiff(names(cnt), "Run_ID")

M <- as.matrix(cnt[, otu_cols, drop = FALSE])
storage.mode(M) <- "numeric"

otu_totals <- tibble(
  otu_id = colnames(M),
  abundance = as.numeric(colSums(M))
)

tax_ab <- otu_totals %>%
  left_join(
    tax_fixed %>% select(-otu_id_raw),
    by = "otu_id"
  )

ranks <- intersect(
  c("Kingdom","Phylum","Class","Order","Family","Genus"),
  names(tax_ab)
)

coverage <- bind_rows(
  lapply(ranks, function(rk) {
    assigned <- !is.na(tax_ab[[rk]]) & tax_ab[[rk]] != ""

    tibble(
      rank = rk,
      OTUs_assigned = sum(assigned),
      OTUs_total = nrow(tax_ab),
      unweighted_OTU_fraction = mean(assigned),
      reads_in_assigned_OTUs = sum(tax_ab$abundance[assigned]),
      total_reads = sum(tax_ab$abundance),
      abundance_weighted_fraction =
        sum(tax_ab$abundance[assigned]) / sum(tax_ab$abundance)
    )
  })
)

write_csv(
  coverage,
  file.path(OUT, "02_FIXED_TAXONOMY_COVERAGE_WEIGHTED_BY_ABUNDANCE.csv")
)

# ------------------------------------------------------------
# Genus aggregation after ID repair
# ------------------------------------------------------------

make_genus_label <- function(row) {

  ranks0 <- c("Kingdom","Phylum","Class","Order","Family","Genus")
  ranks0 <- ranks0[ranks0 %in% names(row)]

  vals <- as.character(row[ranks0])
  names(vals) <- ranks0

  if (
    "Genus" %in% names(vals) &&
    !is.na(vals["Genus"]) &&
    nzchar(vals["Genus"])
  ) {
    return(vals["Genus"])
  }

  nonmissing <- vals[!is.na(vals) & nzchar(vals)]

  if (length(nonmissing) == 0) {
    return("Unclassified_Bacteria")
  }

  lr <- names(nonmissing)[length(nonmissing)]
  lv <- nonmissing[length(nonmissing)]

  paste0(
    "Unclassified_",
    lr,
    "_",
    gsub("[^A-Za-z0-9]+", "_", lv)
  )
}

tax_fixed$genus_label <- apply(
  tax_fixed,
  1,
  make_genus_label
)

genus_for_otu <- tax_fixed$genus_label[
  match(otu_cols, tax_fixed$otu_id)
]

genus_for_otu[
  is.na(genus_for_otu) | genus_for_otu == ""
] <- "Unclassified_Bacteria"

# sample x OTU -> genus x sample for rowsum
genus_counts_t <- rowsum(
  t(M),
  group = genus_for_otu,
  reorder = FALSE
)

genus_counts <- t(genus_counts_t)
genus_rel <- genus_counts / rowSums(genus_counts)

rownames(genus_counts) <- cnt$Run_ID
rownames(genus_rel) <- cnt$Run_ID

write_csv(
  as.data.frame(genus_counts, check.names = FALSE) %>%
    rownames_to_column("Run_ID"),
  file.path(
    DATA_DIR,
    "CRA002354_genus_counts_NONCHIMERIC_TAXONOMY_ID_FIXED.csv"
  )
)

write_csv(
  as.data.frame(genus_rel, check.names = FALSE) %>%
    rownames_to_column("Run_ID"),
  file.path(
    DATA_DIR,
    "CRA002354_genus_relative_abundance_NONCHIMERIC_TAXONOMY_ID_FIXED.csv"
  )
)

genus_cov <- coverage %>%
  filter(rank == "Genus")

ready <- tibble(
  mapping_mode = mapping_mode,
  taxonomy_rows_mapped_fraction = mapped_fraction,
  abundance_weighted_genus_assignment =
    ifelse(nrow(genus_cov) == 1, genus_cov$abundance_weighted_fraction, NA_real_),
  taxonomy_ready_for_exploratory_genus_analysis =
    (
      mapped_fraction >= 0.95 &&
      nrow(genus_cov) == 1 &&
      genus_cov$abundance_weighted_fraction >= 0.50
    )
)

write_csv(
  ready,
  file.path(OUT, "03_STEP94B1E_TAXONOMY_READINESS.csv")
)

writeLines(
  c(
    "STEP94B1E TAXONOMY ID REPAIR",
    "",
    paste0("Mapping mode: ", mapping_mode),
    paste0("Mapped fraction: ", mapped_fraction),
    paste0(
      "Abundance-weighted genus assignment: ",
      ifelse(nrow(genus_cov) == 1,
             signif(genus_cov$abundance_weighted_fraction, 5),
             "NA")
    ),
    paste0(
      "Ready for exploratory genus analysis: ",
      ready$taxonomy_ready_for_exploratory_genus_analysis
    ),
    "",
    "This taxonomy repair does not alter OTU ecology results."
  ),
  file.path(OUT, "04_STEP94B1E_INTERPRETATION.txt")
)

writeLines(
  c(
    paste0("Completed: ", Sys.time()),
    "STEP94B1E COMPLETE"
  ),
  file.path(OUT, "_STEP94B1E_COMPLETE.ok")
)

cat("STEP94B1E COMPLETE\n")
