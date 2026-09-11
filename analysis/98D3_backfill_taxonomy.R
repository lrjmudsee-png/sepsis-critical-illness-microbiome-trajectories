# Backfill SILVA 138.2 taxonomy for the SURVEIL single-end seqtab in batches
suppressPackageStartupMessages({
  library(dada2); library(dplyr); library(tibble); library(readr)
})
st <- readRDS("E:/sepsis_project/data/PRJNA1125274/03_dada2/PRJNA1125274_full_seqtab_se.rds")
silva_train <- "E:/sepsis_project/data/_V2_ANALYSIS_READY/00_REFERENCE/SILVA_138.2_DADA2_SSU/silva_nr99_v138.2_toGenus_trainset.fa.gz"
seqs <- colnames(st)
cat("total ASVs:", length(seqs), "\n")
oklen <- nchar(seqs) >= 50
cat("ASVs with length>=50:", sum(oklen), "\n")

tax <- matrix(NA_character_, nrow = length(seqs), ncol = 6,
              dimnames = list(NULL, c("Kingdom","Phylum","Class","Order","Family","Genus")))
tax <- as.data.frame(tax, stringsAsFactors = FALSE)
BATCH <- 15000
idx <- which(oklen)
nb <- ceiling(length(idx) / BATCH)
cat("batches:", nb, "\n")
for (b in seq_len(nb)) {
  i <- idx[((b - 1) * BATCH + 1):min(b * BATCH, length(idx))]
  cat("batch", b, "n =", length(i), "\n")
  ta <- dada2::assignTaxonomy(seqs[i], refFasta = silva_train,
                              multithread = 1, tryRC = TRUE, verbose = FALSE)
  tax[i, colnames(ta)] <- ta
  gc()
}
cat("taxonomy assignment complete\n")
saveRDS(tax, "E:/sepsis_project/data/PRJNA1125274/03_dada2/PRJNA1125274_taxa_final_se.rds")
cat("taxa rds saved\n")

# update analysis object
obj_path <- "E:/sepsis_project/data/PRJNA1125274/03_dada2/PRJNA1125274_analysis_object_external_validation.rds"
obj <- readRDS(obj_path)
asv_tax <- tax |>
  mutate(ASV_ID = paste0("ASV_", seq_len(nrow(tax))),
         ASV_sequence = seqs) |>
  select(ASV_ID, ASV_sequence, everything())
obj$taxonomy <- as_tibble(asv_tax)
obj$provenance$taxonomy_state <- "COMPLETE SILVA138.2 (batched assignTaxonomy)"
obj$provenance$taxonomy_created <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
saveRDS(obj, obj_path)
cat("analysis object updated with taxonomy\n")
