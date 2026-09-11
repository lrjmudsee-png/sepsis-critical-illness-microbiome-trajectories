# Read existing 98D count object and 98E displacement only; no sequence processing.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr);library(tibble)})
root <- 'E:/sepsis_project'
up <- file.path(root,'results/V2_UPGRADE_20260907')
out <- file.path(up,'98G_PRJNA1125274_PATIENT_RECONCILIATION_20260909')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
source98e <- file.path(root,'code/03_data_processing/98E_PRJNA1125274_external_validation.R')
# Evaluate only the two statistical function definitions, never the script body.
ast <- parse(source98e)
for(e in ast) {
  if(is.call(e) && as.character(e[[1]]) %in% c('<-','=') &&
     is.symbol(e[[2]]) && as.character(e[[2]]) %in% c('J_paired','primary_contrast')) eval(e)
}
stopifnot(exists('primary_contrast'),exists('J_paired'))
obj <- readRDS(file.path(root,'data/PRJNA1125274/03_dada2/PRJNA1125274_analysis_object_external_validation.rds'))
meta <- as.data.frame(obj$metadata)
meta$run_id <- as.character(meta$run_accession)
meta$depth <- rowSums(obj$counts)[match(meta$run_id,rownames(obj$counts))]
meta$passes_depth <- meta$depth >= 2000
write_excel_csv(meta[,c('run_id','patient_id','hospital','timepoint','depth','passes_depth')],
                file.path(out,'analysis_object_run_accounting.csv'))

# Orphans are not in paired tests, but their counts were used before feature filtering.
fw <- meta[meta$passes_depth & meta$timepoint %in% c('T0','T1','T2'),]
cnt <- obj$counts[fw$run_id,,drop=FALSE]
global_keep <- colSums(cnt)>=20 & colSums(cnt>0)>=3
orphans <- c('SNO_4','SNO_34')
cnt2 <- cnt[!fw$patient_id %in% orphans,,drop=FALSE]
restricted_keep <- colSums(cnt2)>=20 & colSums(cnt2>0)>=3
filter_audit <- data.frame(
  samples_in_object=nrow(obj$counts),asvs_in_object=ncol(obj$counts),
  samples_passing_depth=nrow(fw),orphan_depth_passing_samples=sum(fw$patient_id %in% orphans),
  asvs_original_filter=sum(global_keep),
  asvs_after_hypothetical_orphan_removal=sum(restricted_keep),
  asvs_filter_membership_changed=sum(global_keep!=restricted_keep),
  removal_applied=FALSE,
  interpretation='Orphan records do not enter paired tests, but participated in global ASV filtering. No removal or identity correction was applied.')
write_excel_csv(filter_audit,file.path(out,'PRJNA1125274_orphan_feature_filter_influence.csv'))

efile <- file.path(up,'98E_PRJNA1125274_EXTERNAL_VALIDATION')
disp <- read_csv(file.path(efile,'PRJNA1125274_external_validation_sample_displacement.csv'),show_col_types=FALSE)
old <- read_csv(file.path(efile,'PRJNA1125274_external_validation_primary_contrast.csv'),show_col_types=FALSE)
comparison <- list()
diffs <- list()
for(i in seq_len(nrow(old))) {
  row <- old[i,]
  is_sens <- row$sensitivity_set=='EXCLUDING_DATE_ANOMALY'
  recalc <- primary_contrast(disp,row$metric,is_sens)
  numeric_fields <- names(row)[vapply(row,is.numeric,logical(1))]
  errs <- vapply(numeric_fields,function(f) abs(row[[f]]-recalc[[f]]),numeric(1))
  stopifnot(all(errs < 1e-10,na.rm=TRUE))
  sub <- disp[disp$timepoint %in% c('T1','T2'),]
  if(is_sens) sub<-sub[sub$mapping_status!='LETTER_ORDER_TIME_CODE_DATE_ANOMALY',]
  pair <- sub |> select(patient_id,timepoint,val=all_of(row$metric)) |>
    distinct(patient_id,timepoint,.keep_all=TRUE) |>
    pivot_wider(names_from=timepoint,values_from=val) |>
    filter(!is.na(T1),!is.na(T2))
  stopifnot(!any(pair$patient_id %in% orphans))
  d <- pair$T2-pair$T1
  pair$paired_difference <- d
  pair$metric <- row$metric
  pair$sensitivity_set <- row$sensitivity_set
  diffs[[i]] <- pair
  ci <- mean(d)+c(-1,1)*qt(.975,length(d)-1)*sd(d)/sqrt(length(d))
  val <- tibble(metric=row$metric,sensitivity_set=row$sensitivity_set,
    comparison_type='BEFORE_VS_RECONCILIATION_NO_IDENTITY_CORRECTION',
    mapping_changed=FALSE,sequence_processing_rerun=FALSE,
    verification='98E_primary_contrast_reexecuted_on_frozen_displacement',
    n_complete_before=row$n_patients,n_complete_after=recalc$n_patients,
    mean_contrast_before=row$mean_paired_diff,mean_contrast_after=recalc$mean_paired_diff,
    mean_contrast_ci_low_before=ci[1],mean_contrast_ci_high_before=ci[2],
    mean_contrast_ci_low_after=ci[1],mean_contrast_ci_high_after=ci[2],
    cohen_dz_before=row$cohen_dz,cohen_dz_after=recalc$cohen_dz,
    hedges_gz_before=row$hedges_gz,hedges_gz_after=recalc$hedges_gz,
    gz_ci_low_before=row$gz_ci_low,gz_ci_high_before=row$gz_ci_high,
    gz_ci_low_after=recalc$gz_ci_low,gz_ci_high_after=recalc$gz_ci_high,
    p_value_before=row$paired_t_p,p_value_after=recalc$paired_t_p,
    wilcoxon_p_before=row$wilcoxon_p,wilcoxon_p_after=recalc$wilcoxon_p,
    delta_n=recalc$n_patients-row$n_patients,
    delta_mean=recalc$mean_paired_diff-row$mean_paired_diff,
    delta_hedges_gz=recalc$hedges_gz-row$hedges_gz,
    delta_p=recalc$paired_t_p-row$paired_t_p,
    max_abs_difference_from_saved_98E=max(errs,na.rm=TRUE),
    note='No corrected mapping is justified. After means verified unchanged mapping, not a fabricated 132-patient remapping. Mean-difference CI computed from frozen differences; gz CI uses original 98E approximation.')
  comparison[[i]] <- val
}
write_excel_csv(bind_rows(comparison),file.path(out,'PRJNA1125274_primary_results_before_after_correction.csv'))
write_excel_csv(bind_rows(diffs),file.path(out,'PRJNA1125274_primary_patient_differences_verified.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'98G2_R_sessionInfo.txt'))
print(filter_audit)
print(bind_rows(comparison) |> select(metric,sensitivity_set,n_complete_before,mean_contrast_after,hedges_gz_after,p_value_after,max_abs_difference_from_saved_98E))
