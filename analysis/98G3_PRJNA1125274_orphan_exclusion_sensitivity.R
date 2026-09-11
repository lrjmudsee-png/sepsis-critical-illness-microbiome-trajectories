# Exploratory sensitivity, not a corrected patient map.
# Re-execute ONLY 98E statistics from existing 98D object with 3 orphan records omitted.
# Does not read FASTQ, execute 98D, or rewrite the count/taxonomy object.
script <- 'E:/sepsis_project/code/03_data_processing/98E_PRJNA1125274_external_validation.R'
recon_out <- 'E:/sepsis_project/results/V2_UPGRADE_20260907/98G_PRJNA1125274_PATIENT_RECONCILIATION_20260909'
sensitivity_out <- file.path(recon_out,'98E_ORPHAN_EXCLUSION_SENSITIVITY')
ast <- parse(script)
for(e in ast) {
  assignment <- is.call(e) && as.character(e[[1]]) %in% c('<-','=') && is.symbol(e[[2]])
  name <- if(assignment) as.character(e[[2]]) else ''
  if(name=='OUT') {
    OUT <- sensitivity_out
  } else {
    eval(e, envir=.GlobalEnv)
    if(name=='meta') {
      meta <- meta[!meta$patient_id %in% c('SNO_4','SNO_34'),,drop=FALSE]
      counts <- counts[rownames(counts) %in% meta$run_accession,,drop=FALSE]
    }
  }
}
old <- readr::read_csv(file.path(ROOT,'results/V2_UPGRADE_20260907/98E_PRJNA1125274_EXTERNAL_VALIDATION/PRJNA1125274_external_validation_primary_contrast.csv'),show_col_types=FALSE)
new <- readr::read_csv(file.path(OUT,'PRJNA1125274_external_validation_primary_contrast.csv'),show_col_types=FALSE)
comparison <- dplyr::inner_join(old,new,by=c('metric','sensitivity_set'),suffix=c('_original','_orphan_exclusion')) |>
  dplyr::mutate(analysis_label='EXPLORATORY_ORPHAN_EXCLUSION_NOT_IDENTITY_CORRECTION',
                excluded_patient_ids='SNO_4;SNO_34',excluded_runs='SRR29453250;SRR29453247;SRR29453297',
                delta_n=n_patients_orphan_exclusion-n_patients_original,
                delta_mean=mean_paired_diff_orphan_exclusion-mean_paired_diff_original,
                delta_hedges_gz=hedges_gz_orphan_exclusion-hedges_gz_original,
                delta_p=paired_t_p_orphan_exclusion-paired_t_p_original,
                sequence_processing_rerun=FALSE)
readr::write_excel_csv(comparison,file.path(recon_out,'PRJNA1125274_orphan_exclusion_sensitivity_comparison.csv'))
writeLines(capture.output(sessionInfo()),file.path(OUT,'sessionInfo.txt'))
print(comparison |> dplyr::select(metric,sensitivity_set,n_patients_orphan_exclusion,hedges_gz_orphan_exclusion,paired_t_p_orphan_exclusion,delta_mean))
