#!/usr/bin/env Rscript
# Timestamped post-hoc specimen governance audit, 2026-09-15. Frozen primary sets are not redefined.
options(stringsAsFactors=FALSE)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2) stop('Usage: Rscript REPRODUCE_SPECIMEN_SENSITIVITY.R INPUT_DIRECTORY NEW_OUTPUT_DIRECTORY')
if(!requireNamespace('metafor',quietly=TRUE)) stop('Install the documented metafor environment before running; no automatic installation is performed')
input <- normalizePath(args[1],mustWork=TRUE,winslash='/'); out <- args[2]
if(file.exists(out)) stop('Refusing to overwrite output')
dir.create(out,recursive=TRUE)
read <- function(name) read.csv(file.path(input,name),fileEncoding='UTF-8-BOM',check.names=FALSE)
flag <- function(x) toupper(as.character(x)) %in% c('TRUE','T','1')
spec <- read('specimen_primary_pairs.csv'); pd <- read('discovery_paired_differences.csv'); bp <- read('direction_pairs.csv')
stopifnot(nrow(spec)==15,!anyDuplicated(spec$patient_id))
if(anyDuplicated(paste(pd$project,pd$patient_id,pd$scheme))) stop('Duplicate discovery patient/geometry key')
summ <- function(x,project,scheme,subset,continuity=FALSE) {
  x <- x[is.finite(x)]; n <- length(x); m <- mean(x); s <- sd(x); se <- s/sqrt(n)
  if(n<2 || s==0) stop('Insufficient varying observations')
  dz <- m/s; J <- 1-3/(4*(n-1)-1); g <- J*dz; gse <- J*sqrt(1/n+dz^2/(2*n)); crit <- qt(.975,n-1)
  data.frame(project=project,scheme=scheme,subset=subset,n=n,mean=m,sd=s,se=se,ci_low=m-crit*se,ci_high=m+crit*se,hedges_gz=g,hedges_gz_se=gse,hedges_ci_low=g-crit*gse,hedges_ci_high=g+crit*gse,t_p=2*pt(-abs(m/se),n-1),wilcoxon_p=suppressWarnings(wilcox.test(x,exact=FALSE,correct=continuity)$p.value))
}
sets <- list(FROZEN_ALL_15=spec$patient_id,
 SAME_EARLY_LATE_11=spec$patient_id[flag(spec$early_late_same_specimen)],
 SAME_ALL_AVAILABLE_10=spec$patient_id[flag(spec$all_available_same_specimen)],
 SAME_ALL_WITH_EXPLICIT_ANCHOR_9=spec$patient_id[flag(spec$all_available_same_specimen)&flag(spec$anchor_available)],
 STOOL_ALL_AVAILABLE_5=spec$patient_id[flag(spec$all_available_stool)],
 RECTAL_SWAB_ALL_AVAILABLE_5=spec$patient_id[flag(spec$all_available_same_specimen)&spec$early_specimen=='RECTAL_SWAB'])
cohort <- list(); pools <- list(); contrasts <- list()
for(scheme in c('bray','aitchison_CZM')) {
  fixed <- lapply(c('PRJNA691455','PRJNA851469'),function(project) summ(pd$paired_diff[pd$project==project&pd$scheme==scheme],project,scheme,'FROZEN_UNCHANGED'))
  for(label in names(sets)) {
    x <- pd[pd$project=='PRJNA516701' & pd$scheme==scheme & pd$patient_id %in% sets[[label]],]
    stopifnot(nrow(x)==length(sets[[label]]))
    s <- summ(x$paired_diff,'PRJNA516701',scheme,label); cohort[[length(cohort)+1]] <- s
    if(label %in% c('FROZEN_ALL_15','SAME_ALL_AVAILABLE_10','SAME_ALL_WITH_EXPLICIT_ANCHOR_9')) {
      ce <- do.call(rbind,c(fixed,list(s)))
      fit <- metafor::rma(yi=ce$hedges_gz,sei=ce$hedges_gz_se,method='REML',test='knha')
      pools[[length(pools)+1]] <- data.frame(scheme=scheme,subset=label,k=fit$k,n_PRJNA516701=s$n,pooled_hedges_gz=as.numeric(fit$b),ci_low=fit$ci.lb,ci_high=fit$ci.ub,p_value=fit$pval,tau2=fit$tau2,I2_percent=fit$I2)
    }
    contrasts[[length(contrasts)+1]] <- data.frame(project=x$project,patient_id=x$patient_id,scheme=x$scheme,subset=label,paired_diff=x$paired_diff)
  }
  ce <- do.call(rbind,fixed); fit <- metafor::rma(yi=ce$hedges_gz,sei=ce$hedges_gz_se,method='REML',test='knha')
  pools[[length(pools)+1]] <- data.frame(scheme=scheme,subset='EXCLUDE_PRJNA516701_K2',k=fit$k,n_PRJNA516701=0,pooled_hedges_gz=as.numeric(fit$b),ci_low=fit$ci.lb,ci_high=fit$ci.ub,p_value=fit$pval,tau2=fit$tau2,I2_percent=fit$I2)
}
balance <- function(P,C) ifelse(is.finite(P)&is.finite(C)&P+C>0,log((P+.5)/(C+.5)),NA_real_)
b <- bp[bp$project=='PRJNA516701'&flag(bp$in_frozen_primary_pair_population),]
b$delta <- balance(b$late_P_reads,b$late_C_reads)-balance(b$early_P_reads,b$early_C_reads)
for(label in names(sets)) {
  s <- summ(b$delta[b$patient_id %in% sets[[label]]],'PRJNA516701','SECONDARY_FAMILY_BALANCE',label,TRUE)
  cohort[[length(cohort)+1]] <- s
}
write.csv(do.call(rbind,cohort),file.path(out,'specimen_cohort_sensitivity_statistics.csv'),row.names=FALSE)
write.csv(do.call(rbind,pools),file.path(out,'specimen_meta_sensitivity_statistics.csv'),row.names=FALSE)
write.csv(do.call(rbind,contrasts),file.path(out,'specimen_sensitivity_patient_differences_LOCAL.csv'),row.names=FALSE)
writeLines(c('Post-hoc governance sensitivity; same fixed distances and count sums; no sequence processing.',
 'SAME_ALL_AVAILABLE includes the frozen no-explicit-anchor patient; SAME_ALL_WITH_EXPLICIT_ANCHOR excludes this patient.',
 'All sensitivity p values and intervals are nominal, not new primary tests. Low sample sizes and swab soiling remain limitations.',capture.output(sessionInfo())),file.path(out,'specimen_sensitivity_session_info.txt'))
cat('Specimen sensitivity complete; original populations preserved.\n')
