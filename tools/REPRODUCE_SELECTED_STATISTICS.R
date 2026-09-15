#!/usr/bin/env Rscript
# Portable statistical subset. No package installation, network, FASTQ or sequence processing.
options(stringsAsFactors=FALSE)
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2L) stop('Usage: Rscript REPRODUCE_SELECTED_STATISTICS.R INPUT_DIRECTORY NEW_OUTPUT_DIRECTORY')
inputs <- normalizePath(args[1],winslash='/',mustWork=TRUE)
out <- args[2]
if (dir.exists(out) || file.exists(out)) stop('Refusing to overwrite output: ',out)
dir.create(out,recursive=TRUE)
read_input <- function(name) read.csv(file.path(inputs,name),check.names=FALSE,fileEncoding='UTF-8-BOM')
flag <- function(x) toupper(as.character(x)) %in% c('TRUE','T','1','YES')
summarize <- function(x,id,continuity=TRUE) {
  x <- x[is.finite(x)]
  if(length(x)<2L || sd(x)==0) stop('Insufficient varying observations: ',id)
  n <- length(x); m <- mean(x); s <- sd(x); se <- s/sqrt(n); dz <- m/s
  J <- 1-3/(4*(n-1)-1)
  gz <- J*dz; gzse <- J*sqrt(1/n+dz^2/(2*n))
  ci <- m+c(-1,1)*qt(.975,n-1)*se
  w <- suppressWarnings(wilcox.test(x,mu=0,exact=FALSE,correct=continuity))
  data.frame(id=id,n=n,mean=m,sd=s,se=se,ci_low=ci[1],ci_high=ci[2],cohen_dz=dz,hedges_gz=gz,hedges_gz_se=gzse,hedges_ci_low=gz-qt(.975,n-1)*gzse,hedges_ci_high=gz+qt(.975,n-1)*gzse,t_p=2*pt(-abs(m/se),n-1),wilcoxon_p=unname(w$p.value))
}
meta <- function(tab,id) {
  yi <- tab$mean; vi <- tab$sd^2/tab$n; k <- nrow(tab)
  fn <- function(tau) {
    w <- 1/(vi+tau); m <- sum(w*yi)/sum(w)
    .5*(sum(log(vi+tau))+log(sum(w))+sum(w*(yi-m)^2))
  }
  upper <- max(1,var(yi)*100,max(vi)*100)
  fit <- optimize(fn,c(0,upper),tol=1e-12)
  tau <- if(fn(0)<=fit$objective+1e-10) 0 else fit$minimum
  w <- 1/(vi+tau); m <- sum(w*yi)/sum(w)
  conventional <- sqrt(1/sum(w)); hk <- sqrt(sum(w*(yi-m)^2)/(k-1)/sum(w))
  se <- max(conventional,hk); ci <- m+c(-1,1)*qt(.975,k-1)*se
  wf <- 1/vi; mf <- sum(wf*yi)/sum(wf); Q <- sum(wf*(yi-mf)^2)
  data.frame(id=id,k=k,pooled_mean=m,tau2=tau,se=se,ci_low=ci[1],ci_high=ci[2],p_value=2*pt(-abs(m/se),k-1),I2=100*max(0,(Q-(k-1))/Q))
}
checks <- list()
check_num <- function(id,observed,expected,tolerance=2e-8) {
  delta <- abs(observed-expected)
  pass <- is.finite(delta) && delta<=tolerance
  checks[[length(checks)+1L]] <<- data.frame(check=id,observed=observed,expected=expected,absolute_difference=delta,tolerance=tolerance,status=if(pass)'PASS' else 'FAIL')
}

# Reconstruct external paired distance contrasts by keys, rather than reading stored pair effects.
d <- read_input('external_sample_displacement.csv')
p <- read_input('external_population_membership.csv')
ref <- read_input('external_frozen_results.csv')
if (any(duplicated(paste(d$patient_id,d$timepoint,sep='::')))) stop('Duplicate external patient/time key')
distance_results <- list()
for (pop in c('PRESPECIFIED_PRIMARY_DATE_CLEAN','COMPLETE_CASE_SENSITIVITY')) {
  field <- if(pop=='PRESPECIFIED_PRIMARY_DATE_CLEAN') 'in_prespecified_primary' else 'in_complete_case_sensitivity'
  ids <- unique(p$patient_id[flag(p[[field]])])
  for(metric in c('bray_from_T0','aitchison_from_T0_CZM')) {
    t1 <- d[d$timepoint=='T1' & d$patient_id %in% ids,c('patient_id',metric)]
    t2 <- d[d$timepoint=='T2' & d$patient_id %in% ids,c('patient_id',metric)]
    paired <- merge(t1,t2,by='patient_id',suffixes=c('_early','_late'),all=TRUE)
    if(nrow(paired)!=length(ids) || any(!is.finite(paired[[2]]) | !is.finite(paired[[3]]))) stop('Incomplete external paired inputs')
    res <- summarize(paired[[3]]-paired[[2]],paste(pop,metric,sep='::'),continuity=FALSE)
    distance_results[[length(distance_results)+1L]] <- res
    row <- ref[ref$analysis_id==pop & ref$metric==metric,]
    if(nrow(row)!=1L) stop('Missing frozen external result ',pop,' ',metric)
    mapping <- c(n='n_patients',mean='mean_paired_difference',sd='sd_paired_difference',ci_low='mean_paired_difference_ci_low',ci_high='mean_paired_difference_ci_high',cohen_dz='cohen_dz',hedges_gz='hedges_gz',hedges_gz_se='hedges_gz_se',hedges_ci_low='hedges_gz_ci_low',hedges_ci_high='hedges_gz_ci_high',t_p='paired_t_p_value',wilcoxon_p='wilcoxon_p_value')
    for(field in names(mapping)) check_num(paste(res$id,field,sep='::'),res[[field]],row[[mapping[field]]])
  }
}
write.csv(do.call(rbind,distance_results),file.path(out,'external_distance_statistics.csv'),row.names=FALSE)

# Recalculate balance from raw frozen family sums; retain original pair population and double-zero exclusions.
pairs <- read_input('direction_pairs.csv')
if(any(duplicated(paste(pairs$project,pairs$patient_id)))) stop('Duplicate direction patient key')
balance <- function(P,C,pc=.5) ifelse(is.finite(P) & is.finite(C) & P+C>0,log((P+pc)/(C+pc)),NA_real_)
early <- balance(pairs$early_P_reads,pairs$early_C_reads)
late <- balance(pairs$late_P_reads,pairs$late_C_reads)
pairs$recomputed_delta <- late-early
balance_results <- list(); bref <- read_input('balance_frozen_cohort_effects.csv')
for(project in c('PRJNA691455','PRJNA851469','PRJNA516701','PRJNA1125274')) {
  eligible <- pairs$project==project & flag(pairs$in_frozen_primary_pair_population)
  res <- summarize(pairs$recomputed_delta[eligible],project)
  balance_results[[project]] <- res
  row <- bref[bref$project==project,]
  for(field in c('n','mean','sd','se','ci_low','ci_high','cohen_dz','t_p','wilcoxon_p')) {
    ref_field <- if(field=='n') 'n_analyzed' else field
    check_num(paste('BALANCE',project,field,sep='::'),res[[field]],row[[ref_field]])
  }
}
write.csv(do.call(rbind,balance_results),file.path(out,'direction_balance_statistics.csv'),row.names=FALSE)
pooled <- meta(do.call(rbind,balance_results[c('PRJNA691455','PRJNA851469','PRJNA516701')]),'NATURAL_HISTORY_REML_HK')
mref <- read_input('balance_frozen_meta.csv'); row <- mref[mref$analysis=='PRIMARY_PC0_5',]
for(field in c('pooled_mean','tau2','ci_low','ci_high','p_value','I2')) {
  ref_field <- c(pooled_mean='pooled_mean',tau2='tau2_REML',ci_low='ci_low',ci_high='ci_high',p_value='p_value',I2='I2_percent')[field]
  check_num(paste('META',field,sep='::'),pooled[[field]],row[[ref_field]])
}
primary <- data.frame(test_id=c('NATURAL_HISTORY_REML_HK','PRJNA1125274_EXTERNAL_24'),raw_p=c(pooled$p_value,balance_results[['PRJNA1125274']]$t_p))
primary$holm_p <- p.adjust(primary$raw_p,method='holm')
href <- read_input('balance_frozen_holm.csv')
for(i in seq_len(nrow(primary))) check_num(primary$test_id[i],primary$holm_p[i],href$holm_p[match(primary$test_id[i],href$test_id)],0)
write.csv(pooled,file.path(out,'balance_meta_analysis.csv'),row.names=FALSE)
write.csv(primary,file.path(out,'balance_holm_family.csv'),row.names=FALSE)

# Ordinary same-cohort healthy-reference contrast; the two-level bootstrap is not rerun by this subset.
eligible <- pairs$project=='PRJNA851469' & flag(pairs$in_frozen_primary_pair_population)
hres <- summarize(pairs$late_H[eligible]-pairs$early_H[eligible],'PRJNA851469_HEALTHY_REFERENCE_ORDINARY')
h <- read_input('healthy_frozen_result.csv')
for(field in c('n','mean','sd','ci_low','ci_high','t_p')) check_num(paste('HEALTHY',field,sep='::'),hres[[field]],h[[if(field=='n') 'n_analyzed' else field]])
write.csv(hres,file.path(out,'healthy_reference_ordinary_statistics.csv'),row.names=FALSE)
checks <- do.call(rbind,checks)
write.csv(checks,file.path(out,'reproduction_checks.csv'),row.names=FALSE)
writeLines(c('Selected-statistics scope only. Distances and family sums are frozen upstream inputs.','No ASV-distance reconstruction, sequence processing, clinical crosswalk or two-level bootstrap is performed.',capture.output(sessionInfo())),file.path(out,'session_info_and_scope.txt'))
if(any(checks$status!='PASS')) stop('Reproduction failed: see reproduction_checks.csv')
cat(nrow(checks),' numerical checks PASS\n')
