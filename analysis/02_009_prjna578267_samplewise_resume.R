options(stringsAsFactors = FALSE)

pkgs <- c('readr','dplyr','stringr','tibble','purrr')
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if(length(miss)) install.packages(miss, repos='https://cloud.r-project.org')
if(!requireNamespace('dada2', quietly=TRUE)) {
  if(!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org')
  BiocManager::install('dada2', ask=FALSE, update=FALSE)
}

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tibble); library(purrr); library(dada2)
})

PROJECT <- 'PRJNA578267'
ROOT <- 'E:/sepsis_project'
CANON <- file.path(ROOT,'data','_V2_ANALYSIS_READY')
PWORK <- file.path(CANON,'01_PROJECTS',PROJECT,'05_work','step84B_full')
FILT <- file.path(PWORK,'filtered')
CPROOT <- file.path(PWORK,'samplewise_checkpoints')
CF <- file.path(CPROOT,'DADA_F'); CR <- file.path(CPROOT,'DADA_R'); CM <- file.path(CPROOT,'MERGE')
dir.create(CF,recursive=TRUE,showWarnings=FALSE); dir.create(CR,recursive=TRUE,showWarnings=FALSE); dir.create(CM,recursive=TRUE,showWarnings=FALSE)
OUT <- file.path(ROOT,'results','V2_24B3_PRJNA578267_SAMPLEWISE')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
LOG <- file.path(OUT,'_STEP84B3_runtime.log')
ERR <- file.path(OUT,'_STEP84B3_FATAL_ERROR.txt')
if(file.exists(ERR)) unlink(ERR)

logmsg <- function(x){ line <- paste0(Sys.time(),' | ',x); cat(line,'\n'); cat(line,'\n',file=LOG,append=TRUE) }
safe_csv <- function(p) suppressMessages(read_csv(p,show_col_types=FALSE,progress=FALSE,name_repair='unique'))
extract_run <- function(x) str_extract(basename(x),regex('(SRR|ERR|DRR|CRR)[0-9]+',ignore_case=TRUE))
valid_rds <- function(p){ if(!file.exists(p)) return(FALSE); tryCatch({readRDS(p); TRUE}, error=function(e) FALSE) }
atomic_saveRDS <- function(obj,path){ tmp <- paste0(path,'.tmp'); if(file.exists(tmp)) unlink(tmp); saveRDS(obj,tmp); if(file.exists(path)) unlink(path); if(!file.rename(tmp,path)) stop('Could not save checkpoint: ',path) }
getN <- function(x) sum(getUniques(x))

main <- function(){
  manifest <- safe_csv(file.path(CANON,'00_FREEZE','V2_FINAL_FROZEN_16S_sequence_manifest.csv')) |> filter(project==PROJECT)
  inv <- safe_csv(file.path(ROOT,'results','V2_23B_READ_QUALITY_PRIMER_AUDIT','V2_STEP83B_canonical_raw_inventory.csv')) |>
    filter(project==PROJECT,mate %in% c('R1','R2')) |> mutate(run_id=extract_run(file_name))
  if(nrow(manifest)!=207) stop('Expected 207 frozen Runs; got ',nrow(manifest))

  r1 <- inv |> filter(mate=='R1') |> select(run_id,R1=file_path) |> distinct(run_id,.keep_all=TRUE)
  r2 <- inv |> filter(mate=='R2') |> select(run_id,R2=file_path) |> distinct(run_id,.keep_all=TRUE)
  pairs <- inner_join(r1,r2,by='run_id') |> inner_join(manifest |> select(run_id,patient_id,sample_id,time_raw,time_day),by='run_id') |>
    arrange(patient_id,time_day,time_raw,run_id)
  if(nrow(pairs)!=207) stop('Expected 207 exact pairs; got ',nrow(pairs))

  pairs <- pairs |> mutate(
    filtF=file.path(FILT,paste0(run_id,'_R1_filt.fastq.gz')),
    filtR=file.path(FILT,paste0(run_id,'_R2_filt.fastq.gz'))
  )
  if(!all(file.exists(pairs$filtF)) || !all(file.exists(pairs$filtR))) stop('Filtered FASTQ set incomplete.')

  errF_path <- file.path(PWORK,paste0(PROJECT,'_full_errF.rds'))
  errR_path <- file.path(PWORK,paste0(PROJECT,'_full_errR.rds'))
  if(!valid_rds(errF_path) || !valid_rds(errR_path)) stop('Valid full-cohort error models missing.')
  errF <- readRDS(errF_path); errR <- readRDS(errR_path)
  write_excel_csv(pairs,file.path(OUT,'PRJNA578267_samplewise_run_order.csv'),na='')
  logmsg('INPUT GUARDS PASSED: 207 pairs + filtered FASTQ + error models')

  run_dir <- function(files,runs,err,cpdir,label){
    timing <- vector('list',length(runs))
    for(i in seq_along(runs)){
      run <- runs[i]; cp <- file.path(cpdir,paste0(run,'.rds'))
      if(valid_rds(cp)){
        logmsg(paste(label,i,'/207',run,'CHECKPOINT EXISTS -> SKIP'))
        timing[[i]] <- tibble(direction=label,sample_index=i,run_id=run,status='REUSED',elapsed_sec=0)
        next
      }
      logmsg(paste(label,i,'/207',run,'START'))
      t0 <- proc.time()[['elapsed']]
      derep <- derepFastq(files[i],verbose=FALSE); names(derep) <- run
      dobj <- dada(derep,err=err,pool=FALSE,multithread=FALSE)
      elapsed <- proc.time()[['elapsed']] - t0
      atomic_saveRDS(dobj,cp)
      logmsg(paste(label,i,'/207',run,'COMPLETE',round(elapsed,1),'sec'))
      timing[[i]] <- tibble(direction=label,sample_index=i,run_id=run,status='COMPUTED',elapsed_sec=round(elapsed,3))
      rm(derep,dobj); gc(verbose=FALSE)
    }
    bind_rows(timing)
  }

  tf <- run_dir(pairs$filtF,pairs$run_id,errF,CF,'FORWARD')
  write_excel_csv(tf,file.path(OUT,'PRJNA578267_forward_sample_timing.csv'),na='')
  tr <- run_dir(pairs$filtR,pairs$run_id,errR,CR,'REVERSE')
  write_excel_csv(tr,file.path(OUT,'PRJNA578267_reverse_sample_timing.csv'),na='')

  logmsg('SAMPLE-WISE MERGE START')
  stats <- vector('list',207)
  for(i in seq_len(207)){
    run <- pairs$run_id[i]; cp <- file.path(CM,paste0(run,'.rds'))
    if(valid_rds(cp)){
      merger <- readRDS(cp); logmsg(paste('MERGE',i,'/207',run,'CHECKPOINT EXISTS -> SKIP'))
    } else {
      logmsg(paste('MERGE',i,'/207',run,'START'))
      dF <- readRDS(file.path(CF,paste0(run,'.rds'))); dR <- readRDS(file.path(CR,paste0(run,'.rds')))
      derepF <- derepFastq(pairs$filtF[i],verbose=FALSE); derepR <- derepFastq(pairs$filtR[i],verbose=FALSE)
      names(derepF) <- run; names(derepR) <- run
      merger <- mergePairs(dF,derepF,dR,derepR,minOverlap=12,maxMismatch=0,verbose=FALSE)
      atomic_saveRDS(merger,cp)
      logmsg(paste('MERGE',i,'/207',run,'COMPLETE'))
      rm(dF,dR,derepF,derepR); gc(verbose=FALSE)
    }
    filtF_count <- sum(getUniques(derepFastq(pairs$filtF[i],verbose=FALSE)))
    merged_count <- if(is.null(merger)||!nrow(merger)) 0 else sum(merger$abundance,na.rm=TRUE)
    stats[[i]] <- tibble(sample_index=i,run_id=run,patient_id=pairs$patient_id[i],sample_id=pairs$sample_id[i],time_raw=pairs$time_raw[i],filtered_F=filtF_count,merged=merged_count,merge_of_filtered_F_pct=round(ifelse(filtF_count>0,100*merged_count/filtF_count,NA_real_),2))
    rm(merger); gc(verbose=FALSE)
  }
  stats <- bind_rows(stats)
  write_excel_csv(stats,file.path(OUT,'PRJNA578267_samplewise_merge_tracking.csv'),na='')

  logmsg('FINAL SEQUENCE TABLE BUILD START')
  mergers <- lapply(pairs$run_id,function(run) readRDS(file.path(CM,paste0(run,'.rds')))); names(mergers) <- pairs$run_id
  seqtab <- makeSequenceTable(mergers)
  if(!is.matrix(seqtab)||nrow(seqtab)==0||ncol(seqtab)==0) stop('Invalid/empty pre-chimera seqtab.')
  nochim <- removeBimeraDenovo(seqtab,method='consensus',multithread=FALSE,verbose=TRUE)
  if(!is.matrix(nochim)||nrow(nochim)==0||ncol(nochim)==0) stop('Invalid/empty nonchim seqtab.')

  pre_path <- file.path(PWORK,paste0(PROJECT,'_seqtab_prechim_step84B3_samplewise.rds'))
  no_path <- file.path(PWORK,paste0(PROJECT,'_seqtab_nochim_step84B3_samplewise.rds'))
  atomic_saveRDS(seqtab,pre_path); atomic_saveRDS(nochim,no_path)

  nc <- rowSums(nochim); aligned <- setNames(rep(0,207),pairs$run_id); aligned[names(nc)] <- nc
  final_track <- stats |> mutate(nonchim=as.numeric(aligned[run_id]),nonchim_of_merged_pct=round(ifelse(merged>0,100*nonchim/merged,NA_real_),2))
  write_excel_csv(final_track,file.path(OUT,'PRJNA578267_FINAL_samplewise_tracking.csv'),na='')

  lens <- nchar(colnames(nochim))
  asvqc <- tibble(asv_sequence=colnames(nochim),asv_length=lens,total_abundance=colSums(nochim),within_expected_length=lens>=350 & lens<=500)
  write_excel_csv(asvqc,file.path(OUT,'PRJNA578267_FINAL_ASV_length_QC.csv'),na='')
  total_m <- sum(final_track$merged,na.rm=TRUE); total_n <- sum(final_track$nonchim,na.rm=TRUE)
  summary <- tibble(project=PROJECT,runs=207,samplewise_checkpointing=TRUE,zero_merged_samples=sum(final_track$merged==0),zero_nonchim_samples=sum(final_track$nonchim==0),median_merge_of_filtered_F_pct=round(median(final_track$merge_of_filtered_F_pct,na.rm=TRUE),2),total_merged_reads=total_m,total_nonchim_reads=total_n,nonchim_of_merged_pct=round(ifelse(total_m>0,100*total_n/total_m,NA_real_),2),ASVs_prechim=ncol(seqtab),ASVs_nochim=ncol(nochim),median_ASV_length=median(lens),ASVs_350_500bp_pct=round(100*mean(asvqc$within_expected_length),2),status='SAMPLEWISE_FULL_DADA2_COMPLETE_PENDING_FINAL_REVIEW')
  write_excel_csv(summary,file.path(OUT,'PRJNA578267_FINAL_samplewise_summary.csv'),na='')
  writeLines(c(paste0('Completed: ',Sys.time()),'Project: PRJNA578267','Mode: sample-wise checkpointed DADA2',paste0('Final nochim table: ',no_path),'Status: COMPLETE_PENDING_FINAL_REVIEW'),file.path(PWORK,'_STEP84B3_SAMPLEWISE_COMPLETE.ok'))
  logmsg('STEP84B3 COMPLETE')
  print(summary,n=Inf,width=Inf)
}

tryCatch(main(),error=function(e){ msg <- c(paste0('STEP84B3 FATAL ERROR: ',Sys.time()),paste0('Message: ',conditionMessage(e)),paste0('Call: ',paste(deparse(conditionCall(e)),collapse=' '))); writeLines(msg,ERR); message(paste(msg,collapse='\n')); quit(save='no',status=1,runLast=FALSE) })
