# ============================================================
# Sepsis V2 - Step 81
# FINAL master merge + longitudinal eligibility audit + cohort freeze
# ============================================================

options(stringsAsFactors = FALSE)

pkgs <- c("readr","dplyr","tidyr","stringr","purrr","tibble")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing) > 0) install.packages(missing, repos="https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(tibble)
})

ROOT <- "E:/sepsis_project"
STEP73_ROOT <- file.path(ROOT,"results","V2_13_full_master_metadata")
STEP77_ROOT <- file.path(ROOT,"results","V2_17_finalize_68229_adjudicate_1125274")
STEP80F_ROOT <- file.path(ROOT,"results","V2_20F_exact_final_mapping")
OUT <- file.path(ROOT,"results","V2_21_FINAL_COHORT_FREEZE")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

LOG <- file.path(OUT,"_STEP81_runtime_checkpoints.txt")
ERR <- file.path(OUT,"_STEP81_FATAL_ERROR.txt")
if (file.exists(ERR)) unlink(ERR)
ck <- function(x) cat(paste0(x,": ",Sys.time(),"\n"),file=LOG,append=TRUE)

clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)|x==""|tolower(x)%in%c("na","nan","n/a","null","none")] <- NA_character_
  x
}
clean_num <- function(x) suppressWarnings(as.numeric(clean_chr(x)))

safe_csv <- function(path) {
  if (length(path)==0 || is.na(path) || !file.exists(path)) return(NULL)
  suppressMessages(read_csv(path,show_col_types=FALSE,progress=FALSE,name_repair="unique"))
}

find_latest <- function(folder,pattern) {
  if (!dir.exists(folder)) return(NA_character_)
  x <- list.files(folder,pattern=pattern,full.names=TRUE)
  if (!length(x)) return(NA_character_)
  x[which.max(file.info(x)$mtime)]
}

get_chr <- function(df,candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (!length(hit)) return(rep(NA_character_,nrow(df)))
  clean_chr(df[[hit[1]]])
}
get_num <- function(df,candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (!length(hit)) return(rep(NA_real_,nrow(df)))
  clean_num(df[[hit[1]]])
}
combine_severity <- function(sofa,apache) {
  sofa <- clean_chr(sofa); apache <- clean_chr(apache)
  out <- rep(NA_character_,length(sofa))
  for (i in seq_along(out)) {
    z <- character()
    if (!is.na(sofa[i])) z <- c(z,paste0("SOFA=",sofa[i]))
    if (!is.na(apache[i])) z <- c(z,paste0("APACHEII=",apache[i]))
    if (length(z)) out[i] <- paste(z,collapse=";")
  }
  out
}
count_sequence_files <- function(project) {
  roots <- c(file.path(ROOT,"data",project))
  if (project=="PRJNA533528") roots <- c(roots,"F:/sepsis/PRJNA533528")
  files <- character()
  for (r in roots) if (dir.exists(r)) files <- c(files,list.files(r,recursive=TRUE,full.names=TRUE))
  files <- unique(files)
  sum(str_detect(basename(files),regex("\\.(fastq|fq|fa|fasta)\\.gz$",ignore_case=TRUE)))
}

main <- function() {
  ck("STEP81 STARTED")

  step73_path <- find_latest(STEP73_ROOT,"^V2_FULL_master_metadata_[0-9]{8}_[0-9]{6}\\.csv$")
  p682_path <- file.path(STEP77_ROOT,"PRJEB68229_FINAL_patient_time_map.csv")
  p112_path <- file.path(STEP77_ROOT,"PRJNA1125274_STRICT_PROVISIONAL_patient_time_map.csv")
  cra_path <- file.path(STEP80F_ROOT,"CRA002354_MASTER_READY_FINAL.csv")
  p851_path <- file.path(STEP80F_ROOT,"PRJNA851469_MASTER_READY_FINAL.csv")

  inputs <- tibble(
    source_stage=c("STEP73_BASE","STEP77_PRJEB68229_FINAL","STEP77_PRJNA1125274_PROVISIONAL","STEP80F_CRA002354_FINAL","STEP80F_PRJNA851469_FINAL"),
    path=c(step73_path,p682_path,p112_path,cra_path,p851_path)
  ) |>
    mutate(
      exists=!is.na(path)&file.exists(path),
      bytes=map_dbl(path,~if(!is.na(.x)&&file.exists(.x)) file.info(.x)$size else NA_real_),
      md5=map_chr(path,~if(!is.na(.x)&&file.exists(.x)) unname(tools::md5sum(.x)) else NA_character_)
    )

  write_excel_csv(inputs,file.path(OUT,"V2_STEP81_input_provenance.csv"),na="")
  if (!all(inputs$exists)) stop("One or more Step81 inputs are missing; see V2_STEP81_input_provenance.csv")

  base <- safe_csv(step73_path); p682 <- safe_csv(p682_path); p112 <- safe_csv(p112_path)
  cra <- safe_csv(cra_path); p851 <- safe_csv(p851_path)
  if (any(vapply(list(base,p682,p112,cra,p851),is.null,logical(1)))) stop("Failed to read one or more Step81 inputs.")
  ck("INPUTS READ")

  if (nrow(base)<1900) stop(paste0("Step73 master unexpectedly small: ",nrow(base)))
  if (nrow(cra)!=131 || sum(!is.na(clean_chr(cra$run_id)))!=131 || n_distinct(clean_chr(cra$run_id),na.rm=TRUE)!=131)
    stop("CRA002354 Step80F guard failed: expected 131 rows / 131 exact unique Runs.")
  if (nrow(p851)!=123 || sum(!is.na(clean_chr(p851$run_id)))!=119)
    stop("PRJNA851469 Step80F guard failed: expected 123 rows / 119 exact Runs.")
  if (nrow(p682)!=190 || n_distinct(get_chr(p682,c("Patient_ID","patient_id")),na.rm=TRUE)!=96)
    stop("PRJEB68229 guard failed: expected 190 samples / 96 patients.")
  if (nrow(p112)<280 || n_distinct(get_chr(p112,c("Patient_ID_Provisional","patient_candidate","patient_id")),na.rm=TRUE)<120)
    stop("PRJNA1125274 provisional map guard failed.")
  ck("INPUT VERSION GUARDS PASSED")

  replace_projects <- c("PRJEB68229","PRJNA1125274","CRA002354","PRJNA851469")
  base <- base |>
    filter(!(project %in% replace_projects)) |>
    select(-any_of(c("patient_uid","sample_uid","run_uid","n_distinct_timepoints","longitudinal_ge2","longitudinal_ge3",
                     "sequence_n_distinct_timepoints","sequence_longitudinal_ge2","sequence_longitudinal_ge3","analysis_module",
                     "freeze_metadata_include","primary_16s_sequence_include","shotgun_external_include","freeze_status",
                     "freeze_source_stage","freeze_version","time_key"))) |>
    mutate(
      baseline_sepsis_status=NA_character_, sample_sepsis_status=NA_character_,
      infection_source_standard=NA_character_, outcome_28d=NA_character_, sample_class=NA_character_,
      mapping_status=NA_character_, sequence_analysis_include=NA
    )

  new682 <- tibble(
    project="PRJEB68229",
    patient_id=get_chr(p682,c("Patient_ID","patient_id")),
    sample_id=get_chr(p682,c("Sample_ID","sample_id","Sample_Alias")),
    biosample=get_chr(p682,c("BioSample","biosample")),
    run_id=get_chr(p682,c("Run_ID","run_id")), experiment_id=NA_character_,
    time_raw=get_chr(p682,c("Time_Raw","time_raw")), time_day=get_num(p682,c("Time_Day","time_day")),
    time_class=get_chr(p682,c("Time_Standard","time_class")), time_order=get_num(p682,c("Time_Order","time_order")),
    time_anchor=get_chr(p682,c("Time_Anchor","time_anchor")), collection_date=get_chr(p682,c("Collection_Date","collection_date")),
    body_site="fecal/stool", phenotype="critical_illness_ICU", control_type=NA_character_, infection_group=NA_character_,
    outcome=NA_character_, antibiotics=NA_character_, severity=NA_character_, intervention=NA_character_, intervention_arm=NA_character_,
    organ_dysfunction_type=NA_character_, organ_dysfunction_status=NA_character_, sex=NA_character_, age=NA_character_,
    gestational_age=NA_character_, birthweight=NA_character_,
    primary_sample_flag=ifelse(get_chr(p682,c("Time_Raw"))=="D1","YES","NO"),
    analysis_role="ICU_BACKGROUND_LONGITUDINAL", dataset_split="EXTERNAL_VALIDATION",
    cohort_role="ICU_BACKGROUND_LONGITUDINAL", analysis_tier="V2_LONGITUDINAL_VALIDATION",
    sequence_status="PUBLIC_RUN_MAPPING_READY_SEQUENCE_PROCESSING_PENDING", data_modality="16S",
    metadata_source="Step77_PRJEB68229_FINAL_patient_time_map",
    label_source="ENA submission structure + public supplement", label_confidence="HIGH", review_status="FINALIZED_STEP77",
    qc_status=ifelse(get_chr(p682,c("Collection_Date_QC"))=="SOURCE_DATE_ORDER_ANOMALY","PASS_MAPPING_DATE_ANOMALY_RETAINED","PASS_FINAL_MAPPING"),
    notes=get_chr(p682,c("Outcome_Note")), legacy_module=NA_character_,
    baseline_sepsis_status=NA_character_, sample_sepsis_status=NA_character_, infection_source_standard=NA_character_, outcome_28d=NA_character_,
    sample_class="ICU_PATIENT", mapping_status=get_chr(p682,c("Mapping_Status","Mapping_Confidence")),
    sequence_analysis_include=!is.na(get_chr(p682,c("Run_ID")))
  )

  cra_sofa <- get_chr(cra,c("baseline_sofa","SOFA (Sequential organ failure assessment)","SOFA"))
  cra_apache <- get_chr(cra,c("baseline_apache_ii","APACHE II (Acute physiology and chronic health evaluation)","APACHE II"))

  newCRA <- tibble(
    project="CRA002354", patient_id=get_chr(cra,c("patient_id")), sample_id=get_chr(cra,c("sample_id")),
    biosample=get_chr(cra,c("archive_sample_accession","biosample")), run_id=get_chr(cra,c("run_id")),
    experiment_id=get_chr(cra,c("experiment_id")), time_raw=get_chr(cra,c("time_raw")), time_day=get_num(cra,c("time_day")),
    time_class=get_chr(cra,c("time_raw")), time_order=get_num(cra,c("time_order")),
    time_anchor="ICU longitudinal fecal collection day", collection_date=NA_character_, body_site="fecal/stool",
    phenotype=get_chr(cra,c("baseline_sepsis_status","sepsis_status")), control_type=NA_character_,
    infection_group=get_chr(cra,c("infection_source","infection_group")), outcome=get_chr(cra,c("outcome_28d","outcome")),
    antibiotics=get_chr(cra,c("antibiotics_during_icu","Antibiotics use (during ICU stay)","Antibiotics use")),
    severity=combine_severity(cra_sofa,cra_apache), intervention=NA_character_, intervention_arm=NA_character_,
    organ_dysfunction_type=NA_character_, organ_dysfunction_status=NA_character_, sex=get_chr(cra,c("Gender","sex")),
    age=get_chr(cra,c("Age (year)","age")), gestational_age=NA_character_, birthweight=NA_character_,
    primary_sample_flag=ifelse(get_chr(cra,c("sample_id"))==get_chr(cra,c("first_sample_id")),"YES","NO"),
    analysis_role="CORE_SEPSIS_LONGITUDINAL", dataset_split="V2_CORE_DISCOVERY_AND_CLINICAL_TRAJECTORY",
    cohort_role="CORE_SEPSIS_LONGITUDINAL_WITH_CLINICAL_OUTCOMES", analysis_tier="V2_CORE_LONGITUDINAL",
    sequence_status="GSA_RUN_MAPPING_EXACT_131_OF_131_SEQUENCE_PROCESSING_PENDING", data_modality="16S",
    metadata_source="Step80F_CRA002354_MASTER_READY_FINAL", label_source="Supplementary Tables S1/S2 + GSA technical workbook",
    label_confidence="HIGH", review_status="FINALIZED_STEP80F", qc_status="PASS_FINAL_MAPPING_131_OF_131",
    notes="Patient-level baseline sepsis/shock and 28-day outcome retained; sample-level sepsis/shock trajectory retained.",
    legacy_module=NA_character_, baseline_sepsis_status=get_chr(cra,c("baseline_sepsis_status")),
    sample_sepsis_status=get_chr(cra,c("sample_sepsis_status")), infection_source_standard=get_chr(cra,c("infection_source")),
    outcome_28d=get_chr(cra,c("outcome_28d")), sample_class="SEPSIS_PATIENT",
    mapping_status=get_chr(cra,c("run_mapping_status","patient_mapping_status")), sequence_analysis_include=TRUE
  )

  p851_class <- get_chr(p851,c("sample_class")); p851_run <- get_chr(p851,c("run_id"))
  new851 <- tibble(
    project="PRJNA851469", patient_id=get_chr(p851,c("patient_id")), sample_id=get_chr(p851,c("sample_id")),
    biosample=NA_character_, run_id=p851_run, experiment_id=NA_character_, time_raw=get_chr(p851,c("time_raw")),
    time_day=get_num(p851,c("time_day")), time_class=get_chr(p851,c("time_raw")), time_order=get_num(p851,c("time_order")),
    time_anchor=ifelse(p851_class=="ICU_PATIENT","ICU admission / Day 3 / Day 7",NA_character_), collection_date=NA_character_,
    body_site="rectal_swab", phenotype=ifelse(p851_class=="HEALTHY_CONTROL","healthy","critical_illness_ICU"),
    control_type=ifelse(p851_class=="HEALTHY_CONTROL","healthy_control",NA_character_), infection_group=NA_character_, outcome=NA_character_,
    antibiotics=NA_character_, severity=NA_character_, intervention=NA_character_, intervention_arm=NA_character_,
    organ_dysfunction_type=NA_character_, organ_dysfunction_status=NA_character_, sex=NA_character_, age=NA_character_,
    gestational_age=NA_character_, birthweight=NA_character_,
    primary_sample_flag=ifelse(get_num(p851,c("time_day"))==1,"YES","NO"),
    analysis_role=ifelse(p851_class=="ICU_PATIENT","ICU_BACKGROUND_LONGITUDINAL","HEALTHY_STATIC_CONTROL"),
    dataset_split="EXTERNAL_VALIDATION",
    cohort_role=ifelse(p851_class=="ICU_PATIENT","ICU_BACKGROUND_LONGITUDINAL","HEALTHY_STATIC_CONTROL_WITHIN_PRJNA851469"),
    analysis_tier=ifelse(p851_class=="ICU_PATIENT","V2_LONGITUDINAL_VALIDATION","V2_STATIC_CONTROL"),
    sequence_status=ifelse(!is.na(p851_run),"STRICT_EXACT_RUN_MAPPING_SEQUENCE_PROCESSING_PENDING","SEQUENCE_UNMAPPED_EXCLUDED_FROM_SEQUENCE_LAYER"),
    data_modality="16S", metadata_source="Step80F_PRJNA851469_MASTER_READY_FINAL",
    label_source="Nature Medicine Supplementary Table 2 + exact ENA mapping",
    label_confidence=ifelse(!is.na(p851_run),"HIGH","METADATA_HIGH_RUN_UNMAPPED"), review_status="FINALIZED_STEP80F",
    qc_status=ifelse(!is.na(p851_run),"PASS_STRICT_EXACT_MAPPING","PASS_METADATA_RUN_UNMAPPED"),
    notes=ifelse(!is.na(p851_run),NA_character_,"Retained in metadata master but excluded from sequence layer because no strict Run mapping."),
    legacy_module=NA_character_, baseline_sepsis_status=NA_character_, sample_sepsis_status=NA_character_,
    infection_source_standard=NA_character_, outcome_28d=NA_character_, sample_class=p851_class,
    mapping_status=get_chr(p851,c("run_mapping_status")), sequence_analysis_include=!is.na(p851_run)
  )

  p112_patient <- get_chr(p112,c("Patient_ID_Provisional","patient_candidate","patient_id"))
  p112_sample <- get_chr(p112,c("sample_alias","Sample_Alias","sample_id","sample_accession"))
  p112_run <- get_chr(p112,c("run_id","Run_ID","run_accession"))
  new112 <- tibble(
    project="PRJNA1125274", patient_id=p112_patient, sample_id=p112_sample,
    biosample=get_chr(p112,c("sample_accession","BioSample")), run_id=p112_run, experiment_id=NA_character_,
    time_raw=get_chr(p112,c("Time_Raw","alias_letter")), time_day=NA_real_, time_class=get_chr(p112,c("Time_Standard_Candidate")),
    time_order=get_num(p112,c("Time_Order_Candidate")), time_anchor="PROVISIONAL study-stage mapping",
    collection_date=get_chr(p112,c("collection_date")), body_site="rectal_swab", phenotype="sepsis_longitudinal_candidate",
    control_type=NA_character_, infection_group=NA_character_, outcome=NA_character_, antibiotics=NA_character_, severity=NA_character_,
    intervention=NA_character_, intervention_arm=NA_character_, organ_dysfunction_type=NA_character_, organ_dysfunction_status=NA_character_,
    sex=NA_character_, age=NA_character_, gestational_age=NA_character_, birthweight=NA_character_,
    primary_sample_flag=ifelse(get_chr(p112,c("Time_Raw","alias_letter"))=="A","YES","NO"),
    analysis_role="DEFERRED_SEPSIS_LONGITUDINAL_CANDIDATE", dataset_split="METADATA_ONLY_DEFERRED",
    cohort_role="CORE_SEPSIS_LONGITUDINAL_CANDIDATE_AFTER_RAW_REPAIR", analysis_tier="V2_DEFERRED_PROVISIONAL",
    sequence_status="CORRECT_289_RUN_SET_NOT_READY_FOR_PRIMARY_ANALYSIS", data_modality="16S",
    metadata_source="Step77_PRJNA1125274_STRICT_PROVISIONAL_MAP",
    label_source="High-confidence study-design inference; not explicit ENA patient/time fields",
    label_confidence="PROVISIONAL_HIGH_CONFIDENCE_INFERENCE", review_status="DEFERRED_NOT_FROZEN_FOR_SEQUENCE",
    qc_status="PASS_PROVISIONAL_METADATA_ONLY",
    notes="Included in metadata registry only. Do not run primary DADA2/trajectory analysis until the correct 289-Run raw set and mapping are finalized.",
    legacy_module=NA_character_, baseline_sepsis_status=NA_character_, sample_sepsis_status=NA_character_,
    infection_source_standard=NA_character_, outcome_28d=NA_character_, sample_class="SEPSIS_PATIENT_PROVISIONAL",
    mapping_status=get_chr(p112,c("Patient_Mapping_Status","Time_Mapping_Status")), sequence_analysis_include=FALSE
  )

  ck("NEW COHORTS STANDARDIZED")

  full <- bind_rows(base,new682,newCRA,new851,new112) |>
    mutate(
      project=clean_chr(project), patient_id=clean_chr(patient_id), sample_id=clean_chr(sample_id),
      biosample=clean_chr(biosample), run_id=clean_chr(run_id), experiment_id=clean_chr(experiment_id),
      time_raw=clean_chr(time_raw), time_class=clean_chr(time_class), cohort_role=clean_chr(cohort_role),
      analysis_tier=clean_chr(analysis_tier), sequence_status=clean_chr(sequence_status), data_modality=clean_chr(data_modality),
      freeze_source_stage=case_when(
        project=="PRJEB68229" ~ "STEP77_FINAL",
        project=="PRJNA1125274" ~ "STEP77_PROVISIONAL",
        project=="CRA002354" ~ "STEP80F_FINAL",
        project=="PRJNA851469" ~ "STEP80F_FINAL",
        TRUE ~ "STEP73_BASE"
      )
    )

  full <- full |>
    mutate(
      analysis_module=case_when(
        project %in% c("PRJEB33360","PRJNA691455","CRA002354") ~ "CORE_SEPSIS_LONGITUDINAL",
        project=="PRJEB82425" ~ "ICU_INFECTION_LONGITUDINAL_EXTERNAL",
        project %in% c("PRJNA516701","PRJNA595346","PRJEB68229","PRJNA851469") ~ "ICU_BACKGROUND_LONGITUDINAL",
        project %in% c("PRJNA578267","PRJEB67798") ~ "NONSEPSIS_LONGITUDINAL_CONTROL",
        project %in% c("PRJNA430161","PRJNA1166732") ~ "INTERVENTION_LONGITUDINAL_SUPPORT",
        project=="PRJNA912621" ~ "ORGAN_DYSFUNCTION_LONGITUDINAL_SUPPORT",
        project %in% c("PRJNA797231","PRJNA978257","PRJNA1010969") ~ "STATIC_SUPPORT",
        project=="PRJNA884103" ~ "SHOTGUN_BSI_LONGITUDINAL_EXTERNAL",
        project=="PRJNA1125274" ~ "DEFERRED_SEPSIS_LONGITUDINAL_CANDIDATE",
        TRUE ~ "OTHER_REVIEWED_SUPPORT"
      ),
      patient_uid=ifelse(is.na(patient_id),NA_character_,paste(project,patient_id,sep="::")),
      sample_uid=ifelse(is.na(sample_id),NA_character_,paste(project,sample_id,sep="::")),
      run_uid=ifelse(is.na(run_id),NA_character_,paste(project,run_id,sep="::"))
    )

  run_conflicts <- full |>
    filter(!is.na(run_uid)) |>
    group_by(run_uid) |>
    summarise(n_rows=n(),patients=n_distinct(patient_id,na.rm=TRUE),samples=n_distinct(sample_id,na.rm=TRUE),
              timepoints=n_distinct(coalesce(time_class,time_raw),na.rm=TRUE),.groups="drop") |>
    filter(patients>1|samples>1|timepoints>1)
  dup_run <- full |> filter(!is.na(run_uid)) |> count(run_uid,name="n") |> filter(n>1)
  write_excel_csv(run_conflicts,file.path(OUT,"V2_FROZEN_run_conflict_QC.csv"),na="")
  write_excel_csv(dup_run,file.path(OUT,"V2_FROZEN_duplicate_run_QC.csv"),na="")
  if (nrow(run_conflicts)>0 || nrow(dup_run)>0)
    stop(paste0("Cannot freeze: conflicts=",nrow(run_conflicts)," duplicates=",nrow(dup_run)))
  ck("RUN DUPLICATE AND CONFLICT QC PASSED")

  full <- full |> mutate(time_key=coalesce(clean_chr(time_class),clean_chr(time_raw)))
  meta_pt <- full |> filter(!is.na(patient_uid),!is.na(time_key)) |> distinct(patient_uid,time_key) |> count(patient_uid,name="n_distinct_timepoints")
  full <- full |> left_join(meta_pt,by="patient_uid") |>
    mutate(n_distinct_timepoints=coalesce(n_distinct_timepoints,0L),
           longitudinal_ge2=n_distinct_timepoints>=2,
           longitudinal_ge3=n_distinct_timepoints>=3)

  full <- full |>
    mutate(
      freeze_metadata_include=TRUE,
      primary_16s_sequence_include=case_when(
        data_modality!="16S" ~ FALSE,
        project=="PRJNA1125274" ~ FALSE,
        is.na(run_id) ~ FALSE,
        sequence_analysis_include %in% FALSE ~ FALSE,
        TRUE ~ TRUE
      ),
      shotgun_external_include=data_modality=="SHOTGUN",
      freeze_status=case_when(
        project=="PRJNA1125274" ~ "FROZEN_METADATA_ONLY_PROVISIONAL_SEQUENCE_DEFERRED",
        data_modality=="SHOTGUN" ~ "FROZEN_SHOTGUN_EXTERNAL_SEPARATE_LAYER",
        data_modality=="16S" & primary_16s_sequence_include ~ "FROZEN_16S_SEQUENCE_LAYER",
        data_modality=="16S" & !primary_16s_sequence_include ~ "FROZEN_METADATA_SEQUENCE_UNAVAILABLE_OR_UNMAPPED",
        TRUE ~ "FROZEN_METADATA"
      )
    )

  seq_pt <- full |> filter(primary_16s_sequence_include,!is.na(patient_uid),!is.na(time_key)) |>
    distinct(patient_uid,time_key) |> count(patient_uid,name="sequence_n_distinct_timepoints")
  full <- full |> left_join(seq_pt,by="patient_uid") |>
    mutate(sequence_n_distinct_timepoints=coalesce(sequence_n_distinct_timepoints,0L),
           sequence_longitudinal_ge2=sequence_n_distinct_timepoints>=2,
           sequence_longitudinal_ge3=sequence_n_distinct_timepoints>=3)

  patient_summary <- full |> filter(!is.na(patient_uid)) |>
    group_by(project,patient_uid,patient_id,analysis_module) |>
    summarise(
      n_rows=n(),n_samples=n_distinct(sample_uid,na.rm=TRUE),n_runs=n_distinct(run_uid,na.rm=TRUE),
      n_timepoints=n_distinct(time_key[!is.na(time_key)]),
      n_sequence_timepoints=n_distinct(time_key[primary_16s_sequence_include & !is.na(time_key)]),
      metadata_longitudinal_ge2=n_timepoints>=2,metadata_longitudinal_ge3=n_timepoints>=3,
      sequence_longitudinal_ge2=n_sequence_timepoints>=2,sequence_longitudinal_ge3=n_sequence_timepoints>=3,
      .groups="drop"
    )

  cohort_summary <- full |> group_by(project,analysis_module,data_modality) |>
    summarise(
      rows=n(),patients=n_distinct(patient_uid,na.rm=TRUE),samples=n_distinct(sample_uid,na.rm=TRUE),
      mapped_runs=n_distinct(run_uid,na.rm=TRUE),missing_run_rows=sum(is.na(run_id)),
      metadata_patients_ge2=n_distinct(patient_uid[longitudinal_ge2],na.rm=TRUE),
      metadata_patients_ge3=n_distinct(patient_uid[longitudinal_ge3],na.rm=TRUE),
      sequence_patients_ge2=n_distinct(patient_uid[sequence_longitudinal_ge2],na.rm=TRUE),
      sequence_patients_ge3=n_distinct(patient_uid[sequence_longitudinal_ge3],na.rm=TRUE),
      primary_16s_rows=sum(primary_16s_sequence_include),outcome_nonmissing=sum(!is.na(outcome)),
      infection_group_nonmissing=sum(!is.na(infection_group)),.groups="drop"
    )

  module_summary <- full |> group_by(analysis_module) |>
    summarise(
      projects=n_distinct(project),patients=n_distinct(patient_uid,na.rm=TRUE),samples=n_distinct(sample_uid,na.rm=TRUE),
      mapped_runs=n_distinct(run_uid,na.rm=TRUE),metadata_patients_ge2=n_distinct(patient_uid[longitudinal_ge2],na.rm=TRUE),
      sequence_patients_ge2=n_distinct(patient_uid[sequence_longitudinal_ge2],na.rm=TRUE),
      primary_16s_rows=sum(primary_16s_sequence_include),.groups="drop"
    )

  registry <- cohort_summary |>
    mutate(
      freeze_decision=case_when(
        project=="PRJNA1125274" ~ "METADATA_ONLY_DEFERRED",
        data_modality=="SHOTGUN" ~ "SHOTGUN_EXTERNAL_SEPARATE_LAYER",
        analysis_module=="STATIC_SUPPORT" ~ "STATIC_SUPPORT",
        TRUE ~ "INCLUDE"
      ),
      next_sequence_action=case_when(
        project=="PRJNA1125274" ~ "Finish correct 289-Run raw download/mapping before DADA2",
        project=="PRJNA884103" ~ "Keep shotgun analysis separate; raw processing deferred",
        project=="CRA002354" ~ "Process cohort-specific 16S reads / ASV table",
        project=="PRJNA851469" ~ "Use only 119 strict exact-mapped samples in sequence layer",
        TRUE ~ "Use existing ASV or run cohort-specific sequence processing as indicated by sequence_status"
      )
    )

  manifest16s <- full |> filter(primary_16s_sequence_include) |>
    select(project,patient_id,patient_uid,sample_id,sample_uid,biosample,experiment_id,run_id,run_uid,
           time_raw,time_class,time_day,time_order,analysis_module,cohort_role,sequence_status,freeze_status) |>
    arrange(project,patient_id,time_order,time_day,sample_id)

  manifest16s_longitudinal <- full |> filter(primary_16s_sequence_include,sequence_longitudinal_ge2) |>
    select(project,patient_id,patient_uid,sample_id,sample_uid,run_id,run_uid,time_raw,time_class,time_day,time_order,
           analysis_module,cohort_role,sequence_n_distinct_timepoints) |>
    arrange(project,patient_id,time_order,time_day,sample_id)

  deferred <- full |> filter(project=="PRJNA1125274" | !primary_16s_sequence_include | data_modality=="SHOTGUN") |>
    select(project,patient_id,sample_id,run_id,time_raw,analysis_module,data_modality,sequence_status,freeze_status,notes)

  projects <- sort(unique(full$project[!is.na(full$project)]))
  seq_readiness <- tibble(project=projects,local_sequence_files=map_int(projects,count_sequence_files)) |>
    left_join(cohort_summary |> select(project,data_modality,mapped_runs,missing_run_rows,primary_16s_rows),by="project") |>
    mutate(local_sequence_file_status=case_when(
      local_sequence_files>0 ~ "LOCAL_SEQUENCE_FILES_FOUND",
      primary_16s_rows>0 ~ "NO_RAW_FILES_DETECTED_CHECK_EXISTING_ASV_OR_DOWNLOAD",
      data_modality=="SHOTGUN" ~ "SHOTGUN_SEPARATE_LAYER",
      TRUE ~ "NO_SEQUENCE_ACTION_IN_PRIMARY_LAYER"
    ))

  global_summary <- tibble(
    metric=c("Rows","Projects","Patients","Samples","Mapped_Runs","Metadata_Patients_GE2","Metadata_Patients_GE3",
             "Primary_16S_Sequence_Rows","Primary_16S_Patients_GE2","Primary_16S_Patients_GE3","Shotgun_Rows",
             "Deferred_PRJNA1125274_Rows","Missing_Run_Rows","Duplicate_Run_UID","Run_Conflict_Rows"),
    value=c(
      nrow(full),n_distinct(full$project),n_distinct(full$patient_uid,na.rm=TRUE),n_distinct(full$sample_uid,na.rm=TRUE),
      n_distinct(full$run_uid,na.rm=TRUE),n_distinct(full$patient_uid[full$longitudinal_ge2],na.rm=TRUE),
      n_distinct(full$patient_uid[full$longitudinal_ge3],na.rm=TRUE),sum(full$primary_16s_sequence_include),
      n_distinct(full$patient_uid[full$sequence_longitudinal_ge2],na.rm=TRUE),
      n_distinct(full$patient_uid[full$sequence_longitudinal_ge3],na.rm=TRUE),
      sum(full$data_modality=="SHOTGUN",na.rm=TRUE),sum(full$project=="PRJNA1125274"),sum(is.na(full$run_id)),
      nrow(dup_run),nrow(run_conflicts)
    )
  )

  stamp <- format(Sys.time(),"%Y%m%d_%H%M%S")
  freeze_version <- paste0("STEP81_",stamp)
  full <- full |> mutate(freeze_version=freeze_version)

  full_path <- file.path(OUT,paste0("V2_FROZEN_master_metadata_",stamp,".csv"))
  patient_path <- file.path(OUT,paste0("V2_FROZEN_patient_summary_",stamp,".csv"))
  cohort_path <- file.path(OUT,paste0("V2_FROZEN_cohort_summary_",stamp,".csv"))
  global_path <- file.path(OUT,paste0("V2_FROZEN_global_summary_",stamp,".csv"))

  write_excel_csv(full,full_path,na="")
  write_excel_csv(patient_summary,patient_path,na="")
  write_excel_csv(cohort_summary,cohort_path,na="")
  write_excel_csv(global_summary,global_path,na="")
  write_excel_csv(module_summary,file.path(OUT,"V2_FROZEN_analysis_module_summary.csv"),na="")
  write_excel_csv(registry,file.path(OUT,"V2_FROZEN_project_registry.csv"),na="")
  write_excel_csv(manifest16s,file.path(OUT,"V2_FROZEN_16S_sequence_manifest.csv"),na="")
  write_excel_csv(manifest16s_longitudinal,file.path(OUT,"V2_FROZEN_16S_longitudinal_manifest.csv"),na="")
  write_excel_csv(deferred,file.path(OUT,"V2_FROZEN_deferred_registry.csv"),na="")
  write_excel_csv(seq_readiness,file.path(OUT,"V2_FROZEN_sequence_readiness_by_project.csv"),na="")

  readme <- c(
    "SEPSIS V2 - STEP81 FINAL COHORT FREEZE",
    paste0("Freeze version: ",freeze_version),
    paste0("Created: ",Sys.time()),"",
    "UNIQUE METADATA ENTRY POINT",
    paste0("Use this file for all downstream analyses: ",full_path),"",
    "FREEZE RULES",
    "- PRJNA851469: only 119/123 strict exact Run mappings enter the sequence layer; 4 unmapped rows remain metadata-only.",
    "- CRA002354: 131/131 exact Sample -> Experiment -> CRR mapping.",
    "- PRJEB68229: D1/D3 are preserved raw labels; use ordinal S1/S2 time order, not literal day values.",
    "- PRJNA1125274: metadata retained as provisional/deferred; excluded from primary 16S sequence analysis until the correct 289-Run raw set is finalized.",
    "- Shotgun projects are kept separate from the 16S main layer.",
    "- Different 16S cohorts must be processed separately; do not merge raw ASV feature IDs across different amplicon regions.","",
    "NEXT STEP",
    "Step82 should use V2_FROZEN_16S_sequence_manifest.csv and V2_FROZEN_sequence_readiness_by_project.csv to decide which cohorts already have ASV tables and which require DADA2/raw-read processing.",
    "Do not return to candidate/review-only metadata files unless a documented conflict is discovered."
  )
  writeLines(readme,file.path(OUT,"README_STEP81_FINAL_FREEZE.txt"),useBytes=TRUE)

  ck("STEP81 COMPLETE")

  cat("\n============================================================\n")
  cat("SEPSIS V2 - STEP81 FINAL COHORT FREEZE COMPLETE\n")
  cat("============================================================\n\n")
  cat("GLOBAL SUMMARY:\n"); print(global_summary,n=Inf,width=Inf)
  cat("\nPROJECT REGISTRY:\n"); print(registry,n=Inf,width=Inf)
  cat("\nANALYSIS MODULE SUMMARY:\n"); print(module_summary,n=Inf,width=Inf)
  cat("\nFrozen master:\n",full_path,"\n",sep="")
  cat("\nOutput folder:\n",OUT,"\n",sep="")
  cat("============================================================\n")
}

tryCatch(
  main(),
  error=function(e) {
    msg <- c(
      paste0("STEP81 FATAL ERROR: ",Sys.time()),
      paste0("Message: ",conditionMessage(e)),
      paste0("Call: ",paste(deparse(conditionCall(e)),collapse=" "))
    )
    writeLines(msg,ERR,useBytes=TRUE)
    ck("STEP81 FAILED")
    message(paste(msg,collapse="\n"))
    quit(save="no",status=1,runLast=FALSE)
  }
)
