
options(stringsAsFactors = FALSE)

# Sepsis V2 - Step 65 metadata coverage audit
# Designed for R 4.4.0 on Windows

pkgs <- c("readr","readxl","dplyr","stringr","purrr","tibble","xml2")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(xml2)
})

DATA_ROOT <- "E:/sepsis_project/data"
OUT_ROOT  <- "E:/sepsis_project/results/V2_05_metadata_coverage"
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

PROJECTS <- c(
  "PRJEB37289","PRJEB67798","PRJEB82425",
  "PRJNA516701","PRJNA578267","PRJNA595346","PRJNA884103"
)

patterns <- list(
  patient_id = c("\\bpatient\\b","\\bsubject\\b","participant","patient.?id","subject.?id"),
  sample_id = c("sample.?id","sample_accession","biosample","sequence.?name","sample.?name","specimen"),
  run_id = c("run_accession","sra.?run","\\brun.?id\\b","\\bsrr[0-9]+","\\berr[0-9]+","\\bcrr[0-9]+"),
  time = c("time.?point","timepoint","\\bday\\b","\\bdol\\b","collection.?date","sampling.?date","visit","baseline","discharge","intubation","icu.?day","hospital.?day"),
  infection_group = c("\\bsepsis\\b","septic.?shock","\\bvap\\b","pneumonia","infected","bacteremia","bacteraemia","bloodstream.?infection","\\bbsi\\b"),
  infection_source = c("infection.?source","infection.?site","site.?of.?infection","focus.?of.?infection","pulmonary.?infection","urinary.?infection","abdominal.?origin","bloodstream.?infection"),
  outcome = c("\\boutcome\\b","mortality","\\bdeath\\b","survival","survivor","deceased","icu.?length.?of.?stay","hospital.?length.?of.?stay"),
  antibiotics = c("antibiotic","antimicrobial","\\babx\\b","cefazolin","clindamycin","tobramycin","colistin","cefotaxime","antibiotic.?score"),
  severity = c("\\bsofa\\b","\\bapache\\b","severity","shock","vasopressor","mechanical.?ventilation","ventilation.?duration"),
  demographics = c("\\bage\\b","\\bsex\\b","gender","birth.?weight","gestational.?age","\\bbmi\\b"),
  body_site = c("body.?site","sample.?type","rectal","stool","faecal","fecal","tracheal","endotracheal","gut","lung")
)

detect_domains <- function(text) {
  text <- tolower(paste(text, collapse = " "))
  vapply(patterns, function(pats) {
    any(vapply(pats, function(p) str_detect(text, regex(p, ignore_case = TRUE)), logical(1)))
  }, logical(1))
}

read_tabular <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    x <- tryCatch(read_csv(path, n_max = 3000, show_col_types = FALSE, progress = FALSE), error=function(e) NULL)
    return(if (is.null(x)) list() else list(list(sheet="", data=x)))
  }
  if (ext %in% c("tsv","txt")) {
    x <- tryCatch(read_tsv(path, n_max = 3000, show_col_types = FALSE, progress = FALSE), error=function(e) NULL)
    if (!is.null(x) && ncol(x) > 1) return(list(list(sheet="", data=x)))
    x <- tryCatch(read_delim(path, delim = NULL, n_max = 3000, show_col_types = FALSE, progress = FALSE), error=function(e) NULL)
    return(if (is.null(x)) list() else list(list(sheet="", data=x)))
  }
  if (ext %in% c("xlsx","xls")) {
    sh <- tryCatch(excel_sheets(path), error=function(e) character())
    out <- list()
    for (s in sh) {
      x <- tryCatch(read_excel(path, sheet=s, n_max=3000), error=function(e) NULL)
      if (!is.null(x)) out[[length(out)+1]] <- list(sheet=s, data=x)
    }
    return(out)
  }
  list()
}

docx_text <- function(path) {
  td <- tempfile("docx_")
  dir.create(td)
  on.exit(unlink(td, recursive=TRUE, force=TRUE), add=TRUE)
  ok <- tryCatch({ unzip(path, exdir=td); TRUE }, error=function(e) FALSE)
  if (!ok) return("")
  xmls <- list.files(file.path(td,"word"), pattern="\\.xml$", full.names=TRUE, recursive=TRUE)
  paste(map_chr(xmls, function(x) {
    tryCatch(xml_text(read_xml(x)), error=function(e) "")
  }), collapse=" ")
}

text_file <- function(path) {
  tryCatch(paste(readLines(path, warn=FALSE, encoding="UTF-8"), collapse=" "),
           error=function(e) "")
}

files_for_project <- function(project) {
  base <- file.path(DATA_ROOT, project, "00_metadata")
  if (!dir.exists(base)) return(character())
  x <- list.files(base, recursive=TRUE, full.names=TRUE)
  x[file.exists(x) & !dir.exists(x)]
}

domains <- names(patterns)
file_rows <- list()
detected_src <- setNames(lapply(PROJECTS, function(x) setNames(replicate(length(domains), character(), simplify=FALSE), domains)), PROJECTS)
structured_src <- detected_src
unparsed <- setNames(replicate(length(PROJECTS), character(), simplify=FALSE), PROJECTS)

for (project in PROJECTS) {
  fs <- files_for_project(project)
  message("Scanning ", project, ": ", length(fs), " files")

  for (path in fs) {
    ext <- tolower(tools::file_ext(path))
    rel <- sub(paste0("^", gsub("\\\\","/", normalizePath(file.path(DATA_ROOT, project), winslash="/", mustWork=FALSE)), "/?"), "",
               normalizePath(path, winslash="/", mustWork=FALSE))

    det <- setNames(rep(FALSE, length(domains)), domains)
    struc <- det
    status <- "UNPARSED"
    evidence <- ""

    if (ext %in% c("csv","tsv","txt","xlsx","xls")) {
      tabs <- read_tabular(path)
      if (length(tabs) > 0) status <- "TABULAR_PARSED" else status <- "TABULAR_READ_FAILED"

      ev <- character()
      for (item in tabs) {
        df <- item$data
        if (is.null(df) || ncol(df)==0) next
        hd <- detect_domains(names(df))
        vals <- character()
        for (j in seq_len(min(ncol(df),100))) {
          x <- unique(na.omit(as.character(df[[j]])))
          vals <- c(vals, head(x[x!=""],20))
        }
        vd <- detect_domains(vals)
        for (d in domains) {
          det[d] <- det[d] || hd[d] || vd[d]
          struc[d] <- struc[d] || hd[d]
        }
        ev <- c(ev, paste0(ifelse(item$sheet=="","",paste0("sheet=",item$sheet,"; ")),
                           "columns=", paste(head(names(df),40), collapse=";")))
      }
      evidence <- paste(ev, collapse=" || ")
    }

    if (ext %in% c("md","rmd","r","py","json","xml")) {
      status <- "TEXT_PARSED"
      x <- detect_domains(text_file(path))
      det <- det | x
    }

    if (ext == "docx") {
      status <- "DOCX_TEXT_PARSED"
      x <- detect_domains(docx_text(path))
      det <- det | x
    }

    if (ext %in% c("pdf","rdata","rds","biom","zip","gz","tgz","tar")) {
      status <- "PRESENT_NOT_PARSED"
      unparsed[[project]] <- c(unparsed[[project]], rel)
    }

    fn <- detect_domains(basename(path))
    det <- det | fn

    for (d in domains) {
      if (det[d]) detected_src[[project]][[d]] <- c(detected_src[[project]][[d]], rel)
      if (struc[d]) structured_src[[project]][[d]] <- c(structured_src[[project]][[d]], rel)
    }

    file_rows[[length(file_rows)+1]] <- tibble(
      Project=project,
      RelativePath=rel,
      FileName=basename(path),
      Extension=ext,
      SizeKB=round(file.info(path)$size/1024,2),
      ParseStatus=status,
      DetectedDomains=paste(names(det)[det], collapse=";"),
      StructuredDomains=paste(names(struc)[struc], collapse=";"),
      EvidencePreview=evidence
    )
  }
}

file_inventory <- if (length(file_rows)>0) bind_rows(file_rows) else tibble()

coverage_rows <- list()
manual_rows <- list()
priority_rows <- list()

CORE <- c("patient_id","sample_id","time")
CLINICAL <- c("infection_group","infection_source","outcome","antibiotics","severity")

for (project in PROJECTS) {
  row <- list(Project=project)
  n_struct <- 0
  n_detect <- 0

  for (d in domains) {
    detected <- length(unique(detected_src[[project]][[d]])) > 0
    structured <- length(unique(structured_src[[project]][[d]])) > 0

    st <- if (structured) "STRUCTURED" else if (detected) "PRESENT_UNSTRUCTURED" else "NOT_FOUND"
    if (structured) n_struct <- n_struct + 1
    if (detected) n_detect <- n_detect + 1

    row[[paste0(d,"_status")]] <- st
    row[[paste0(d,"_sources")]] <- paste(head(unique(detected_src[[project]][[d]]),10), collapse=" | ")
  }

  core_status <- vapply(CORE, function(d) row[[paste0(d,"_status")]], character(1))
  ready <- if (all(core_status=="STRUCTURED")) "YES" else if (all(core_status!="NOT_FOUND")) "MANUAL_LINKAGE_REVIEW" else "NO"

  row$CoreMappingReady <- ready
  row$StructuredDomainCount <- n_struct
  row$DetectedDomainCount <- n_detect
  row$UnparsedPotentialFiles <- paste(head(unique(unparsed[[project]]),20), collapse=" | ")

  coverage_rows[[length(coverage_rows)+1]] <- as_tibble(row)

  missing_core <- CORE[vapply(CORE, function(d) row[[paste0(d,"_status")]]=="NOT_FOUND", logical(1))]
  missing_clin <- CLINICAL[vapply(CLINICAL, function(d) row[[paste0(d,"_status")]]=="NOT_FOUND", logical(1))]
  unstruct <- c(CORE,CLINICAL)[vapply(c(CORE,CLINICAL), function(d) row[[paste0(d,"_status")]]=="PRESENT_UNSTRUCTURED", logical(1))]

  if (length(missing_core)>0) {
    pri <- "HIGH"
    act <- "Core patient/sample/time linkage is still missing. Inspect existing supplement/repository first; contact authors only if public data truly lack it."
  } else if (length(unstruct)>0) {
    pri <- "MEDIUM"
    act <- "No new data collection yet. Existing files appear to contain the information but require parsing/linkage."
  } else if (length(missing_clin)>0) {
    pri <- "LOW"
    act <- "Core longitudinal mapping is available. Collect missing clinical covariates only if required by the planned model."
  } else {
    pri <- "NONE"
    act <- "No manual collection currently required."
  }

  if (project=="PRJEB67798" && row$patient_id_status!="STRUCTURED") {
    pri <- "HIGH"
    act <- "Parse existing supplement first. If patient-level linkage is still absent, author contact is likely required because the publication states public source data exclude patient data."
  }

  manual_rows[[length(manual_rows)+1]] <- tibble(
    Project=project,
    Priority=pri,
    MissingCore=paste(missing_core,collapse=";"),
    MissingClinical=paste(missing_clin,collapse=";"),
    PresentButUnstructured=paste(unstruct,collapse=";"),
    RecommendedAction=act
  )

  score <- ifelse(ready=="YES",100,0) + n_struct*10 + n_detect
  next_step <- if (ready=="YES") "BUILD_SAMPLE_PATIENT_TIME_MAP" else if (length(missing_core)==0) "PARSE_EXISTING_RESOURCES" else "FILL_CORE_MAPPING_GAPS"

  priority_rows[[length(priority_rows)+1]] <- tibble(
    Project=project,
    MappingReadiness=ready,
    CoverageScore=score,
    RecommendedNextStep=next_step
  )
}

coverage <- bind_rows(coverage_rows)
manual <- bind_rows(manual_rows)
priority <- bind_rows(priority_rows) |> arrange(desc(CoverageScore), Project)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

p1 <- file.path(OUT_ROOT, paste0("V2_metadata_coverage_audit_",stamp,".csv"))
p2 <- file.path(OUT_ROOT, paste0("V2_metadata_file_inventory_",stamp,".csv"))
p3 <- file.path(OUT_ROOT, paste0("V2_manual_collection_required_",stamp,".csv"))
p4 <- file.path(OUT_ROOT, paste0("V2_mapping_priority_",stamp,".csv"))

write_excel_csv(coverage,p1,na="")
write_excel_csv(file_inventory,p2,na="")
write_excel_csv(manual,p3,na="")
write_excel_csv(priority,p4,na="")

cat("\n================ STEP 65 COMPLETE ================\n")
cat("R version:", R.version.string, "\n\n")
print(manual |> select(Project,Priority,MissingCore,MissingClinical,PresentButUnstructured), n=Inf, width=Inf)
cat("\nOutputs:\n",p1,"\n",p2,"\n",p3,"\n",p4,"\n")
cat("==================================================\n")
