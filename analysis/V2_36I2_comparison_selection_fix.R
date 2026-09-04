
suppressPackageStartupMessages({
library(readr)
library(dplyr)
library(stringr)
})

ROOT <- "E:/sepsis_project"
SOURCE <- file.path(ROOT,"results",
"V2_29A2_STEP89A_TAXONOMIC_PAIRED_TRAJECTORIES_FIXED")

OUT <- file.path(ROOT,"results",
"V2_36I2_CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY_FIXED")

dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

projects <- c("PRJNA691455","PRJNA851469","PRJNA516701")
ranks <- c("GENUS","FAMILY")

getfile <- function(p,r){
file.path(
SOURCE,
ifelse(r=="GENUS","02_GENUS","03_FAMILY"),
paste0(p,ifelse(r=="GENUS",
"_genus_paired_CLR_results.csv",
"_family_paired_CLR_results.csv"))
)
}

choose <- function(x){

d <- x %>%
distinct(
subset_type,
required_anchor,
early_label,
late_label
) %>%
mutate(
explicit_primary=str_detect(
tolower(paste(subset_type,required_anchor)),
"primary"
),
longitudinal=str_detect(
tolower(paste(required_anchor,early_label,late_label)),
"day1|baseline|admission"
)
)

d$selected <- FALSE
reason <- "manual_review_needed"

i <- which(d$explicit_primary)

if(length(i)==1){
d$selected[i] <- TRUE
reason <- "explicit_primary"
}else{
i <- which(d$longitudinal)
if(length(i)>=1){
d$selected[i[1]] <- TRUE
reason <- "main_longitudinal_match"
}
}

d$selection_reason <- reason
d
}

out <- list()

for(p in projects){
for(r in ranks){

f <- getfile(p,r)

if(file.exists(f)){

x <- read_csv(f,show_col_types=FALSE)

a <- choose(x) %>%
mutate(project=p,rank=r,source_file=f)

out[[paste(p,r)]] <- a
}
}
}

audit <- bind_rows(out)

write_csv(
audit,
file.path(OUT,"01_COMPARISON_SELECTION_AUDIT.csv")
)

write_csv(
audit %>% filter(selected),
file.path(OUT,"02_SELECTED_PRIMARY_COMPARISONS.csv")
)

writeLines(
c(
paste0("Completed: ",Sys.time()),
"STEP36I2 COMPLETE"
),
file.path(OUT,"_STEP36I2_COMPLETE.txt")
)

cat("STEP36I2 COMPLETE\n")
