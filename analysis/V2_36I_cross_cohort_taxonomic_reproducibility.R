
# ============================================================
# V2_36I CROSS-COHORT TAXONOMIC REPRODUCIBILITY
#
# Central question:
# Are ecosystem-level longitudinal displacement directions more
# reproducible across cohorts than individual taxonomic effects?
#
# CRITICAL DESIGN RULE:
#   DO NOT merge ASV/OTU abundance matrices across cohorts.
#   This analysis uses already-computed WITHIN-COHORT paired CLR
#   taxonomic effects from frozen Step89A2 outputs, then harmonizes
#   taxon labels/effect directions only at the effect level.
#
# Primary cohorts:
#   PRJNA691455
#   PRJNA851469
#   PRJNA516701
#
# Secondary/context:
#   PRJEB82425   supportive external ICU infection cohort
#   PRJNA578267  non-sepsis longitudinal control
#
# R 4.4.x
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

set.seed(20260825)

ROOT <- "E:/sepsis_project"
RESULTS <- file.path(ROOT, "results")

SOURCE_DIR <- file.path(
  RESULTS,
  "V2_29A2_STEP89A_TAXONOMIC_PAIRED_TRAJECTORIES_FIXED"
)

OUT <- file.path(
  RESULTS,
  "V2_36I_CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY"
)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

primary_projects <- c(
  "PRJNA691455",
  "PRJNA851469",
  "PRJNA516701"
)

secondary_projects <- c(
  "PRJEB82425",
  "PRJNA578267"
)

all_projects <- c(primary_projects, secondary_projects)
ranks <- c("GENUS","FAMILY")

# Conservative taxon eligibility for cross-cohort comparison.
MIN_PAIRS <- 8
MIN_PREVALENCE_EITHER <- 0.20
N_PERM <- 5000

# ------------------------------------------------------------
# 1. Locate exact frozen Step89A paired-CLR result files
# ------------------------------------------------------------
rank_folder <- c(
  GENUS="02_GENUS",
  FAMILY="03_FAMILY"
)

file_for <- function(project, rank) {
  suffix <- ifelse(rank=="GENUS",
                   "_genus_paired_CLR_results.csv",
                   "_family_paired_CLR_results.csv")
  file.path(
    SOURCE_DIR,
    rank_folder[[rank]],
    paste0(project, suffix)
  )
}

registry <- expand_grid(
  project=all_projects,
  rank=ranks
) %>%
  mutate(
    file=map2_chr(project,rank,file_for),
    exists=file.exists(file)
  )

write_csv(
  registry,
  file.path(OUT,"00_INPUT_REGISTRY.csv")
)

missing_primary <- registry %>%
  filter(project %in% primary_projects, !exists)

if(nrow(missing_primary)>0) {
  stop(
    "Missing primary Step89A paired CLR file(s): ",
    paste(missing_primary$file, collapse=" ; ")
  )
}

# ------------------------------------------------------------
# 2. Canonicalize taxon labels
# ------------------------------------------------------------
canonicalize_taxon <- function(x, rank) {
  x <- as.character(x)
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_trim(x)

  if(rank=="GENUS") {
    m <- str_match(x, "Genus__([^|;]+)")
    out <- ifelse(!is.na(m[,2]), m[,2], x)
  } else {
    m <- str_match(x, "Family__([^|;]+)")
    out <- ifelse(!is.na(m[,2]), m[,2], x)
  }

  out <- str_trim(out)
  out <- str_replace_all(out, "^\\[|\\]$", "")
  out
}

is_named_taxon <- function(x) {
  lx <- tolower(str_trim(as.character(x)))
  !is.na(lx) &
    lx!="" &
    !str_detect(
      lx,
      "unclassified|uncultured|unknown|unassigned|^na$|^none$|metagenome"
    )
}

# ------------------------------------------------------------
# 3. Read and validate each within-cohort result table
# ------------------------------------------------------------
required_cols <- c(
  "project","tax_rank","subset_type","required_anchor",
  "early_label","late_label","taxon","n_pairs",
  "prevalence_early","prevalence_late",
  "mean_clr_difference","paired_effect_dz",
  "wilcoxon_p","direction","fdr"
)

read_effect_file <- function(project,rank) {

  f <- file_for(project,rank)

  if(!file.exists(f)) return(tibble())

  x <- read_csv(f, show_col_types=FALSE, progress=FALSE)

  miss <- setdiff(required_cols,names(x))
  if(length(miss)>0) {
    stop(
      project," ",rank," missing columns: ",
      paste(miss,collapse=", ")
    )
  }

  # A Step89A result file should normally represent one paired comparison.
  comparison_registry <- x %>%
    distinct(
      project, tax_rank, subset_type, required_anchor,
      early_label, late_label
    )

  # If several comparisons somehow exist, prefer one explicitly marked primary.
  if(nrow(comparison_registry)>1) {

    primary_rows <- x %>%
      filter(
        str_detect(
          tolower(paste(subset_type,required_anchor)),
          "primary"
        )
      )

    if(nrow(primary_rows)>0) {
      x <- primary_rows
    } else {
      stop(
        project," ",rank,
        " contains >1 paired-comparison definition and none is explicitly PRIMARY. ",
        "Inspect 01_COMPARISON_DEFINITION_REGISTRY.csv before proceeding."
      )
    }
  }

  x %>%
    mutate(
      canonical_taxon=canonicalize_taxon(taxon,rank),
      named_taxon=is_named_taxon(canonical_taxon),
      effect=as.numeric(mean_clr_difference),
      prevalence_max=pmax(
        as.numeric(prevalence_early),
        as.numeric(prevalence_late),
        na.rm=TRUE
      ),
      eligible=
        is.finite(effect) &
        n_pairs>=MIN_PAIRS &
        prevalence_max>=MIN_PREVALENCE_EITHER
    )
}

effects_list <- list()
comparison_defs <- list()

for(prj in all_projects) {
  for(rk in ranks) {

    f <- file_for(prj,rk)
    if(!file.exists(f)) next

    x <- read_effect_file(prj,rk)

    effects_list[[paste(prj,rk,sep="__")]] <- x

    comparison_defs[[paste(prj,rk,sep="__")]] <-
      x %>%
      distinct(
        project,tax_rank,subset_type,required_anchor,
        early_label,late_label,n_pairs
      )
  }
}

comparison_registry <- bind_rows(comparison_defs)
write_csv(
  comparison_registry,
  file.path(OUT,"01_COMPARISON_DEFINITION_REGISTRY.csv")
)

all_effects <- bind_rows(effects_list)

write_csv(
  all_effects,
  file.path(OUT,"02_HARMONIZED_WITHIN_COHORT_TAXON_EFFECTS.csv")
)

# ------------------------------------------------------------
# 4. Pairwise taxonomic concordance
# ------------------------------------------------------------
pairwise_one <- function(project_a,project_b,rank,named_only=TRUE) {

  a <- all_effects %>%
    filter(
      project==project_a,
      tax_rank==rank,
      eligible
    )

  b <- all_effects %>%
    filter(
      project==project_b,
      tax_rank==rank,
      eligible
    )

  if(named_only) {
    a <- a %>% filter(named_taxon)
    b <- b %>% filter(named_taxon)
  }

  a <- a %>%
    group_by(canonical_taxon) %>%
    summarise(
      effect_a=mean(effect,na.rm=TRUE),
      dz_a=mean(paired_effect_dz,na.rm=TRUE),
      .groups="drop"
    )

  b <- b %>%
    group_by(canonical_taxon) %>%
    summarise(
      effect_b=mean(effect,na.rm=TRUE),
      dz_b=mean(paired_effect_dz,na.rm=TRUE),
      .groups="drop"
    )

  m <- inner_join(a,b,by="canonical_taxon") %>%
    filter(is.finite(effect_a),is.finite(effect_b))

  n <- nrow(m)

  if(n<5) {
    return(tibble(
      project_a=project_a,
      project_b=project_b,
      rank=rank,
      named_only=named_only,
      n_shared=n,
      spearman_rho=NA_real_,
      spearman_p=NA_real_,
      sign_agreement=NA_real_,
      top_k=NA_integer_,
      top_jaccard=NA_real_,
      perm_p_rho=NA_real_,
      perm_p_sign=NA_real_,
      perm_p_top_jaccard=NA_real_
    ))
  }

  ct <- suppressWarnings(
    cor.test(
      m$effect_a,m$effect_b,
      method="spearman",
      exact=FALSE
    )
  )

  rho_obs <- unname(ct$estimate)
  rho_p <- ct$p.value

  nonzero <- sign(m$effect_a)!=0 & sign(m$effect_b)!=0
  sign_obs <- if(sum(nonzero)>0) {
    mean(sign(m$effect_a[nonzero])==sign(m$effect_b[nonzero]))
  } else NA_real_

  k <- max(5, ceiling(.20*n))
  k <- min(k,20,n)

  top_a <- m %>%
    arrange(desc(abs(effect_a))) %>%
    slice_head(n=k) %>%
    pull(canonical_taxon)

  top_b <- m %>%
    arrange(desc(abs(effect_b))) %>%
    slice_head(n=k) %>%
    pull(canonical_taxon)

  jac_obs <- length(intersect(top_a,top_b)) /
    length(union(top_a,top_b))

  perm_rho <- rep(NA_real_,N_PERM)
  perm_sign <- rep(NA_real_,N_PERM)
  perm_jac <- rep(NA_real_,N_PERM)

  for(i in seq_len(N_PERM)) {

    eb <- sample(m$effect_b,replace=FALSE)

    perm_rho[i] <- suppressWarnings(
      cor(m$effect_a,eb,method="spearman")
    )

    nz <- sign(m$effect_a)!=0 & sign(eb)!=0
    if(sum(nz)>0) {
      perm_sign[i] <- mean(
        sign(m$effect_a[nz])==sign(eb[nz])
      )
    }

    tmp <- tibble(
      taxon=m$canonical_taxon,
      ea=m$effect_a,
      eb=eb
    )

    ta <- tmp %>%
      arrange(desc(abs(ea))) %>%
      slice_head(n=k) %>%
      pull(taxon)

    tb <- tmp %>%
      arrange(desc(abs(eb))) %>%
      slice_head(n=k) %>%
      pull(taxon)

    perm_jac[i] <- length(intersect(ta,tb)) /
      length(union(ta,tb))
  }

  perm_p_rho <- (
    1+sum(abs(perm_rho)>=abs(rho_obs),na.rm=TRUE)
  ) / (
    1+sum(is.finite(perm_rho))
  )

  perm_p_sign <- (
    1+sum(perm_sign>=sign_obs,na.rm=TRUE)
  ) / (
    1+sum(is.finite(perm_sign))
  )

  perm_p_jac <- (
    1+sum(perm_jac>=jac_obs,na.rm=TRUE)
  ) / (
    1+sum(is.finite(perm_jac))
  )

  tibble(
    project_a=project_a,
    project_b=project_b,
    rank=rank,
    named_only=named_only,
    n_shared=n,
    spearman_rho=rho_obs,
    spearman_p=rho_p,
    sign_agreement=sign_obs,
    top_k=k,
    top_jaccard=jac_obs,
    perm_p_rho=perm_p_rho,
    perm_p_sign=perm_p_sign,
    perm_p_top_jaccard=perm_p_jac
  )
}

primary_pairs <- combn(primary_projects,2,simplify=FALSE)

primary_pairwise <- bind_rows(
  lapply(
    ranks,
    function(rk) bind_rows(
      lapply(
        primary_pairs,
        function(pp) pairwise_one(pp[1],pp[2],rk,TRUE)
      )
    )
  )
)

write_csv(
  primary_pairwise,
  file.path(OUT,"03_PRIMARY3_PAIRWISE_TAXONOMIC_CONCORDANCE.csv")
)

# Sensitivity including unclassified labels where canonically matchable.
primary_pairwise_alltaxa <- bind_rows(
  lapply(
    ranks,
    function(rk) bind_rows(
      lapply(
        primary_pairs,
        function(pp) pairwise_one(pp[1],pp[2],rk,FALSE)
      )
    )
  )
)

write_csv(
  primary_pairwise_alltaxa,
  file.path(OUT,"04_PRIMARY3_PAIRWISE_SENSITIVITY_ALL_TAXA.csv")
)

# ------------------------------------------------------------
# 5. Three-way shared taxa across all three primary cohorts
# ------------------------------------------------------------
threeway_one <- function(rank,named_only=TRUE) {

  xs <- lapply(primary_projects,function(prj) {

    x <- all_effects %>%
      filter(
        project==prj,
        tax_rank==rank,
        eligible
      )

    if(named_only) x <- x %>% filter(named_taxon)

    x %>%
      group_by(canonical_taxon) %>%
      summarise(
        effect=mean(effect,na.rm=TRUE),
        .groups="drop"
      ) %>%
      rename(!!prj := effect)
  })

  m <- reduce(xs,inner_join,by="canonical_taxon")

  if(nrow(m)==0) return(tibble())

  mat <- as.matrix(
    m[,primary_projects,drop=FALSE]
  )

  sign_mat <- sign(mat)

  m$all_same_direction <- apply(
    sign_mat,1,function(z) {
      z <- z[z!=0]
      length(z)==3 && length(unique(z))==1
    }
  )

  m$n_positive <- rowSums(sign_mat>0)
  m$n_negative <- rowSums(sign_mat<0)
  m$rank <- rank
  m$named_only <- named_only

  m
}

threeway <- bind_rows(
  lapply(
    ranks,
    function(rk) threeway_one(rk,TRUE)
  )
)

write_csv(
  threeway,
  file.path(OUT,"05_PRIMARY3_SHARED_TAXA_DIRECTION_PATTERNS.csv")
)

threeway_summary <- threeway %>%
  group_by(rank) %>%
  summarise(
    n_taxa_shared_all3=n(),
    n_same_direction_all3=sum(all_same_direction),
    proportion_same_direction_all3=mean(all_same_direction),
    n_all_positive=sum(n_positive==3),
    n_all_negative=sum(n_negative==3),
    .groups="drop"
  )

write_csv(
  threeway_summary,
  file.path(OUT,"06_PRIMARY3_THREEWAY_DIRECTION_SUMMARY.csv")
)

# ------------------------------------------------------------
# 6. Secondary/contextual comparisons
# ------------------------------------------------------------
context_pairs <- list(
  c("PRJNA691455","PRJEB82425"),
  c("PRJNA851469","PRJEB82425"),
  c("PRJNA516701","PRJEB82425"),
  c("PRJNA691455","PRJNA578267"),
  c("PRJNA851469","PRJNA578267"),
  c("PRJNA516701","PRJNA578267")
)

context_results <- bind_rows(
  lapply(
    ranks,
    function(rk) bind_rows(
      lapply(
        context_pairs,
        function(pp) {
          if(
            all(
              registry$exists[
                registry$project %in% pp &
                registry$rank==rk
              ]
            )
          ) {
            pairwise_one(pp[1],pp[2],rk,TRUE)
          } else {
            tibble()
          }
        }
      )
    )
  )
)

write_csv(
  context_results,
  file.path(OUT,"07_CONTEXT_EXTERNAL_AND_CONTROL_CONCORDANCE.csv")
)

# ------------------------------------------------------------
# 7. Compare ecological-level reproducibility with taxonomic level
# ------------------------------------------------------------
# Frozen ecological evidence:
# robust progressive displacement in all 3 primary cohorts.
ecological <- tibble(
  project=primary_projects,
  ecological_direction=c(
    "PROGRESSIVE",
    "PROGRESSIVE",
    "PROGRESSIVE"
  )
)

ecological_direction_concordance <- 1.0

rank_summary <- primary_pairwise %>%
  group_by(rank) %>%
  summarise(
    n_pairwise_comparisons=n(),
    median_spearman_rho=median(spearman_rho,na.rm=TRUE),
    min_spearman_rho=min(spearman_rho,na.rm=TRUE),
    max_spearman_rho=max(spearman_rho,na.rm=TRUE),
    median_sign_agreement=median(sign_agreement,na.rm=TRUE),
    median_top_jaccard=median(top_jaccard,na.rm=TRUE),
    pairwise_rho_perm_significant=sum(perm_p_rho<0.05,na.rm=TRUE),
    pairwise_sign_perm_significant=sum(perm_p_sign<0.05,na.rm=TRUE),
    .groups="drop"
  )

synthesis <- rank_summary %>%
  mutate(
    ecological_primary_direction_concordance=
      ecological_direction_concordance,
    interpretation_class=case_when(
      abs(median_spearman_rho)<0.20 ~ "LOW_TAXONOMIC_EFFECT_CONCORDANCE",
      abs(median_spearman_rho)<0.50 ~ "MODEST_TAXONOMIC_EFFECT_CONCORDANCE",
      TRUE ~ "HIGH_TAXONOMIC_EFFECT_CONCORDANCE"
    )
  )

write_csv(
  synthesis,
  file.path(OUT,"08_ECOLOGICAL_VS_TAXONOMIC_SYNTHESIS.csv")
)

# ------------------------------------------------------------
# 8. Evidence decision
# Transparent, descriptive thresholds only.
# ------------------------------------------------------------
genus_class <- synthesis %>%
  filter(rank=="GENUS") %>%
  pull(interpretation_class)

family_class <- synthesis %>%
  filter(rank=="FAMILY") %>%
  pull(interpretation_class)

if(length(genus_class)==0) genus_class <- NA_character_
if(length(family_class)==0) family_class <- NA_character_

both_not_high <- all(
  c(genus_class,family_class) !=
    "HIGH_TAXONOMIC_EFFECT_CONCORDANCE",
  na.rm=TRUE
)

tier <- if(
  ecological_direction_concordance==1 &&
  both_not_high
) {
  "ECOLOGICAL_DIRECTION_MORE_REPRODUCIBLE_THAN_TAXONOMIC_EFFECTS"
} else {
  "TAXONOMIC_CONCORDANCE_NOT_CLEarly_LOWER_THAN_ECOLOGICAL_DIRECTION"
}

decision <- tibble(
  primary_ecological_direction_concordance=
    ecological_direction_concordance,
  genus_class=genus_class[1],
  family_class=family_class[1],
  evidence_tier=tier
)

write_csv(
  decision,
  file.path(OUT,"09_EVIDENCE_DECISION.csv")
)

# ------------------------------------------------------------
# 9. Candidate figure
# ------------------------------------------------------------
plot_dat <- primary_pairwise %>%
  mutate(
    pair=paste(project_a,project_b,sep=" vs "),
    rank=factor(rank,levels=c("GENUS","FAMILY"))
  )

p1 <- ggplot(
  plot_dat,
  aes(x=pair,y=spearman_rho,shape=rank)
) +
  geom_hline(yintercept=0,linetype=2) +
  geom_point(size=3) +
  coord_flip() +
  labs(
    x=NULL,
    y="Spearman correlation of within-cohort CLR effect vectors",
    title="Cross-cohort taxonomic effect concordance"
  ) +
  theme_classic(base_size=10)

ggsave(
  file.path(OUT,"10_Figure_candidate_pairwise_effect_concordance.pdf"),
  p1,width=7.0,height=4.8
)
ggsave(
  file.path(OUT,"10_Figure_candidate_pairwise_effect_concordance.png"),
  p1,width=7.0,height=4.8,dpi=300
)

# ------------------------------------------------------------
# 10. Manuscript-safe summary
# ------------------------------------------------------------
fmt <- function(x,d=3) {
  ifelse(is.na(x),"NA",formatC(x,format="f",digits=d))
}

lines <- c(
  "V2 STEP36I CROSS-COHORT TAXONOMIC REPRODUCIBILITY",
  "",
  "Primary comparison: PRJNA691455, PRJNA851469, PRJNA516701.",
  "All taxonomic effects were estimated within cohort in frozen Step89A paired CLR analyses.",
  "No ASV/OTU abundance matrices were merged across cohorts.",
  "",
  paste0(
    "Ecological displacement direction concordance across the three primary cohorts: ",
    fmt(ecological_direction_concordance,2),
    " (3/3 progressive)."
  )
)

for(i in seq_len(nrow(synthesis))) {
  z <- synthesis[i,]
  lines <- c(
    lines,
    paste0(
      z$rank,
      ": median pairwise taxonomic effect-vector rho=",
      fmt(z$median_spearman_rho),
      "; median sign agreement=",
      fmt(z$median_sign_agreement),
      "; median top-effect Jaccard=",
      fmt(z$median_top_jaccard),
      "; class=",
      z$interpretation_class,
      "."
    )
  )
}

lines <- c(
  lines,
  "",
  paste0("EVIDENCE TIER: ",tier),
  "",
  "Interpretation guardrails:",
  "1. This analysis compares reproducibility at different biological resolutions; ecological and taxonomic metrics are not numerically interchangeable.",
  "2. Do not claim that all patients or cohorts share identical ecological trajectories.",
  "3. Do not claim a universal sepsis-specific taxonomic signature.",
  "4. Low taxonomic concordance should be framed as heterogeneous taxonomic routes accompanying a more reproducible ecological displacement pattern.",
  "5. The three primary cohorts are the inferential focus; PRJEB82425 and PRJNA578267 are contextual sensitivity comparisons.",
  "6. No previously frozen V2 result was modified."
)

writeLines(
  lines,
  file.path(OUT,"11_MANUSCRIPT_SAFE_SUMMARY.txt")
)

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP36I COMPLETE",
    paste0("Evidence tier: ",tier),
    "No previously frozen V2 inference modified."
  ),
  file.path(OUT,"_STEP36I_COMPLETE.txt")
)

cat("STEP36I COMPLETE\n")
cat("Evidence tier:",tier,"\n")
