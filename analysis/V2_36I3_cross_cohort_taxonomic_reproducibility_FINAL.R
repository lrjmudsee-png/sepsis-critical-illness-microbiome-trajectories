
# ============================================================
# V2_36I3 CROSS-COHORT TAXONOMIC REPRODUCIBILITY - FINAL FIX
#
# Primary rule resolved from Step36I2 audit:
#   ALL_PAIRED                = PRIMARY
#   COMMON_ANCHOR_SENSITIVITY = SENSITIVITY
#
# This matches the frozen V2 longitudinal evidence framework.
#
# Important:
# - Effects are estimated WITHIN cohort only (from Step89A2).
# - No ASV/OTU abundance matrices are pooled across cohorts.
# - Taxon labels/effect directions are harmonized only AFTER
#   within-cohort estimation.
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
RESULTS <- file.path(ROOT,"results")

SOURCE_DIR <- file.path(
  RESULTS,
  "V2_29A2_STEP89A_TAXONOMIC_PAIRED_TRAJECTORIES_FIXED"
)

OUT <- file.path(
  RESULTS,
  "V2_36I3_CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY_FINAL"
)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

primary_projects <- c(
  "PRJNA691455",
  "PRJNA851469",
  "PRJNA516701"
)

context_projects <- c(
  "PRJEB82425",
  "PRJNA578267"
)

all_projects <- c(primary_projects,context_projects)

ranks <- c("GENUS","FAMILY")

subset_definitions <- c(
  PRIMARY="ALL_PAIRED",
  SENSITIVITY="COMMON_ANCHOR_SENSITIVITY"
)

MIN_PAIRS <- 8
MIN_PREVALENCE <- 0.20
N_PERM <- 5000

rank_dir <- c(
  GENUS="02_GENUS",
  FAMILY="03_FAMILY"
)

get_file <- function(project,rank){
  suffix <- ifelse(
    rank=="GENUS",
    "_genus_paired_CLR_results.csv",
    "_family_paired_CLR_results.csv"
  )

  file.path(
    SOURCE_DIR,
    rank_dir[[rank]],
    paste0(project,suffix)
  )
}

canonicalize_taxon <- function(x,rank){

  x <- str_trim(as.character(x))

  if(rank=="GENUS"){
    mm <- str_match(x,"Genus__([^|;]+)")
    out <- ifelse(!is.na(mm[,2]),mm[,2],x)
  }else{
    mm <- str_match(x,"Family__([^|;]+)")
    out <- ifelse(!is.na(mm[,2]),mm[,2],x)
  }

  out <- str_trim(out)
  out <- str_replace_all(out,"^\\[|\\]$","")
  out
}

is_named <- function(x){

  lx <- tolower(str_trim(as.character(x)))

  !is.na(lx) &
    lx != "" &
    !str_detect(
      lx,
      "unclassified|uncultured|unknown|unassigned|metagenome|^na$|^none$"
    )
}

# ------------------------------------------------------------
# 1. Input registry + exact comparison audit
# ------------------------------------------------------------
registry <- expand_grid(
  project=all_projects,
  rank=ranks
) %>%
  mutate(
    file=map2_chr(project,rank,get_file),
    exists=file.exists(file)
  )

write_csv(
  registry,
  file.path(OUT,"00_INPUT_REGISTRY.csv")
)

missing_primary <- registry %>%
  filter(
    project %in% primary_projects,
    !exists
  )

if(nrow(missing_primary)>0){
  stop(
    "Missing primary Step89A2 file(s): ",
    paste(missing_primary$file,collapse=" ; ")
  )
}

read_one <- function(project,rank){

  f <- get_file(project,rank)

  if(!file.exists(f)) return(tibble())

  x <- read_csv(
    f,
    show_col_types=FALSE,
    progress=FALSE
  )

  needed <- c(
    "subset_type",
    "required_anchor",
    "early_label",
    "late_label",
    "taxon",
    "n_pairs",
    "prevalence_early",
    "prevalence_late",
    "mean_clr_difference"
  )

  miss <- setdiff(needed,names(x))

  if(length(miss)>0){
    stop(
      project," ",rank,
      " missing required columns: ",
      paste(miss,collapse=", ")
    )
  }

  x %>%
    mutate(
      project=project,
      tax_rank=rank,
      source_file=f
    )
}

raw_list <- list()

for(prj in all_projects){
  for(rk in ranks){
    x <- read_one(prj,rk)
    if(nrow(x)>0){
      raw_list[[paste(prj,rk,sep="__")]] <- x
    }
  }
}

raw <- bind_rows(raw_list)

comparison_audit <- raw %>%
  distinct(
    project,
    tax_rank,
    subset_type,
    required_anchor,
    early_label,
    late_label,
    source_file
  ) %>%
  mutate(
    role=case_when(
      subset_type=="ALL_PAIRED" ~ "PRIMARY",
      subset_type=="COMMON_ANCHOR_SENSITIVITY" ~ "SENSITIVITY",
      TRUE ~ "OTHER_NOT_USED"
    ),
    selected=role %in% c("PRIMARY","SENSITIVITY"),
    selection_reason=case_when(
      role=="PRIMARY" ~
        "Frozen primary rule: ALL_PAIRED",
      role=="SENSITIVITY" ~
        "Frozen sensitivity rule: COMMON_ANCHOR_SENSITIVITY",
      TRUE ~
        "Not part of prespecified Step36I3 comparison"
    )
  )

write_csv(
  comparison_audit,
  file.path(OUT,"01_COMPARISON_SELECTION_AUDIT.csv")
)

# Hard checks for primary 3 cohorts:
# each rank must have exactly one PRIMARY and one SENSITIVITY definition.
check_defs <- comparison_audit %>%
  filter(project %in% primary_projects) %>%
  count(project,tax_rank,role)

for(prj in primary_projects){
  for(rk in ranks){

    np <- check_defs %>%
      filter(
        project==prj,
        tax_rank==rk,
        role=="PRIMARY"
      ) %>%
      pull(n)

    ns <- check_defs %>%
      filter(
        project==prj,
        tax_rank==rk,
        role=="SENSITIVITY"
      ) %>%
      pull(n)

    if(length(np)!=1 || np!=1){
      stop(prj," ",rk," does not have exactly one ALL_PAIRED definition.")
    }

    if(length(ns)!=1 || ns!=1){
      stop(prj," ",rk," does not have exactly one COMMON_ANCHOR_SENSITIVITY definition.")
    }
  }
}

# ------------------------------------------------------------
# 2. Harmonize within-cohort effect tables
# ------------------------------------------------------------
effects <- raw %>%
  filter(
    subset_type %in% unname(subset_definitions)
  ) %>%
  mutate(
    analysis_role=case_when(
      subset_type=="ALL_PAIRED" ~ "PRIMARY",
      subset_type=="COMMON_ANCHOR_SENSITIVITY" ~ "SENSITIVITY"
    ),
    canonical_taxon=map2_chr(
      taxon,
      tax_rank,
      canonicalize_taxon
    ),
    named_taxon=is_named(canonical_taxon),
    effect=as.numeric(mean_clr_difference),
    prevalence_max=pmax(
      as.numeric(prevalence_early),
      as.numeric(prevalence_late),
      na.rm=TRUE
    ),
    eligible=
      is.finite(effect) &
      as.numeric(n_pairs)>=MIN_PAIRS &
      prevalence_max>=MIN_PREVALENCE
  )

write_csv(
  effects,
  file.path(OUT,"02_HARMONIZED_WITHIN_COHORT_EFFECTS.csv")
)

# ------------------------------------------------------------
# 3. Pairwise concordance function
# ------------------------------------------------------------
pairwise_concordance <- function(
  project_a,
  project_b,
  rank,
  analysis_role,
  named_only=TRUE
){

  a <- effects %>%
    filter(
      project==project_a,
      tax_rank==rank,
      .data$analysis_role==analysis_role,
      eligible
    )

  b <- effects %>%
    filter(
      project==project_b,
      tax_rank==rank,
      .data$analysis_role==analysis_role,
      eligible
    )

  if(named_only){
    a <- a %>% filter(named_taxon)
    b <- b %>% filter(named_taxon)
  }

  a <- a %>%
    group_by(canonical_taxon) %>%
    summarise(
      effect_a=mean(effect,na.rm=TRUE),
      .groups="drop"
    )

  b <- b %>%
    group_by(canonical_taxon) %>%
    summarise(
      effect_b=mean(effect,na.rm=TRUE),
      .groups="drop"
    )

  m <- inner_join(
    a,b,
    by="canonical_taxon"
  ) %>%
    filter(
      is.finite(effect_a),
      is.finite(effect_b)
    )

  n <- nrow(m)

  if(n<5){

    return(tibble(
      project_a=project_a,
      project_b=project_b,
      rank=rank,
      analysis_role=analysis_role,
      named_only=named_only,
      n_shared=n,
      spearman_rho=NA_real_,
      spearman_p=NA_real_,
      sign_agreement=NA_real_,
      top_k=NA_integer_,
      top_jaccard=NA_real_,
      permutation_p_rho=NA_real_,
      permutation_p_sign=NA_real_,
      permutation_p_top_jaccard=NA_real_
    ))
  }

  ct <- suppressWarnings(
    cor.test(
      m$effect_a,
      m$effect_b,
      method="spearman",
      exact=FALSE
    )
  )

  rho_obs <- unname(ct$estimate)
  p_rho <- ct$p.value

  nz <- sign(m$effect_a)!=0 &
        sign(m$effect_b)!=0

  sign_obs <- if(sum(nz)>0){
    mean(
      sign(m$effect_a[nz]) ==
      sign(m$effect_b[nz])
    )
  }else{
    NA_real_
  }

  k <- min(
    n,
    max(
      5,
      min(
        20,
        ceiling(0.20*n)
      )
    )
  )

  top_a <- m %>%
    arrange(desc(abs(effect_a))) %>%
    slice_head(n=k) %>%
    pull(canonical_taxon)

  top_b <- m %>%
    arrange(desc(abs(effect_b))) %>%
    slice_head(n=k) %>%
    pull(canonical_taxon)

  jac_obs <-
    length(intersect(top_a,top_b)) /
    length(union(top_a,top_b))

  perm_rho <- numeric(N_PERM)
  perm_sign <- numeric(N_PERM)
  perm_jac <- numeric(N_PERM)

  for(i in seq_len(N_PERM)){

    eb <- sample(
      m$effect_b,
      replace=FALSE
    )

    perm_rho[i] <- suppressWarnings(
      cor(
        m$effect_a,
        eb,
        method="spearman"
      )
    )

    nzp <- sign(m$effect_a)!=0 &
           sign(eb)!=0

    perm_sign[i] <- if(sum(nzp)>0){
      mean(
        sign(m$effect_a[nzp]) ==
        sign(eb[nzp])
      )
    }else{
      NA_real_
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

    perm_jac[i] <-
      length(intersect(ta,tb)) /
      length(union(ta,tb))
  }

  perm_p_rho <-
    (1+sum(abs(perm_rho)>=abs(rho_obs),na.rm=TRUE)) /
    (1+sum(is.finite(perm_rho)))

  perm_p_sign <-
    (1+sum(perm_sign>=sign_obs,na.rm=TRUE)) /
    (1+sum(is.finite(perm_sign)))

  perm_p_jac <-
    (1+sum(perm_jac>=jac_obs,na.rm=TRUE)) /
    (1+sum(is.finite(perm_jac)))

  tibble(
    project_a=project_a,
    project_b=project_b,
    rank=rank,
    analysis_role=analysis_role,
    named_only=named_only,
    n_shared=n,
    spearman_rho=rho_obs,
    spearman_p=p_rho,
    sign_agreement=sign_obs,
    top_k=k,
    top_jaccard=jac_obs,
    permutation_p_rho=perm_p_rho,
    permutation_p_sign=perm_p_sign,
    permutation_p_top_jaccard=perm_p_jac
  )
}

primary_pairs <- combn(
  primary_projects,
  2,
  simplify=FALSE
)

# ------------------------------------------------------------
# 4. Primary + common-anchor sensitivity
# ------------------------------------------------------------
primary_pairwise <- bind_rows(
  lapply(
    c("PRIMARY","SENSITIVITY"),
    function(role){
      bind_rows(
        lapply(
          ranks,
          function(rk){
            bind_rows(
              lapply(
                primary_pairs,
                function(pp){
                  pairwise_concordance(
                    pp[1],pp[2],
                    rk,role,TRUE
                  )
                }
              )
            )
          }
        )
      )
    }
  )
)

write_csv(
  primary_pairwise,
  file.path(OUT,"03_PRIMARY3_PAIRWISE_CONCORDANCE.csv")
)

# all-taxa sensitivity
alltaxa_pairwise <- bind_rows(
  lapply(
    c("PRIMARY","SENSITIVITY"),
    function(role){
      bind_rows(
        lapply(
          ranks,
          function(rk){
            bind_rows(
              lapply(
                primary_pairs,
                function(pp){
                  pairwise_concordance(
                    pp[1],pp[2],
                    rk,role,FALSE
                  )
                }
              )
            )
          }
        )
      )
    }
  )
)

write_csv(
  alltaxa_pairwise,
  file.path(OUT,"04_ALL_TAXA_SENSITIVITY_PAIRWISE_CONCORDANCE.csv")
)

# ------------------------------------------------------------
# 5. Three-way shared taxa
# ------------------------------------------------------------
threeway <- list()

for(role in c("PRIMARY","SENSITIVITY")){
  for(rk in ranks){

    xs <- lapply(
      primary_projects,
      function(prj){

        effects %>%
          filter(
            project==prj,
            tax_rank==rk,
            analysis_role==role,
            eligible,
            named_taxon
          ) %>%
          group_by(canonical_taxon) %>%
          summarise(
            effect=mean(effect,na.rm=TRUE),
            .groups="drop"
          ) %>%
          rename(!!prj := effect)
      }
    )

    m <- reduce(
      xs,
      inner_join,
      by="canonical_taxon"
    )

    if(nrow(m)>0){

      mat <- as.matrix(
        m[,primary_projects,drop=FALSE]
      )

      sm <- sign(mat)

      m$all_same_direction <- apply(
        sm,
        1,
        function(z){
          z <- z[z!=0]
          length(z)==3 &&
            length(unique(z))==1
        }
      )

      m$n_positive <- rowSums(sm>0)
      m$n_negative <- rowSums(sm<0)

      m$rank <- rk
      m$analysis_role <- role

      threeway[[paste(role,rk)]] <- m
    }
  }
}

threeway_df <- bind_rows(threeway)

write_csv(
  threeway_df,
  file.path(OUT,"05_PRIMARY3_SHARED_TAXA_DIRECTION_PATTERNS.csv")
)

if(nrow(threeway_df)>0){

  threeway_summary <- threeway_df %>%
    group_by(
      analysis_role,
      rank
    ) %>%
    summarise(
      n_taxa_shared_all3=n(),
      n_same_direction_all3=
        sum(all_same_direction),
      proportion_same_direction_all3=
        mean(all_same_direction),
      n_all_positive=
        sum(n_positive==3),
      n_all_negative=
        sum(n_negative==3),
      .groups="drop"
    )

}else{

  threeway_summary <- tibble(
    analysis_role=character(),
    rank=character(),
    n_taxa_shared_all3=integer(),
    n_same_direction_all3=integer(),
    proportion_same_direction_all3=double(),
    n_all_positive=integer(),
    n_all_negative=integer()
  )
}

write_csv(
  threeway_summary,
  file.path(OUT,"06_THREEWAY_DIRECTION_SUMMARY.csv")
)

# ------------------------------------------------------------
# 6. Contextual external/control comparisons
# ------------------------------------------------------------
context_pairs <- list(
  c("PRJNA691455","PRJEB82425"),
  c("PRJNA851469","PRJEB82425"),
  c("PRJNA516701","PRJEB82425"),
  c("PRJNA691455","PRJNA578267"),
  c("PRJNA851469","PRJNA578267"),
  c("PRJNA516701","PRJNA578267")
)

context_results <- list()

for(role in c("PRIMARY","SENSITIVITY")){
  for(rk in ranks){
    for(pp in context_pairs){

      needed <- registry %>%
        filter(
          project %in% pp,
          rank==rk
        )

      if(
        nrow(needed)==2 &&
        all(needed$exists)
      ){

        rr <- pairwise_concordance(
          pp[1],pp[2],
          rk,role,TRUE
        )

        context_results[[
          paste(
            role,rk,
            pp[1],pp[2],
            sep="__"
          )
        ]] <- rr
      }
    }
  }
}

context_df <- bind_rows(context_results)

write_csv(
  context_df,
  file.path(OUT,"07_CONTEXT_EXTERNAL_CONTROL_CONCORDANCE.csv")
)

# ------------------------------------------------------------
# 7. Primary synthesis
# ------------------------------------------------------------
primary_summary <- primary_pairwise %>%
  filter(analysis_role=="PRIMARY") %>%
  group_by(rank) %>%
  summarise(
    n_pairwise_comparisons=n(),
    median_spearman_rho=
      median(spearman_rho,na.rm=TRUE),
    min_spearman_rho=
      min(spearman_rho,na.rm=TRUE),
    max_spearman_rho=
      max(spearman_rho,na.rm=TRUE),
    median_sign_agreement=
      median(sign_agreement,na.rm=TRUE),
    median_top_jaccard=
      median(top_jaccard,na.rm=TRUE),
    pairwise_rho_perm_significant=
      sum(permutation_p_rho<0.05,na.rm=TRUE),
    pairwise_sign_perm_significant=
      sum(permutation_p_sign<0.05,na.rm=TRUE),
    pairwise_jaccard_perm_significant=
      sum(permutation_p_top_jaccard<0.05,na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(
    ecological_direction_concordance=1.0,
    taxonomic_concordance_class=case_when(
      abs(median_spearman_rho)<0.20 ~
        "LOW",
      abs(median_spearman_rho)<0.50 ~
        "MODEST",
      TRUE ~
        "HIGH"
    )
  )

write_csv(
  primary_summary,
  file.path(OUT,"08_ECOLOGICAL_VS_TAXONOMIC_SYNTHESIS.csv")
)

sens_summary <- primary_pairwise %>%
  filter(analysis_role=="SENSITIVITY") %>%
  group_by(rank) %>%
  summarise(
    median_spearman_rho=
      median(spearman_rho,na.rm=TRUE),
    median_sign_agreement=
      median(sign_agreement,na.rm=TRUE),
    median_top_jaccard=
      median(top_jaccard,na.rm=TRUE),
    .groups="drop"
  )

write_csv(
  sens_summary,
  file.path(OUT,"09_COMMON_ANCHOR_SENSITIVITY_SUMMARY.csv")
)

# ------------------------------------------------------------
# 8. Evidence tier
# ------------------------------------------------------------
classes <- primary_summary$taxonomic_concordance_class

tier <- if(
  nrow(primary_summary)==2 &&
  all(classes %in% c("LOW","MODEST"))
){
  "ECOLOGICAL_DIRECTION_MORE_REPRODUCIBLE_THAN_TAXONOMIC_EFFECTS"
}else{
  "TAXONOMIC_CONCORDANCE_NOT_CLEARLY_LOWER_THAN_ECOLOGICAL_DIRECTION"
}

decision <- tibble(
  primary_ecological_direction_concordance=1.0,
  genus_taxonomic_class=
    primary_summary$taxonomic_concordance_class[
      primary_summary$rank=="GENUS"
    ],
  family_taxonomic_class=
    primary_summary$taxonomic_concordance_class[
      primary_summary$rank=="FAMILY"
    ],
  evidence_tier=tier
)

write_csv(
  decision,
  file.path(OUT,"10_EVIDENCE_DECISION.csv")
)

# ------------------------------------------------------------
# 9. Candidate figure
# ------------------------------------------------------------
plot_df <- primary_pairwise %>%
  mutate(
    pair=paste(
      project_a,
      project_b,
      sep=" vs "
    ),
    analysis_role=factor(
      analysis_role,
      levels=c("PRIMARY","SENSITIVITY")
    )
  )

p <- ggplot(
  plot_df,
  aes(
    x=pair,
    y=spearman_rho,
    shape=rank
  )
) +
  geom_hline(
    yintercept=0,
    linetype=2
  ) +
  geom_point(size=3) +
  facet_wrap(~analysis_role) +
  coord_flip() +
  labs(
    x=NULL,
    y="Spearman correlation of within-cohort CLR effect vectors",
    title="Cross-cohort taxonomic trajectory reproducibility"
  ) +
  theme_classic(base_size=10)

ggsave(
  file.path(
    OUT,
    "11_Figure_candidate_taxonomic_concordance.pdf"
  ),
  p,
  width=8,
  height=5.5
)

ggsave(
  file.path(
    OUT,
    "11_Figure_candidate_taxonomic_concordance.png"
  ),
  p,
  width=8,
  height=5.5,
  dpi=300
)

# ------------------------------------------------------------
# 10. Manuscript-safe summary
# ------------------------------------------------------------
fmt <- function(x,d=3){
  ifelse(
    is.na(x),
    "NA",
    formatC(
      x,
      format="f",
      digits=d
    )
  )
}

lines <- c(
  "V2 STEP36I3 CROSS-COHORT TAXONOMIC REPRODUCIBILITY",
  "",
  "Comparison selection:",
  "ALL_PAIRED = primary.",
  "COMMON_ANCHOR_SENSITIVITY = prespecified sensitivity.",
  "",
  "Primary cohorts:",
  paste(primary_projects,collapse=", "),
  "",
  "Ecological displacement direction: 3/3 primary cohorts progressive (direction concordance = 1.00).",
  ""
)

for(i in seq_len(nrow(primary_summary))){

  z <- primary_summary[i,]

  lines <- c(
    lines,
    paste0(
      z$rank,
      " primary taxonomic effects: median pairwise rho=",
      fmt(z$median_spearman_rho),
      "; median sign agreement=",
      fmt(z$median_sign_agreement),
      "; median top-effect Jaccard=",
      fmt(z$median_top_jaccard),
      "; class=",
      z$taxonomic_concordance_class,
      "."
    )
  )
}

if(nrow(threeway_summary)>0){

  for(i in seq_len(nrow(threeway_summary))){

    z <- threeway_summary[i,]

    lines <- c(
      lines,
      paste0(
        z$analysis_role,
        " ",
        z$rank,
        " three-way shared taxa: n=",
        z$n_taxa_shared_all3,
        "; all-three same-direction proportion=",
        fmt(
          z$proportion_same_direction_all3
        ),
        "."
      )
    )
  }
}

lines <- c(
  lines,
  "",
  paste0(
    "EVIDENCE TIER: ",
    tier
  ),
  "",
  "Interpretation guardrails:",
  "1. Taxonomic effects were estimated within cohort; no ASV/OTU matrices were pooled.",
  "2. Ecological and taxonomic reproducibility are different biological resolutions and are not numerically interchangeable.",
  "3. Do not claim a universal sepsis-specific taxonomic signature.",
  "4. If taxonomic concordance is low/modest, frame this as heterogeneous taxonomic routes accompanying a more reproducible ecological displacement direction.",
  "5. COMMON_ANCHOR_SENSITIVITY is supportive sensitivity, not a second independent discovery analysis.",
  "6. No previously frozen V2 inference was modified."
)

writeLines(
  lines,
  file.path(
    OUT,
    "12_MANUSCRIPT_SAFE_SUMMARY.txt"
  )
)

writeLines(
  c(
    paste0(
      "Completed: ",
      Sys.time()
    ),
    "STEP36I3 COMPLETE",
    paste0(
      "Evidence tier: ",
      tier
    ),
    "No previously frozen V2 inference was modified."
  ),
  file.path(
    OUT,
    "_STEP36I3_COMPLETE.txt"
  )
)

cat("STEP36I3 COMPLETE\n")
cat("Evidence tier:",tier,"\n")
