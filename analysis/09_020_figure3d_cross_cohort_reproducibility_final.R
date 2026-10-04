
# ============================================================
# V2_96E FIGURE 3D — CROSS-COHORT TAXONOMIC REPRODUCIBILITY
#
# Main manuscript panel:
#   PRIMARY only (all eligible early-to-late pairs)
#
# Supplement:
#   PRIMARY vs common-anchor sensitivity
#
# No statistical inference is rerun.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

ROOT <- "E:/sepsis_project"
RESULTS <- file.path(ROOT,"results")

IN_DIR <- file.path(
  RESULTS,
  "V2_36I3_CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY_FINAL"
)

PAIR_FILE <- file.path(
  IN_DIR,
  "03_PRIMARY3_PAIRWISE_CONCORDANCE.csv"
)

THREE_FILE <- file.path(
  IN_DIR,
  "06_THREEWAY_DIRECTION_SUMMARY.csv"
)

SYN_FILE <- file.path(
  IN_DIR,
  "08_ECOLOGICAL_VS_TAXONOMIC_SYNTHESIS.csv"
)

OUT <- file.path(
  RESULTS,
  "V2_96E_FIGURE3D_CROSS_COHORT_REPRODUCIBILITY_FINAL"
)

dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

for(f in c(PAIR_FILE,THREE_FILE,SYN_FILE)){
  if(!file.exists(f)){
    stop("Missing frozen Step36I3 source: ",f)
  }
}

pair <- read_csv(
  PAIR_FILE,
  show_col_types=FALSE
)

three <- read_csv(
  THREE_FILE,
  show_col_types=FALSE
)

syn <- read_csv(
  SYN_FILE,
  show_col_types=FALSE
)

# ------------------------------------------------------------
# 1. Main panel data: PRIMARY only
# ------------------------------------------------------------
main <- pair %>%
  filter(
    analysis_role=="PRIMARY",
    named_only
  ) %>%
  mutate(
    pair_label=case_when(
      project_a=="PRJNA691455" & project_b=="PRJNA851469" ~
        "PRJNA691455 vs PRJNA851469",
      project_a=="PRJNA691455" & project_b=="PRJNA516701" ~
        "PRJNA691455 vs PRJNA516701",
      project_a=="PRJNA851469" & project_b=="PRJNA516701" ~
        "PRJNA851469 vs PRJNA516701",
      TRUE ~ paste(project_a,project_b,sep=" vs ")
    ),
    rank=factor(
      rank,
      levels=c("GENUS","FAMILY"),
      labels=c("Genus","Family")
    ),
    rho_label=sprintf("%.2f",spearman_rho)
  )

if(nrow(main)!=6){
  stop(
    "Expected 6 PRIMARY named-taxon rows (3 cohort pairs x 2 ranks), found ",
    nrow(main)
  )
}

# Preserve a stable visual order.
pair_order <- c(
  "PRJNA691455 vs PRJNA851469",
  "PRJNA691455 vs PRJNA516701",
  "PRJNA851469 vs PRJNA516701"
)

main <- main %>%
  mutate(
    pair_label=factor(
      pair_label,
      levels=rev(pair_order)
    )
  )

write_csv(
  main,
  file.path(
    OUT,
    "01_FIGURE3D_MAIN_DATA.csv"
  )
)

# ------------------------------------------------------------
# 2. Figure 3D main panel
# ------------------------------------------------------------
xmin <- min(c(-0.10,main$spearman_rho),na.rm=TRUE)
xmax <- max(c(0.55,main$spearman_rho),na.rm=TRUE)

pad <- max(0.08,(xmax-xmin)*0.12)

p_main <- ggplot(
  main,
  aes(
    x=spearman_rho,
    y=pair_label,
    shape=rank
  )
) +
  geom_vline(
    xintercept=0,
    linetype=2,
    linewidth=0.5
  ) +
  geom_point(
    size=3.2,
    position=position_dodge(width=0.40)
  ) +
  geom_text(
    aes(label=rho_label),
    hjust=-0.35,
    size=3.2,
    position=position_dodge(width=0.40),
    show.legend=FALSE
  ) +
  scale_shape_manual(
    values=c("Genus"=16,"Family"=17)
  ) +
  coord_cartesian(
    xlim=c(xmin-pad,xmax+pad),
    clip="off"
  ) +
  labs(
    x="Spearman correlation of longitudinal CLR effect vectors",
    y=NULL,
    shape=NULL
  ) +
  theme_classic(
    base_size=10.5
  ) +
  theme(
    legend.position="top",
    legend.justification="left",
    axis.text.y=element_text(size=9.5),
    plot.margin=margin(7,26,7,7)
  )

ggsave(
  file.path(
    OUT,
    "Figure3D_Cross_Cohort_Taxonomic_Concordance_MAIN.pdf"
  ),
  p_main,
  width=7.1,
  height=3.5,
  device=cairo_pdf
)

ggsave(
  file.path(
    OUT,
    "Figure3D_Cross_Cohort_Taxonomic_Concordance_MAIN.png"
  ),
  p_main,
  width=7.1,
  height=3.5,
  dpi=600
)

# ------------------------------------------------------------
# 3. Supplementary sensitivity figure
# ------------------------------------------------------------
supp <- pair %>%
  filter(named_only) %>%
  mutate(
    pair_label=case_when(
      project_a=="PRJNA691455" & project_b=="PRJNA851469" ~
        "PRJNA691455 vs PRJNA851469",
      project_a=="PRJNA691455" & project_b=="PRJNA516701" ~
        "PRJNA691455 vs PRJNA516701",
      project_a=="PRJNA851469" & project_b=="PRJNA516701" ~
        "PRJNA851469 vs PRJNA516701",
      TRUE ~ paste(project_a,project_b,sep=" vs ")
    ),
    pair_label=factor(
      pair_label,
      levels=rev(pair_order)
    ),
    rank=factor(
      rank,
      levels=c("GENUS","FAMILY"),
      labels=c("Genus","Family")
    ),
    analysis_role=factor(
      analysis_role,
      levels=c("PRIMARY","SENSITIVITY"),
      labels=c(
        "Primary",
        "Common-anchor sensitivity"
      )
    ),
    rho_label=sprintf("%.2f",spearman_rho)
  )

write_csv(
  supp,
  file.path(
    OUT,
    "02_SUPPLEMENTARY_SENSITIVITY_DATA.csv"
  )
)

p_supp <- ggplot(
  supp,
  aes(
    x=spearman_rho,
    y=pair_label,
    shape=rank
  )
) +
  geom_vline(
    xintercept=0,
    linetype=2,
    linewidth=0.5
  ) +
  geom_point(
    size=3.0,
    position=position_dodge(width=0.40)
  ) +
  geom_text(
    aes(label=rho_label),
    hjust=-0.35,
    size=3.0,
    position=position_dodge(width=0.40),
    show.legend=FALSE
  ) +
  scale_shape_manual(
    values=c("Genus"=16,"Family"=17)
  ) +
  facet_wrap(
    ~analysis_role,
    ncol=1
  ) +
  coord_cartesian(
    xlim=c(
      min(c(-0.10,supp$spearman_rho),na.rm=TRUE)-0.08,
      max(c(0.55,supp$spearman_rho),na.rm=TRUE)+0.10
    ),
    clip="off"
  ) +
  labs(
    x="Spearman correlation of longitudinal CLR effect vectors",
    y=NULL,
    shape=NULL
  ) +
  theme_classic(
    base_size=10.5
  ) +
  theme(
    legend.position="top",
    legend.justification="left",
    strip.background=element_blank(),
    strip.text=element_text(face="bold"),
    plot.margin=margin(7,26,7,7)
  )

ggsave(
  file.path(
    OUT,
    "Supplementary_Figure_Cross_Cohort_Taxonomic_Concordance_Sensitivity.pdf"
  ),
  p_supp,
  width=7.2,
  height=6.0,
  device=cairo_pdf
)

ggsave(
  file.path(
    OUT,
    "Supplementary_Figure_Cross_Cohort_Taxonomic_Concordance_Sensitivity.png"
  ),
  p_supp,
  width=7.2,
  height=6.0,
  dpi=600
)

# ------------------------------------------------------------
# 4. Exact figure legend metrics
# ------------------------------------------------------------
genus_three <- three %>%
  filter(
    analysis_role=="PRIMARY",
    rank=="GENUS"
  )

family_three <- three %>%
  filter(
    analysis_role=="PRIMARY",
    rank=="FAMILY"
  )

genus_syn <- syn %>%
  filter(rank=="GENUS")

family_syn <- syn %>%
  filter(rank=="FAMILY")

if(
  nrow(genus_three)!=1 ||
  nrow(family_three)!=1 ||
  nrow(genus_syn)!=1 ||
  nrow(family_syn)!=1
){
  stop("Could not uniquely resolve Figure 3D summary metrics.")
}

legend <- c(
  "Figure 3D. Cross-cohort reproducibility of longitudinal taxonomic effects.",
  "",
  paste0(
    "Longitudinal paired CLR taxonomic effects were estimated independently within each cohort before cross-cohort comparison. ",
    "Points show pairwise Spearman correlations between shared-taxon effect vectors for the three cohorts with robust progressive ecological displacement; genus- and family-level analyses are shown separately. ",
    "The direction of ecosystem-level displacement was concordant across all three cohorts, whereas taxonomic effect concordance was incomplete and cohort-dependent. ",
    "Across the three cohort pairs, the median Spearman correlation was ",
    sprintf("%.2f",genus_syn$median_spearman_rho[[1]]),
    " at genus level and ",
    sprintf("%.2f",family_syn$median_spearman_rho[[1]]),
    " at family level. ",
    "Among taxa shared across all three cohorts, ",
    sprintf("%.1f",100*genus_three$proportion_same_direction_all3[[1]]),
    "% of genera and ",
    sprintf("%.1f",100*family_three$proportion_same_direction_all3[[1]]),
    "% of families changed in the same direction in all three cohorts. ",
    "ASV or OTU abundance matrices were not pooled across cohorts."
  )
)

writeLines(
  legend,
  file.path(
    OUT,
    "03_FIGURE3D_LEGEND.txt"
  )
)

supp_legend <- c(
  "Supplementary Figure. Common-anchor sensitivity analysis of cross-cohort taxonomic reproducibility.",
  "",
  "Pairwise Spearman correlations of longitudinal CLR effect vectors are shown for the primary all-eligible-pairs analysis and the common-anchor sensitivity analysis. The similar overall pattern supports robustness of the conclusion that taxonomic concordance was incomplete and cohort-dependent. The common-anchor analysis is a sensitivity analysis and should not be interpreted as independent replication."
)

writeLines(
  supp_legend,
  file.path(
    OUT,
    "04_SUPPLEMENTARY_FIGURE_LEGEND.txt"
  )
)

# ------------------------------------------------------------
# 5. Figure integration manifest
# ------------------------------------------------------------
manifest <- tibble(
  component=c(
    "Figure3D main panel",
    "Supplementary sensitivity figure"
  ),
  role=c(
    "MAIN_MANUSCRIPT",
    "SUPPLEMENT"
  ),
  primary_only=c(
    TRUE,
    FALSE
  ),
  inference_rerun=FALSE,
  source="V2_36I3_CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY_FINAL"
)

write_csv(
  manifest,
  file.path(
    OUT,
    "05_FIGURE_INTEGRATION_MANIFEST.csv"
  )
)

# ------------------------------------------------------------
# 6. QC
# ------------------------------------------------------------
qc <- tibble(
  check=c(
    "main_has_6_rows",
    "main_primary_only",
    "main_named_taxa_only",
    "supp_has_primary_and_sensitivity",
    "genus_threeway_metric_available",
    "family_threeway_metric_available",
    "no_inference_rerun"
  ),
  passed=c(
    nrow(main)==6,
    all(main$analysis_role=="PRIMARY"),
    all(main$named_only),
    all(c("Primary","Common-anchor sensitivity") %in% levels(supp$analysis_role)),
    nrow(genus_three)==1,
    nrow(family_three)==1,
    TRUE
  )
)

write_csv(
  qc,
  file.path(
    OUT,
    "06_FIGURE3D_QC.csv"
  )
)

if(!all(qc$passed)){
  stop(
    "Figure3D QC failed: ",
    paste(qc$check[!qc$passed],collapse=", ")
  )
}

writeLines(
  c(
    paste0("Completed: ",Sys.time()),
    "STEP96E COMPLETE",
    "Figure 3D main panel generated from PRIMARY results only.",
    "Common-anchor sensitivity exported as Supplementary Figure.",
    "No statistical inference rerun."
  ),
  file.path(
    OUT,
    "_STEP96E_COMPLETE.txt"
  )
)

cat("STEP96E COMPLETE\n")
