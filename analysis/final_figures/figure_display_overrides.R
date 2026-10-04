# Display-only revisions. Read frozen estimates; do not run inferential tests.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
base_theme <- function() {
  theme_classic(base_size = 9, base_family = "Arial") +
    theme(plot.title = element_text(size = 10, face = "bold", margin = margin(b = 8)),
          axis.title = element_text(size = 9), axis.text = element_text(size = 8, colour = "black"),
          legend.text = element_text(size = 8), legend.title = element_blank(),
          plot.caption = element_text(size = 8, hjust = 0),
          plot.margin = margin(9, 12, 8, 9), legend.position = "bottom")
}
theme_set(base_theme())
tagged <- function(p) p + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 11), plot.tag.position = c(0, 1))
term_pair <- function(x) gsub(" vs ", "\nvs ", x, fixed = TRUE)
mm8 <- 8 / ggplot2::.pt
reader_forest <- function(dd, title, axis_title, labels, est, lo, hi, ptxt = NULL, col = blue) {
  d <- data.frame(label = labels, est = est, lo = lo, hi = hi)
  d$label <- factor(d$label, levels = rev(unique(labels)))
  xr <- range(c(lo, hi, 0), na.rm = TRUE); span <- diff(xr)
  if (span == 0) span <- 1
  p <- ggplot(d, aes(y = label, x = est)) +
    geom_vline(xintercept = 0, linetype = 2, colour = grey, linewidth = .4) +
    geom_segment(aes(x = lo, xend = hi, yend = label), colour = col, linewidth = .6) +
    geom_point(colour = col, size = 2.2) +
    scale_x_continuous(limits = c(xr[1] - .06 * span, xr[2] + ifelse(is.null(ptxt), .08, .65) * span)) +
    labs(title = title, x = axis_title, y = NULL) + base_theme()
  if (!is.null(ptxt)) {
    d$ptxt <- ptxt; d$text_x <- xr[2] + .09 * span
    p <- p + geom_text(data = d, aes(x = text_x, label = ptxt), hjust = 0, size = mm8, family = "Arial")
  }
  p
}
corrected_plot <- function(stem) {
  if (stem == "Figure_1_Study_Architecture") {
    p <- ggplot() + coord_cartesian(xlim = c(0,12), ylim = c(0,9), expand = FALSE) + theme_void(base_family = "Arial")
    boxes <- data.frame(x0=c(3.2,.15,3.15,6.15,9.15,.3,3.1,1),
      x1=c(8.8,2.85,5.85,8.85,11.85,2.75,5.55,11),
      y0=c(7.9,5.35,5.35,5.35,5.35,2.95,2.95,.45),
      y1=c(8.8,7.25,7.25,7.25,7.25,4.45,4.45,1.8),
      col=c(grey,blue,green,orange,purple,blue,blue,grey),
      fill=c("white","#E8F2F8","#E8F5EF","#FCEFE8","#F2ECF6","white","white","white"))
    p <- p + geom_rect(data=boxes,aes(xmin=x0,xmax=x1,ymin=y0,ymax=y1),fill=boxes$fill,colour=boxes$col,linewidth=.5)
    edges <- data.frame(x=c(rep(6,4),1.5,1.5,1.5,4.3,4.5,7.5,10.5),
      xend=c(1.5,4.5,7.5,10.5,1.5,4.3,3.2,4.8,6,7.5,9.2),
      y=c(rep(7.9,4),5.35,5.35,2.95,2.95,5.35,5.35,5.35),
      yend=c(rep(7.25,4),4.45,4.45,rep(1.8,5)))
    p <- p + geom_segment(data=edges,aes(x=x,y=y,xend=xend,yend=yend),colour=grey,linewidth=.4,
      arrow=grid::arrow(length=grid::unit(1.6,"mm")))
    p <- p + annotate("text",x=6,y=8.45,label="Public gut-microbiome cohorts",fontface="bold",size=3.4,family="Arial") +
      annotate("text",x=6,y=8.13,label="n = 11 cohorts",size=3,family="Arial")
    titles <- c("Discovery evidence","Independent\nlongitudinal validation","External cross-sectional\necological-state comparison","Infection-source\nlongitudinal extension")
    details <- c("7 longitudinal cohorts\n1 static supportive cohort","PRJNA1125274",
                 "PRJNA1010969\nControl / Trauma / Sepsis","CRA002354\nPulmonary vs recorded\nnon-pulmonary")
    p <- p + annotate("text",x=c(1.5,4.5,7.5,10.5),y=6.73,label=titles,fontface="bold",size=2.85,family="Arial") +
      annotate("text",x=c(1.5,4.5,7.5,10.5),y=5.96,label=details,size=2.85,family="Arial") +
      annotate("text",x=1.525,y=4.05,label="Longitudinal discovery",fontface="bold",size=2.85,family="Arial") +
      annotate("text",x=1.525,y=3.47,label="n = 7 cohorts\nwithin-patient change",size=2.85,family="Arial") +
      annotate("text",x=4.325,y=4.05,label="Static support",fontface="bold",size=2.85,family="Arial") +
      annotate("text",x=4.325,y=3.47,label="n = 1 cohort\nPRJNA978257",size=2.85,family="Arial") +
      annotate("text",x=6,y=1.37,label="Longitudinal ecological displacement recurred across multiple cohorts",fontface="bold",size=3.1,family="Arial") +
      annotate("text",x=6,y=.91,label="despite heterogeneous taxonomic trajectories",size=3.1,family="Arial")
    return(p)
  }
  if (grepl("^Figure_2", stem)) {
    d <- f2; d$lab <- factor(paste0(d$project, "  (n = ", d$primary_n, ")"),
                            levels = rev(paste0(d$project, "  (n = ", d$primary_n, ")")))
    return(ggplot(d, aes(y = lab)) +
      geom_vline(xintercept = 0, linetype = 2, colour = grey) +
      geom_segment(aes(x = primary_effect_dz, xend = common_anchor_effect_dz, yend = lab),
                   linetype = 3, colour = grey) +
      geom_point(aes(x = primary_effect_dz, colour = "Primary early-to-late", shape = "Primary early-to-late"), size = 2.6) +
      geom_point(aes(x = common_anchor_effect_dz, colour = "Common-anchor sensitivity", shape = "Common-anchor sensitivity"),
                 fill = "white", size = 2.8, stroke = .7) +
      geom_text(data = d[d$primary_significant, ], aes(x = primary_effect_dz, label = "*"), nudge_y = .22,
                colour = blue, size = 4, family = "Arial") +
      scale_colour_manual(name = "Analysis", breaks = c("Primary early-to-late", "Common-anchor sensitivity"),
                          values = c("Primary early-to-late" = blue, "Common-anchor sensitivity" = orange)) +
      scale_shape_manual(name = "Analysis", breaks = c("Primary early-to-late", "Common-anchor sensitivity"),
                         values = c("Primary early-to-late" = 16, "Common-anchor sensitivity" = 21)) +
      labs(x = expression("Paired standardized change in Bray-Curtis displacement ("*d[z]*")"), y = NULL,
           caption = "* Primary BH-FDR < 0.05. Connecting lines link analyses of the same cohort; they are not CIs.") + base_theme())
  }
  if (grepl("^Figure_3", stem)) {
    da <- f3a %>% pivot_longer(c(Day3, Day7), names_to = "day", values_to = "value")
    da$day <- factor(da$day, levels = c("Day3", "Day7"), labels = c("Day 3", "Day 7"))
    a <- ggplot(da, aes(day, value, group = Patient_true)) +
      geom_line(colour = grey, alpha = .6, linewidth = .4) + geom_point(aes(colour = day), size = 2) +
      scale_colour_manual(values = c(blue, orange), guide = "none") +
      labs(title = "Within-patient ecological\ndisplacement", subtitle = "PRJNA691455; 9 paired patients", x = NULL,
           y = "Bray-Curtis displacement") + base_theme()
    b <- ggplot(f3b, aes(x = "", y = signed_Euclidean)) +
      geom_boxplot(width = .35, outlier.shape = NA, fill = "#E8F2F8", colour = blue) +
      geom_jitter(colour = grey, alpha = .65, size = 1.2,
                  position = position_jitter(width = .08, height = 0, seed = 100)) +
      coord_flip() + labs(title = "Taxonomic trajectory\nheterogeneity", subtitle = "45 comparisons; 10 patients",
                          x = NULL, y = "Signed-trajectory Euclidean distance") + base_theme() +
      theme(plot.title = element_text(margin = margin(l = 12, b = 8)),
            plot.subtitle = element_text(margin = margin(l = 12, b = 8)))
    dc <- f3c; dc$label <- c("latest_EII" = "Exploratory EII", "latest_Simpson_instability" = "Simpson-diversity\ninstability",
                           "latest_Bray" = "Latest Bray-Curtis\ndisplacement")[dc$outcome]
    dc$label <- factor(dc$label, levels = rev(c("Exploratory EII", "Simpson-diversity\ninstability", "Latest Bray-Curtis\ndisplacement")))
    c <- ggplot(dc, aes(spearman_rho, label)) + geom_vline(xintercept = 0, linetype = 2, colour = grey) +
      geom_point(colour = blue, size = 2.3) +
      geom_text(aes(x = 1.05, label = paste0("perm. p = ", vapply(permutation_p_10000, fmt_p, character(1)))),
                hjust = 0, size = mm8, family = "Arial") +
      scale_x_continuous(limits = c(-.05, 1.72), breaks = c(0, .25, .5, .75, 1)) +
      labs(title = "Trajectory heterogeneity\nand ecological measures", subtitle = "Patient-level correlations; n = 10",
           x = expression("Spearman "*rho), y = NULL) + base_theme()
    dd <- f3d; dd$lab <- factor(term_pair(dd$pair_label), levels = rev(unique(term_pair(dd$pair_label))))
    d <- ggplot(dd, aes(spearman_rho, lab, colour = rank, shape = rank)) +
      geom_vline(xintercept = 0, linetype = 2, colour = grey) +
      geom_point(position = position_dodge(width = .4), size = 2.3) +
      scale_colour_manual(values = c(Genus = blue, Family = orange)) + scale_shape_manual(values = c(Genus = 16, Family = 17)) +
      scale_x_continuous(limits = c(-.05, .55), breaks = c(0, .2, .4)) +
      labs(title = "Cross-cohort taxonomic\nconcordance", x = expression("Spearman "*rho), y = NULL) + base_theme()
    return(tagged((a | b) / (c | d)))
  }
  if (grepl("^Figure_4", stem)) {
    da <- read_pack("FIGURE4a_sample_distances.csv") %>% filter(Group %in% c("Trauma", "Sepsis"))
    da$Group <- factor(da$Group, levels = c("Trauma", "Sepsis"))
    med <- data.frame(Group = factor(c("Trauma", "Sepsis"), levels = c("Trauma", "Sepsis")),
                      value = c(f4a$median_B[1], f4a$median_A[1]))
    a <- ggplot(da, aes(Group, Bray_to_Control_Centroid, colour = Group)) +
      geom_boxplot(width = .45, outlier.shape = NA, fill = "white", linewidth = .5) +
      geom_jitter(size = 1.6, alpha = .8,
                  position = position_jitter(width = .10, height = 0, seed = 100)) +
      geom_text(data = med, aes(y = 1.03, label = paste0("median = ", fmt(value))), size = mm8, family = "Arial") +
      scale_colour_manual(values = c(Trauma = grey, Sepsis = orange), guide = "none") +
      scale_x_discrete(labels = c(Trauma = "Trauma\nn = 18", Sepsis = "Sepsis\nn = 18")) +
      scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, .25)) +
      labs(title = "Ecological-state differentiation", subtitle = "PRJNA1010969; Control reference n = 17",
           x = NULL, y = "Bray-Curtis distance from control centroid") + base_theme()
    db <- f4b; db$pc <- factor(format(db$pseudocount, scientific = TRUE, digits = 1),
                             levels = format(db$pseudocount, scientific = TRUE, digits = 1))
    b <- ggplot(db, aes(pc, PERMANOVA_R2)) + geom_point(colour = blue, size = 2.5) +
      geom_text(aes(label = fmt(PERMANOVA_R2)), nudge_y = .0012,
                size = mm8, family = "Arial") +
      scale_y_continuous(limits = c(.077, .089), breaks = c(.078, .082, .086)) +
      labs(title = "Aitchison pseudocount\nsensitivity", subtitle = "Sepsis vs Trauma; n = 36 samples",
           x = "Pseudocount", y = expression("PERMANOVA "*R^2)) + base_theme()
    dc <- f4c; dc$label <- factor(c("Longitudinal vs Sepsis-Control", "Longitudinal vs Sepsis-Trauma"),
                               levels = rev(c("Longitudinal vs Sepsis-Control", "Longitudinal vs Sepsis-Trauma")))
    c <- ggplot(dc, aes(spearman_rho, label)) +
      geom_vline(xintercept = 0, linetype = 2, colour = grey) + geom_point(colour = green, size = 2.5) +
      geom_text(aes(x = .43, label = paste0("n = ", overlapping_genera, "; p = ", vapply(p_value, fmt_p, character(1)))),
                hjust = 0, size = mm8, family = "Arial") +
      scale_x_continuous(limits = c(-.04, .66), breaks = c(0, .1, .2, .3, .4)) +
      labs(title = "Exploratory taxonomic alignment", x = expression("Spearman "*rho), y = NULL) + base_theme()
    b <- b + labs(caption = "PERMANOVA FDR: 0.0002-0.0005\nDispersion FDR: 0.168-0.234")
    return(tagged((a | b) / c + plot_layout(heights = c(2, 1))))
  }
  if (grepl("^Figure_5", stem)) {
    a <- reader_forest(f5a, "Primary Bray-Curtis\ninteraction", expression("Source-by-time interaction "*beta),
       "Pulmonary vs recorded\nnon-pulmonary", f5a$effect, f5a$ci_low, f5a$ci_high,
       paste0("p = ", vapply(f5a$p_value, fmt_p, character(1))), orange)
    b <- reader_forest(f5b, "Source-specific slopes", "Bray-Curtis displacement slope per day",
       paste0(gsub("Recorded non-pulmonary", "Recorded\nnon-pulmonary", f5b$group_label), "\n", f5b$npatients,
              " patients\n", f5b$nobs, " observations"), f5b$slope_per_day, f5b$ci95_low, f5b$ci95_high,
       paste0("p = ", vapply(f5b$p_LRT, fmt_p, character(1))))
    smap <- c(ADJUSTED_SOFA = "SOFA adjustment", ADJUSTED_APACHE = "APACHE adjustment",
              ADJUSTED_LACTATE = "Lactate adjustment", EXCLUDE_OTHER_UNKNOWN = "Exclude\nOTHER_UNKNOWN",
              RAREFIED4000_BRAY = "Rarefied 4,000", SENSITIVITY_MIN2000 = "Minimum depth 2,000")
    sl <- unname(smap[f5c$sensitivity]); stopifnot(!anyNA(sl))
    c <- reader_forest(f5c, "Sensitivity analyses", expression("Source-by-time interaction "*beta), sl,
                      f5c$beta, f5c$ci_low, f5c$ci_high, paste0("p = ", vapply(f5c$p_value, fmt_p, character(1))), green)
    dd <- data.frame(label = c("Unrestricted", "ICU Day 3 or earlier"),
       r2 = c(f5d$primary_PERMANOVA_R2[1], f5d$day3_restricted_PERMANOVA_R2[1]),
       p = c(f5d$primary_PERMANOVA_p[1], f5d$day3_restricted_PERMANOVA_p[1]),
       disp = c(f5d$primary_dispersion_p[1], f5d$day3_restricted_dispersion_p[1]))
    dd$label <- factor(dd$label, levels = dd$label)
    d <- ggplot(dd, aes(label, r2)) + geom_point(colour = blue, size = 2.5) +
      geom_text(aes(y = .040, label = paste0("PERMANOVA\np = ", vapply(p, fmt_p, character(1)),
                                           "\nDispersion\np = ", vapply(disp, fmt_p, character(1)))),
                size = mm8, family = "Arial", lineheight = 1.1) +
      scale_y_continuous(limits = c(0, .065), breaks = c(0, .02, .04, .06)) +
      scale_x_discrete(labels = c("Unrestricted", "ICU Day 3\nor earlier")) +
      labs(title = "Baseline beta-diversity\nrobustness", y = expression("PERMANOVA "*R^2), x = NULL) + base_theme()
    return(tagged((a | b) / (c | d)) + plot_annotation(caption = "Intervals in panels a-c: 95% CIs. p values: likelihood-ratio tests (a-c), permutation tests (d)."))
  }
  if (grepl("Supplementary_Figure_S1_", stem)) {
    build <- function(role) {
      d <- s1 %>% filter(analysis_role == role)
      d$lab <- factor(term_pair(d$pair_label), levels = rev(unique(term_pair(d$pair_label))))
      ggplot(d, aes(spearman_rho, lab, colour = rank, shape = rank)) +
        geom_vline(xintercept = 0, linetype = 2, colour = grey) +
        geom_point(position = position_dodge(.35), size = 2.6) +
        scale_colour_manual(values = c(Genus = blue, Family = orange)) +
        scale_shape_manual(values = c(Genus = 16, Family = 17)) +
        labs(title = role, x = expression("Spearman "*rho*" of longitudinal CLR effect vectors"), y = NULL) + base_theme()
    }
    return(tagged(build("Primary") / build("Common-anchor sensitivity")))
  }
  if (grepl("Supplementary_Figure_S2_", stem)) {
    d <- s2; d$lab <- factor(d$project, levels = rev(d$project))
    return(ggplot(d, aes(y = lab)) +
      geom_segment(aes(x = primary_n, xend = common_anchor_n, yend = lab), colour = light_grey) +
      geom_point(aes(primary_n, colour = "Primary", shape = "Primary"), size = 2.7) +
      geom_point(aes(common_anchor_n, colour = "Common-anchor sensitivity", shape = "Common-anchor sensitivity"),
                 fill = "white", size = 2.7) +
      scale_colour_manual(name = "Analysis", breaks = c("Primary", "Common-anchor sensitivity"),
                          values = c(Primary = blue, "Common-anchor sensitivity" = orange)) +
      scale_shape_manual(name = "Analysis", breaks = c("Primary", "Common-anchor sensitivity"),
                         values = c(Primary = 16, "Common-anchor sensitivity" = 21)) +
      labs(title = "Paired patient counts by analysis", x = "Number of paired patients", y = NULL) + base_theme())
  }
  if (grepl("Supplementary_Figure_S[34]_", stem)) {
    is_bray <- grepl("S3_", stem); d <- if (is_bray) s3 else s4
    pool <- s34 %>% filter(scheme == ifelse(is_bray, "bray", "aitchison_CZM"))
    p <- reader_forest(d, paste0("Primary random-effects synthesis: ", ifelse(is_bray, "Bray-Curtis", "CZM-Aitchison")),
       expression("Hedges "*g[z]*" (95% CI)"), c(d$project, "Random-effects model"),
       c(d$hedges_gz, pool$pooled_hedges_gz), c(d$gz_ci_low, pool$ci_low), c(d$gz_ci_high, pool$ci_high))
    pd <- p$data; pd$txt <- paste0(fmt(pd$est), " [", fmt(pd$lo), ", ", fmt(pd$hi), "]")
    xr <- range(c(pd$lo, pd$hi, 0)); span <- diff(xr); pd$tx <- xr[2] + .10 * span
    return(p + scale_x_continuous(limits = c(xr[1]-.06*span, xr[2]+.8*span)) +
       geom_text(data = pd, aes(x = tx, label = txt), hjust = 0, size = mm8, family = "Arial") +
       labs(caption = paste0("k = ", pool$k_studies, "; \u03c4\u00b2 = ", fmt(pool$tau2), "; I\u00b2 = ", fmt(pool$I2,1), "%")))
  }
  if (grepl("Supplementary_Figure_S5_", stem)) {
    d <- s5; d$lab <- factor(d$label, levels = rev(d$label))
    return(ggplot(d, aes(estimate, lab)) + geom_vline(xintercept = 0, linetype = 2, colour = grey) +
       geom_segment(aes(x = low, xend = high, yend = lab), linewidth = .6, colour = blue) +
       geom_point(aes(shape = type), colour = blue, size = 2.5) +
       scale_shape_manual(values = c("Natural cohort" = 17, "External primary" = 16, "Pooled" = 15),
                          breaks = c("Natural cohort", "External primary", "Pooled"),
                          labels = c("Natural-history cohort", "Independent validation", "Natural-history pooled")) +
       labs(title = "Fixed secondary family-balance change", x = "Mean early-to-late family-balance change", y = NULL,
            caption = "Horizontal intervals are pointwise 95% CIs.") + base_theme())
  }
  if (grepl("Supplementary_Figure_S6_", stem)) {
    d <- s6 %>% pivot_longer(c(`Day-3`, `Day-7`), names_to = "day", values_to = "value")
    idcol <- setdiff(names(s6),c("Day-3","Day-7","change"))[1]
    d$day <- factor(d$day, levels = c("Day-3","Day-7"), labels = c("Day 3","Day 7"))
    return(ggplot(d, aes(day, value, group = .data[[idcol]])) + geom_line(colour = grey, alpha = .6) +
       geom_point(aes(colour = day), size = 2) + scale_colour_manual(values = c(blue,orange), guide = "none") +
       labs(title = "PRJNA851469 distance to healthy reference", x = NULL,
            y = "Mean Bray-Curtis distance to 13 healthy references",
            caption = "14 paired patients. Lines connect repeated observations; they are not CIs.") + base_theme())
  }
  stop("No display override for ", stem)
}
