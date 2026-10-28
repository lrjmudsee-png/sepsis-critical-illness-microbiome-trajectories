# Display layer only: frozen estimates, selection and inference are unchanged.
# Explanatory sentences remain in the manuscript/supplementary legends.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

base_theme <- function() {
  theme_classic(base_size = 9, base_family = "Arial") +
    theme(plot.title = element_text(size = 10, face = "bold", margin = margin(b = 4)),
          plot.subtitle = element_text(size = 8, margin = margin(b = 5), colour = "#4C4C4C"),
          axis.title = element_text(size = 8.5),
          axis.title.x = element_text(margin = margin(t = 5)),
          axis.title.y = element_text(margin = margin(r = 5)),
          axis.text = element_text(size = 8, colour = "#222222"),
          axis.line = element_line(linewidth = .35, colour = "#333333"),
          axis.ticks = element_line(linewidth = .3),
          axis.ticks.length = grid::unit(1.1, "mm"),
          legend.text = element_text(size = 8), legend.title = element_blank(),
          legend.key.height = grid::unit(3.5, "mm"),
          legend.key.width = grid::unit(4, "mm"),
          legend.spacing.x = grid::unit(2, "mm"),
          legend.margin = margin(0,0,0,0), legend.box.margin = margin(1,0,0,0),
          plot.margin = margin(5, 7, 5, 5), legend.position = "bottom",
          plot.tag = element_text(face = "bold", size = 11),
          plot.tag.position = "topleft")
}
theme_set(base_theme())
mm8 <- 8 / ggplot2::.pt
panel <- function(p, tag) p + labs(tag = tag) +
  theme(plot.tag.position = "topleft", plot.tag = element_text(face="bold", size=11))
term_pair <- function(x) gsub(" vs ", "\nvs ", x, fixed = TRUE)

reader_forest <- function(dd, title, axis_title, labels, est, lo, hi, ptxt = NULL, col = blue) {
  d <- data.frame(label = labels, est = est, lo = lo, hi = hi)
  d$label <- factor(d$label, levels = rev(unique(labels)))
  xr <- range(c(lo, hi, 0), na.rm = TRUE); span <- diff(xr)
  if (span == 0) span <- 1
  p <- ggplot(d, aes(y=label, x=est)) +
    geom_vline(xintercept=0, linetype=2, colour="#A0A0A0", linewidth=.35) +
    geom_segment(aes(x=lo,xend=hi,yend=label),colour=col,linewidth=.65) +
    geom_point(colour=col,size=2.3) +
    scale_x_continuous(limits=c(xr[1]-.055*span,xr[2]+ifelse(is.null(ptxt),.07,.70)*span)) +
    scale_y_discrete(expand=expansion(add=.5)) +
    labs(title=title,x=axis_title,y=NULL) + base_theme()
  if (!is.null(ptxt)) {
    d$ptxt <- ptxt; d$text_x <- xr[2]+.07*span
    p <- p + geom_text(data=d,aes(x=text_x,label=ptxt),hjust=0,size=mm8,family="Arial")
  }
  p
}

corrected_plot <- function(stem) {
  if (stem == "Figure_1_Study_Architecture") {
    branches <- c("DISCOVERY","INDEPENDENT_LONGITUDINAL_VALIDATION",
                  "EXTERNAL_CROSS_SECTIONAL_COMPARISON","INFECTION_SOURCE_EXTENSION")
    n_total <- sum(f1$cohort_count[match(branches,f1$branch)])
    nd <- f1$cohort_count[f1$branch=="DISCOVERY_LONGITUDINAL"]
    ns <- f1$cohort_count[f1$branch=="DISCOVERY_STATIC"]
    centres <- c(1.5,4.5,7.5,10.5)
    boxes <- data.frame(x0=centres-1.36,x1=centres+1.36,
                        col=c(blue,"#8C959C","#8C959C","#8C959C"),
                        fill=rep("#F7F8F9",4))
    p <- ggplot() + coord_cartesian(xlim=c(0,12),ylim=c(0,5),expand=FALSE) +
      theme_void(base_family="Arial") + theme(plot.margin=margin(5,4,4,4)) +
      annotate("rect",xmin=3.3,xmax=8.7,ymin=4.04,ymax=4.86,fill="white",colour="#6B6B6B",linewidth=.45) +
      annotate("text",x=6,y=4.56,label="Public gut-microbiome cohorts",fontface="bold",size=3.5,family="Arial") +
      annotate("text",x=6,y=4.22,label=paste0("n = ",n_total," cohorts"),size=2.95,family="Arial") +
      annotate("segment",x=6,xend=6,y=4.04,yend=3.78,colour=grey,linewidth=.4) +
      annotate("segment",x=min(centres),xend=max(centres),y=3.78,yend=3.78,colour=grey,linewidth=.4) +
      geom_segment(data=data.frame(x=centres),aes(x=x,xend=x),y=3.78,yend=3.48,
                   colour=grey,linewidth=.4,arrow=grid::arrow(length=grid::unit(1.4,"mm"))) +
      geom_rect(data=boxes,aes(xmin=x0,xmax=x1),ymin=1.08,ymax=3.47,
                fill=boxes$fill,colour=boxes$col,linewidth=.5)
    titles <- c("Discovery","Independent\nlongitudinal validation",
                "Cross-sectional\ncomparison","Infection-source\nextension")
    details <- c(paste0(nd," longitudinal cohorts\n",ns," static supportive cohort\nPRJNA978257 (static)"),
                 "PRJNA1125274\nWithin-patient change",
                 "PRJNA1010969\nControl / Trauma / Sepsis\nEcological state",
                 "CRA002354\nPulmonary vs recorded\nnon-pulmonary")
    p <- p + annotate("text",x=centres,y=2.96,label=titles,fontface="bold",size=2.95,
                       lineheight=1.12,family="Arial") +
      annotate("text",x=centres,y=1.92,label=details,size=2.78,lineheight=1.3,family="Arial") +
      geom_segment(data=data.frame(x=centres),aes(x=x,xend=x),y=1.08,yend=.77,colour=grey,linewidth=.4) +
      annotate("segment",x=min(centres),xend=max(centres),y=.77,yend=.77,colour=grey,linewidth=.4) +
      annotate("segment",x=6,xend=6,y=.77,yend=.54,colour=grey,linewidth=.4,
               arrow=grid::arrow(length=grid::unit(1.4,"mm"))) +
      annotate("text",x=6,y=.26,label="Within-cohort analyses and cross-cohort synthesis",
               size=3.05,family="Arial")
    return(p)
  }

  if (grepl("^Figure_2",stem)) {
    d <- f2
    d$lab <- factor(paste0(d$project,"  (n = ",d$primary_n,")"),
                    levels=rev(paste0(d$project,"  (n = ",d$primary_n,")")))
    return(ggplot(d,aes(y=lab)) +
      geom_vline(xintercept=0,linetype=2,colour="#999999",linewidth=.35) +
      geom_segment(aes(x=primary_effect_dz,xend=common_anchor_effect_dz,yend=lab),
                   linetype=3,colour="#888888",linewidth=.5) +
      geom_point(aes(x=primary_effect_dz,colour="Primary early-to-late",shape="Primary early-to-late"),size=2.7) +
      geom_point(aes(x=common_anchor_effect_dz,colour="Common-anchor sensitivity",shape="Common-anchor sensitivity"),
                 fill="white",size=2.8,stroke=.7) +
      geom_text(data=d[d$primary_significant,],aes(x=primary_effect_dz,label="*"),nudge_y=.45,colour=blue,size=3.6,family="Arial") +
      scale_colour_manual(name="Analysis",breaks=c("Primary early-to-late","Common-anchor sensitivity"),
                          values=c("Primary early-to-late"=blue,"Common-anchor sensitivity"=orange)) +
      scale_shape_manual(name="Analysis",breaks=c("Primary early-to-late","Common-anchor sensitivity"),
                         values=c("Primary early-to-late"=16,"Common-anchor sensitivity"=21)) +
      scale_y_discrete(expand=expansion(add=.75)) +
      labs(x=expression("Standardized change in Bray-Curtis displacement ("*d[z]*")"),y=NULL) + base_theme())
  }

  if (grepl("^Figure_3",stem)) {
    da <- f3a %>% pivot_longer(c(Day3,Day7),names_to="day",values_to="value")
    da$day <- factor(da$day,levels=c("Day3","Day7"),labels=c("Day 3","Day 7"))
    a <- ggplot(da,aes(day,value,group=Patient_true)) +
      geom_line(colour="#919191",alpha=.65,linewidth=.4) +
      geom_point(aes(colour=day),size=2) + scale_colour_manual(values=c(blue,orange),guide="none") +
      labs(title="Within-patient displacement",subtitle="PRJNA691455; 9 paired patients",
           x=NULL,y="Bray-Curtis displacement") + base_theme()
    b <- ggplot(f3b,aes(x="",y=signed_Euclidean)) +
      geom_boxplot(width=.34,outlier.shape=NA,fill="#E8F2F8",colour=blue,linewidth=.5) +
      geom_jitter(colour=grey,alpha=.65,size=1.3,
                  position=position_jitter(width=.08,height=0,seed=100)) + coord_flip() +
      labs(title="Taxonomic trajectory heterogeneity",subtitle="45 comparisons; 10 patients",
           x=NULL,y="Signed-trajectory Euclidean distance") + base_theme()
    dc <- f3c
    dc$label <- c(latest_EII="Exploratory EII",latest_Simpson_instability="Simpson-diversity\ninstability",
                  latest_Bray="Latest Bray-Curtis\ndisplacement")[dc$outcome]
    dc$label <- factor(dc$label,levels=rev(c("Exploratory EII","Simpson-diversity\ninstability","Latest Bray-Curtis\ndisplacement")))
    c <- ggplot(dc,aes(spearman_rho,label)) +
      geom_vline(xintercept=0,linetype=2,colour="#999999",linewidth=.35) + geom_point(colour=blue,size=2.2) +
      geom_text(aes(x=1.04,label=paste0("p = ",vapply(permutation_p_10000,fmt_p,character(1)))),
                hjust=0,size=mm8,family="Arial") +
      scale_x_continuous(limits=c(-.05,1.62),breaks=c(0,.5,1)) +
      labs(title="Heterogeneity and\necological measures",subtitle="Patient-level correlations; n = 10",
           x=expression("Spearman "*rho),y=NULL) + base_theme()
    dd <- f3d
    dd$lab <- factor(term_pair(dd$pair_label),levels=rev(unique(term_pair(dd$pair_label))))
    d <- ggplot(dd,aes(spearman_rho,lab,colour=rank,shape=rank)) +
      geom_vline(xintercept=0,linetype=2,colour="#999999",linewidth=.35) +
      geom_point(position=position_dodge(width=.4),size=2.4) +
      scale_colour_manual(values=c(Genus=blue,Family=orange)) +
      scale_shape_manual(values=c(Genus=16,Family=17)) +
      scale_x_continuous(limits=c(-.05,.55),breaks=c(0,.1,.2,.3,.4,.5)) +
      labs(title="Cross-cohort taxonomic concordance",x=expression("Spearman "*rho),y=NULL) + base_theme()
    return(wrap_plots(A=panel(a,"a"),B=panel(b,"b"),C=panel(c,"c"),D=panel(d,"d"),
                      design="AB\nAC\nDD",widths=c(1,1.25),heights=c(.83,1.25,1.15)))
  }

  if (grepl("^Figure_4",stem)) {
    da <- read_pack("FIGURE4a_sample_distances.csv") %>% filter(Group %in% c("Trauma","Sepsis"))
    da$Group <- factor(da$Group,levels=c("Trauma","Sepsis"))
    med <- data.frame(Group=factor(c("Trauma","Sepsis"),levels=c("Trauma","Sepsis")),
                      value=c(f4a$median_B[1],f4a$median_A[1]))
    a <- ggplot(da,aes(Group,Bray_to_Control_Centroid,colour=Group)) +
      geom_boxplot(width=.45,outlier.shape=NA,fill="white",linewidth=.5) +
      geom_jitter(size=1.7,alpha=.8,position=position_jitter(width=.10,height=0,seed=100)) +
      geom_text(data=med,aes(y=1.03,label=paste0("Median ",fmt(value))),size=mm8,family="Arial") +
      scale_colour_manual(values=c(Trauma=grey,Sepsis=orange),guide="none") +
      scale_x_discrete(labels=c(Trauma="Trauma\nn = 18",Sepsis="Sepsis\nn = 18")) +
      scale_y_continuous(limits=c(0,1.08),breaks=seq(0,1,.25)) +
      labs(title="Ecological-state differentiation",subtitle="PRJNA1010969; control reference n = 17",
           x=NULL,y="Bray-Curtis distance from control centroid") + base_theme()
    db <- f4b; db$pc <- factor(format(db$pseudocount,scientific=TRUE,digits=1),
                             levels=format(db$pseudocount,scientific=TRUE,digits=1))
    b <- ggplot(db,aes(pc,PERMANOVA_R2)) + geom_point(colour=blue,size=2.5) +
      geom_text(aes(label=fmt(PERMANOVA_R2)),nudge_y=.0018,size=mm8,family="Arial") +
      scale_y_continuous(limits=c(.077,.089),breaks=c(.078,.082,.086)) +
      labs(title="Aitchison pseudocount sensitivity",subtitle="Sepsis vs trauma; n = 36 samples",
           x="Pseudocount",y=expression("PERMANOVA "*R^2)) + base_theme()
    dc <- f4c
    dc$label <- factor(c("Sepsis-control","Sepsis-trauma"),levels=rev(c("Sepsis-control","Sepsis-trauma")))
    c <- ggplot(dc,aes(spearman_rho,label)) +
      geom_vline(xintercept=0,linetype=2,colour="#999999",linewidth=.35) + geom_point(colour=green,size=2.5) +
      geom_text(aes(x=.43,label=paste0("n = ",overlapping_genera,"\np = ",vapply(p_value,fmt_p,character(1)))),
                hjust=0,size=mm8,lineheight=1.12,family="Arial") +
      scale_x_continuous(limits=c(-.04,.72),breaks=c(0,.1,.2,.3,.4)) +
      labs(title="Exploratory taxonomic alignment",subtitle="Longitudinal vs cross-sectional contrasts",
           x=expression("Spearman "*rho),y=NULL) + base_theme()
    return(wrap_plots(A=panel(a,"a"),B=panel(b,"b"),C=panel(c,"c"),
                      design="AB\nAC",widths=c(1,1.05),heights=c(1,1)))
  }

  if (grepl("^Figure_5",stem)) {
    a <- reader_forest(f5a,"Primary Bray-Curtis\ninteraction",expression("Source-by-time interaction "*beta),
                      "Pulmonary vs recorded\nnon-pulmonary",f5a$effect,f5a$ci_low,f5a$ci_high,
                      paste0("p = ",vapply(f5a$p_value,fmt_p,character(1))),orange)
    b <- reader_forest(f5b,"Source-specific slopes","Bray-Curtis displacement\nslope per day",
                      paste0(gsub("Recorded non-pulmonary","Recorded\nnon-pulmonary",f5b$group_label),
                             "\n",f5b$npatients," patients; ",f5b$nobs," obs."),
                      f5b$slope_per_day,f5b$ci95_low,f5b$ci95_high,
                      paste0("p = ",vapply(f5b$p_LRT,fmt_p,character(1))))
    smap <- c(ADJUSTED_SOFA="SOFA adjustment",ADJUSTED_APACHE="APACHE adjustment",
              ADJUSTED_LACTATE="Lactate adjustment",EXCLUDE_OTHER_UNKNOWN="Exclude other/unknown",
              RAREFIED4000_BRAY="Rarefied 4,000",SENSITIVITY_MIN2000="Minimum depth 2,000")
    sl <- unname(smap[f5c$sensitivity]);stopifnot(!anyNA(sl))
    c <- reader_forest(f5c,"Sensitivity analyses",expression("Source-by-time interaction "*beta),sl,
                      f5c$beta,f5c$ci_low,f5c$ci_high,paste0("p = ",vapply(f5c$p_value,fmt_p,character(1))),green)
    dd <- data.frame(label=c("Unrestricted","ICU Day 3 or earlier"),
                     r2=c(f5d$primary_PERMANOVA_R2[1],f5d$day3_restricted_PERMANOVA_R2[1]),
                     p=c(f5d$primary_PERMANOVA_p[1],f5d$day3_restricted_PERMANOVA_p[1]),
                     disp=c(f5d$primary_dispersion_p[1],f5d$day3_restricted_dispersion_p[1]))
    dd$label <- factor(dd$label,levels=dd$label)
    dd$label <- factor(dd$label,levels=rev(levels(dd$label)))
    d <- ggplot(dd,aes(r2,label)) + geom_point(colour=blue,size=2.5) +
      geom_text(aes(x=.055,label=paste0("PERMANOVA p = ",vapply(p,fmt_p,character(1)),
                                      "; dispersion p = ",vapply(disp,fmt_p,character(1)))),
                hjust=0,size=mm8,family="Arial") +
      scale_x_continuous(limits=c(0,.13),breaks=c(0,.01,.02,.03,.04)) +
      scale_y_discrete(expand=expansion(add=.55)) +
      labs(title="Baseline beta-diversity robustness",x=expression("PERMANOVA "*R^2),y=NULL) + base_theme()
    return(wrap_plots(A=panel(a,"a"),B=panel(b,"b"),C=panel(c,"c"),D=panel(d,"d"),
                      design="AB\nCC\nDD",widths=c(1,1.05),heights=c(1.1,1.4,.72)))
  }

  if (grepl("Supplementary_Figure_S1_",stem)) {
    build <- function(role) {
      d <- s1 %>% filter(analysis_role==role)
      d$lab <- factor(term_pair(d$pair_label),levels=rev(unique(term_pair(d$pair_label))))
      ggplot(d,aes(spearman_rho,lab,colour=rank,shape=rank)) +
        geom_vline(xintercept=0,linetype=2,colour="#999999",linewidth=.35) +
        geom_point(position=position_dodge(.35),size=2.5) +
        scale_colour_manual(values=c(Genus=blue,Family=orange)) +
        scale_shape_manual(values=c(Genus=16,Family=17)) +
        labs(title=role,x=expression("Spearman "*rho*" of longitudinal CLR effect vectors"),y=NULL) + base_theme()
    }
    return((panel(build("Primary"),"a") / panel(build("Common-anchor sensitivity"),"b")) +
             plot_layout(guides="collect") & theme(legend.position="bottom"))
  }
  if (grepl("Supplementary_Figure_S2_",stem)) {
    d <- s2; d$lab <- factor(d$project,levels=rev(d$project))
    return(ggplot(d,aes(y=lab)) +
      geom_segment(aes(x=primary_n,xend=common_anchor_n,yend=lab),colour=light_grey,linewidth=.6) +
      geom_point(aes(primary_n,colour="Primary",shape="Primary"),size=2.6) +
      geom_point(aes(common_anchor_n,colour="Common-anchor sensitivity",shape="Common-anchor sensitivity"),
                 fill="white",size=2.8,stroke=.7) +
      scale_colour_manual(name="Analysis",breaks=c("Primary","Common-anchor sensitivity"),
                          values=c(Primary=blue,"Common-anchor sensitivity"=orange)) +
      scale_shape_manual(name="Analysis",breaks=c("Primary","Common-anchor sensitivity"),
                         values=c(Primary=16,"Common-anchor sensitivity"=21)) +
      scale_y_discrete(expand=expansion(add=.45)) +
      labs(title="Paired patient counts by analysis",x="Number of paired patients",y=NULL) + base_theme())
  }
  if (grepl("Supplementary_Figure_S[34]_",stem)) {
    is_bray <- grepl("S3_",stem);d <- if(is_bray)s3 else s4
    pool <- s34 %>% filter(scheme==ifelse(is_bray,"bray","aitchison_CZM"))
    p <- reader_forest(d,paste0("Random-effects synthesis: ",ifelse(is_bray,"Bray-Curtis","CZM-Aitchison")),
                       expression("Hedges "*g[z]*" (95% CI)"),c(d$project,"Random-effects model"),
                       c(d$hedges_gz,pool$pooled_hedges_gz),c(d$gz_ci_low,pool$ci_low),c(d$gz_ci_high,pool$ci_high))
    pd <- p$data;pd$txt <- paste0(fmt(pd$est)," [",fmt(pd$lo),", ",fmt(pd$hi),"]")
    xr <- range(c(pd$lo,pd$hi,0));span <- diff(xr);pd$tx <- xr[2]+.1*span
    return(p + scale_x_continuous(limits=c(xr[1]-.06*span,xr[2]+.75*span)) +
             geom_text(data=pd,aes(x=tx,label=txt),hjust=0,size=mm8,family="Arial") +
             labs(subtitle=paste0("k = ",pool$k_studies,"; \u03c4\u00b2 = ",fmt(pool$tau2),
                                  "; I\u00b2 = ",fmt(pool$I2,1),"%")))
  }
  if (grepl("Supplementary_Figure_S5_",stem)) {
    d <- s5;d$lab <- factor(d$label,levels=rev(d$label))
    return(ggplot(d,aes(estimate,lab)) +
      geom_vline(xintercept=0,linetype=2,colour="#999999",linewidth=.35) +
      geom_segment(aes(x=low,xend=high,yend=lab),linewidth=.65,colour=blue) +
      geom_point(aes(shape=type),colour=blue,size=2.5) +
      scale_shape_manual(values=c("Natural cohort"=17,"External primary"=16,"Pooled"=15),
                         breaks=c("Natural cohort","External primary","Pooled"),
                         labels=c("Natural-history cohort","Independent validation","Natural-history pooled")) +
      labs(title="Fixed secondary family-balance change",x="Mean early-to-late family-balance change",y=NULL) + base_theme())
  }
  if (grepl("Supplementary_Figure_S6_",stem)) {
    d <- s6 %>% pivot_longer(c(`Day-3`,`Day-7`),names_to="day",values_to="value")
    idcol <- setdiff(names(s6),c("Day-3","Day-7","change"))[1]
    d$day <- factor(d$day,levels=c("Day-3","Day-7"),labels=c("Day 3","Day 7"))
    return(ggplot(d,aes(day,value,group=.data[[idcol]])) +
      geom_line(colour=grey,alpha=.6,linewidth=.4) +
      geom_point(aes(colour=day),size=2) + scale_colour_manual(values=c(blue,orange),guide="none") +
      labs(title="PRJNA851469 distance to healthy reference",subtitle="14 paired patients; 13 healthy references",
           x=NULL,y="Mean Bray-Curtis distance to healthy references") + base_theme())
  }
  stop("No display override for ",stem)
}
