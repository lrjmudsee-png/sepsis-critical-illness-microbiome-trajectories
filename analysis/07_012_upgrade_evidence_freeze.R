# ============================================================
# Sepsis V2 Upgrade - Step 98F
# Evidence freeze + integrated summary (answers the 16 questions)
# ============================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
})

ROOT <- "E:/sepsis_project"
UP <- file.path(ROOT, "results", "V2_UPGRADE_20260907")
OUT <- file.path(UP, "98F_EVIDENCE_FREEZE")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

readcsv <- function(p) suppressMessages(read_csv(p, show_col_types = FALSE))

# ---- collect key numbers ----
dir98A <- file.path(UP, "98A_AITCHISON_LONGITUDINAL_ROBUSTNESS")
dir98B <- file.path(UP, "98B_RANDOM_EFFECTS_META_ANALYSIS")
dir98C <- file.path(UP, "98C_PRJNA1125274_METADATA_FREEZE")
dir98E <- file.path(UP, "98E_PRJNA1125274_EXTERNAL_VALIDATION")

dir_reg <- readcsv(file.path(dir98A, "05_QC_AND_CORRELATION",
                             "V2_98A_bray_aitchison_direction_consistency_registry.csv"))
pooled <- readcsv(file.path(dir98B, "02_POOLED_RESULTS",
                            "V2_98B_PRIMARY_pooled_summary.csv"))
concl <- readcsv(file.path(dir98E,
                           "PRJNA1125274_external_validation_conclusion.csv"))
prim  <- readcsv(file.path(dir98E,
                           "PRJNA1125274_external_validation_primary_contrast.csv"))
acct  <- readcsv(file.path(dir98E,
                           "PRJNA1125274_external_validation_accounting.csv"))
beta  <- readcsv(file.path(dir98E,
                           "PRJNA1125274_T0_T1_community_change_permutation.csv"))
qc98c <- readcsv(file.path(dir98C,
                           "PRJNA1125274_metadata_and_run_QC_summary.csv"))

bray_row <- pooled |> filter(scheme == "bray")
ait_row  <- pooled |> filter(scheme == "aitchison_CZM")
bc_conc  <- concl |> filter(metric == "bray_from_T0")
ac_conc  <- concl |> filter(metric == "aitchison_from_T0_CZM")

sprintf1 <- function(x, d = 2) formatC(as.numeric(x), digits = d, format = "f")

# ---- 16 answers ----
answers <- list()

answers[[1]] <- paste0(
  "1. 哪些原 V2 cohort 成功完成 Aitchison 分析？\n",
  "全部 7 个纵向 cohort（PRJEB82425、PRJNA1166732、PRJNA430161、PRJNA516701、",
  "PRJNA578267、PRJNA691455、PRJNA851469）；PRJNA978257 为静态支持队列（无 ≥2 时点）不适用。")

answers[[2]] <- paste0(
  "2. Bray 与 Aitchison 的方向一致率？\n",
  "配对对比方向一致 6/7 cohort；唯一不一致为非脓毒症外科对照 PRJNA578267",
  "（两种指标均为非显著、接近零：Bray DECREASE p=0.083，Aitchison 基本无方向）。",
  "3 个 primary 自然病程队列（691455/851469/516701）Bray 与 Aitchison 方向 3/3 一致（INCREASE）。")

answers[[3]] <- paste0(
  "3. 哪些 cohort 出现不一致？\n",
  "PRJNA578267（nonsepsis 对照，两者均非显著）；注意 PRJNA691455 的 Aitchison 配对效应为阳性但",
  "不显著（n=9，CI 宽；Bray 显著 p=0.021）——属精度问题而非方向不一致，meta 中已反映。")

answers[[4]] <- paste0(
  "4. primary Meta-analysis 纳入哪些 cohort？为什么？\n",
  "PRJNA691455（CORE_SEPSIS_LONGITUDINAL）、PRJNA851469 与 PRJNA516701",
  "（ICU_BACKGROUND_LONGITUDINAL）：自然病程、有严格/基本共同 baseline anchor、early-to-late estimand 可比。",
  "依据 META_ELIGIBILITY_FREEZE.csv（设计驱动，先于查看新显著性冻结）。PRJEB82425 仅入敏感性池；",
  "PRJNA578267（非脓毒症对照）与干预队列 1166732/430161 单独报告、不并入自然病程效应。")

answers[[5]] <- paste0(
  "5. Bray pooled effect？\n",
  "Hedges g_z = ", sprintf1(bray_row$pooled_hedges_gz), "  (HKSJ 95% CI ",
  sprintf1(bray_row$ci_low), "–", sprintf1(bray_row$ci_high),
  "; tau^2=0, I^2=", sprintf1(bray_row$I2, 1), "%, Q p=", sprintf1(bray_row$Q_p, 3), ").")

answers[[6]] <- paste0(
  "6. Aitchison pooled effect？\n",
  "CZM: Hedges g_z = ", sprintf1(ait_row$pooled_hedges_gz), "  (HKSJ 95% CI ",
  sprintf1(ait_row$ci_low), "–", sprintf1(ait_row$ci_high), "); 伪计数 0.5/1 为 ",
  sprintf1(pooled$pooled_hedges_gz[pooled$scheme == "aitchison_PC0.5"]), " / ",
  sprintf1(pooled$pooled_hedges_gz[pooled$scheme == "aitchison_PC1"]),
  "（方向一致；CI 因异质性较宽）。")

answers[[7]] <- paste0(
  "7. Heterogeneity？\n",
  "Bray: tau^2=0, I^2=0%（Q p=", sprintf1(bray_row$Q_p, 3), "）。",
  "Aitchison(CZM): tau^2=", sprintf1(ait_row$tau2, 3), ", I^2=",
  sprintf1(ait_row$I2, 1), "%（Cochran Q p=", sprintf1(ait_row$Q_p, 3),
  "）——主要由 PRJNA691455 的小效应 vs 851469 的大效应驱动。")

answers[[8]] <- paste0(
  "8. leave-one-out 后是否稳定？\n",
  "是。Bray 移除任一共队列后 pooled g_z 仍在 0.63–0.87（全部为正）；Aitchison(CZM) 在 0.49–1.10",
  "（全部为正；移除 PRJNA851469 时效应最低，提示其对总体量级贡献最大）。方向结论不变。")

answers[[9]] <- paste0(
  "9. patient-level bootstrap 后是否稳定？\n",
  "是。1000 次患者块自举：Bray median g_z 0.79（95%CI 0.58–1.06），Aitchison(CZM) 0.85（0.49–1.26），",
  "两种指标 pooled 效应 100% 为正。注意：k=3 的 HKSJ 在自举内显著性频率偏低（Bray 57%、Aitchison 3%）",
  "反映的是小研究数+异质性下的保守推断，方向与量级稳定。")

answers[[10]] <- paste0(
  "10. PRJNA1125274 最终确认多少患者具有 T0/T1？\n",
  "字母结构层面 50 例（T0 即字母 A、T1 即字母 B）；通过深度≥2000 且两时点均存在、用于",
  "T0→T1 配对分析的为 ", beta$n_patients[1], " 例；T0 样本共 ", acct$runs_T0, " 个 run。")

answers[[11]] <- paste0(
  "11. 多少患者具有完整 T0/T1/T2？\n",
  "字母结构层面 37 例；因 1 例 T0 样本数据完整性剔除（SRR29453361/SAL_78）及深度/样本缺失，",
  "实际可进行三时点配对主终点分析的为 ", prim$n_patients[1], " 例。")

pb <- prim |> filter(metric == "bray_from_T0", sensitivity_set == "ALL_COMPLETE_T0T1T2")
answers[[12]] <- paste0(
  "12. SURVEIL 的 Bray trajectory 与原 V2 一致？\n",
  "是。Δ=D(T0,T2)−D(T0,T1) 配对对比：均值 +", sprintf1(pb$mean_paired_diff),
  "，Hedges g_z=", sprintf1(bc_conc$hedges_gz),
  "（95%CI ", sprintf1(bc_conc$gz_ci_low), "–", sprintf1(bc_conc$gz_ci_high),
  "；p_t=", sprintf1(bc_conc$paired_t_p, 2), "），方向 INCREASE，与 discovery 一致。")

answers[[13]] <- paste0(
  "13. SURVEIL 的 Aitchison trajectory 是否一致？\n",
  "是。CZM: Hedges g_z=", sprintf1(ac_conc$hedges_gz), "（95%CI ",
  sprintf1(ac_conc$gz_ci_low), "–", sprintf1(ac_conc$gz_ci_high),
  "；p_t=", sprintf1(ac_conc$paired_t_p, 3), "），方向 INCREASE；",
  "伪计数 0.5/1 敏感性同样显著（g_z 0.64–0.65）。")

answers[[14]] <- paste0(
  "14. 外部验证最终分级？\n",
  "两种指标均为 Situation A → independent external replication **supported**",
  "（Bray 与 Aitchison 均与 discovery 方向一致且 CI 不含 0）。")

answers[[15]] <- paste0(
  "15. 是否发现会改变原论文核心结论的问题？\n",
  "未发现反向证据；新证据整体支持原核心结论。需如实披露的局限/备注：",
  "(1) SURVEIL 的 run↔时点映射为文档化推断（字母均为 A→B→C 子序列佐证，12/118 例 BioSample 日期异常），",
  "论文未提供逐样本映射表；(2) 1 个 run（SAL_78 T0）因反复下载损坏按数据完整性剔除；",
  "(3) 该队列采用 R1 单端 16S 分析（2×250 无法可靠拼接），区域分辨率与 discovery 队列不同但估量一致；",
  "(4) discovery meta 中 Aitchison I^2≈75%，主要提示队列间效应量异质，方向稳健。")

answers[[16]] <- paste0(
  "16. 建议哪些 Figure/Table 升级？\n",
  "- 主图：新增/并列 Aitchison 几何下的纵向 displacement 图（与 Bray 同框架）；",
  "新增 98B forest 图（Bray 与 Aitchison 分列、标 leave-one-out 区间）；",
  "新增 PRJNA1125274 外部验证图（T0→T1→T2 患者轨迹 + 主对比 Δ 森林/配对图）。",
  "- 表：cohort 级证据表扩展 dz/Hedges g_z（两 metric×零值方案）；pooled 表；外部验证主终点表。",
  "- 补充材料：零处理与流行度预过滤敏感性、日期异常/医院分层、SE-R1 方案说明、SAP。")

writeLines(paste(answers, collapse = "\n\n"),
           file.path(OUT, "V2_UPGRADE_16_questions_summary.md"))

# ---- compact numeric freeze csv ----
freeze <- bind_rows(
  pooled |> select(scheme, pooled_hedges_gz, ci_low, ci_high, tau2, I2, Q_p) |>
    mutate(scope = "PRIMARY_META"),
  concl |> select(metric, hedges_gz, gz_ci_low, gz_ci_high, paired_t_p, direction) |>
    mutate(scope = "PRJNA1125274_EXTERNAL_VALIDATION")
)
write_csv(freeze, file.path(OUT, "V2_UPGRADE_key_effect_freeze.csv"))

cat("98F done. 16-question summary and effect freeze written to:\n", OUT, "\n")
