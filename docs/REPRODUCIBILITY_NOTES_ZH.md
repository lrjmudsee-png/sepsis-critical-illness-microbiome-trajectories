# V2 复现代码包 —— 入口说明

> 本包为 V2（Sepsis_V2）分析的可复现性骨架，对标 V1 的 RC1 复现包。
> 原则：原始数据不随包分发（260 GB），通过公开 accession 获取；代码完整；关键数字可核对。
> 状态：2026-09-04 已完成代码静态语法修复与基础复现文件补齐。请以 `00_REPRODUCE_V2_README.md` 和 `V2-SUBMISSION/V2_投稿收尾与修复清单_20260904.md` 为最新入口。

## 1. 包结构

```
03_代码\
├── 00_REPRODUCE_V2_README.md   ← 最新入口
├── 01_RUN_ORDER.csv            ← 脚本 → results 步骤 → 用途 全映射（48 行）
├── 02_results步骤状态表.csv    ← results 全部 200 个步骤目录的状态扫描
├── 03_ENVIRONMENT.md           ← 已核对的软件环境摘要
├── 04_DATA_ACCESS.md           ← 数据 accession 与冻结检查点
└── 05_MANIFEST_SHA256.csv      ← 投稿包 SHA256 清单
```

代码本体在 `harness\V2\代码快照_03_data_processing\`（561 文件：R/py/ps1 + 压缩包 zip 备份）。

## 2. 运行顺序（摘要，详见 01_RUN_ORDER.csv）

1. **STEP60–81（元数据→队列冻结）**：60(ps1)→61(py)→63/64/65(py/R)→66–69(R)→70/70B→71/71B→72→73→74→75→76(B/C/D)→77→78A/A2（78A3/A4/A4B 在 code\02_data_download）→79→80(B–F)→81(B/C/D)
   - 终点：`results\V2_21D_FINAL_FREEZE_CLINICAL_REPAIRED\V2_FINAL_FROZEN_master_metadata_*.csv`（唯一元数据入口）
2. **STEP82–87（序列处理）**：82(B/C)→83A（构建 `data\_V2_ANALYSIS_READY` 硬链接工作区）→83B/C→84A/A2/B(B3/B4A/B4B)/C→85A(A2/A3)/B(B2)/C→86B→87A/A2→87B（分析对象冻结，785 观测）
3. **STEP88–93（统计分析）**：88A/A2→88B→89A/A2→89B→90A→90B→91A→91B(B2/B3/B4)→92A→92B(B2)→92C→92D→93A–93X4（49 个 R）
4. **STEP94（CRA002354 特殊路线）**：94A(A2/A3/A4)→94B1(A–E)→94B2(A–D)→94B3(A–C)
5. **STEP95–97（论文）**：95A(A2)→95B→95C→95D2/D3→95E(E2)→95F(F2/F3)→95G→96B→96C→96D(D2/D3/D4)→96E→96F→96F4–F11（⚠️无脚本，见下）→96G→96H→96I(I2)→96J3→96K(K2/K3)→97A（⚠️无脚本）

运行方式：每步 `RUN_STEP*.ps1`（87B 起有运行器）或直接 Rscript；日志与 `.ok` 完成标记在对应 results 目录。

## 3. 环境要求

- R 4.4.0；包：DADA2 1.34.0、vegan 2.7.3、lme4 2.0.1、lmerTest 3.2.1；VSEARCH 2.31.0（`tools\vsearch-2.31.0-win-x86_64\`）
- 参考库：SILVA 138.2（`SILVA_138_2\silva_nr99_v138.2_toGenus_trainset.fa.gz`、`silva_v138.2_assignSpecies.fa.gz`）
- 下载工具：aria2（ENA manifest 断点续传）
- 已补：`03_ENVIRONMENT.md`（最终包 sessionInfo 与本地统计环境摘要）

## 4. 数据获取（详见 04_DATA_ACCESS.md）

10 个分析队列 accession（与稿件 Table 1 一致）：
- NCBI SRA：PRJNA691455、PRJNA516701、PRJNA851469、PRJNA578267、PRJNA430161、PRJNA1166732、PRJNA978257、PRJNA1010969
- ENA：PRJEB82425
- NGDC-GSA：CRA002354（FASTA 无质量分，走 VSEARCH 97% OTU 路线）

## 5. 期望值核对（复现验证点）

以 `Manuscript_v5_2_SUBMISSION_DRAFT_FIXED2.txt` 摘要为准：
dz=0.765/1.095/0.623（FDR 0.0365/0.0123/0.0365）；0.725→0.904；mean increase 0.113；
cosine similarity=0.162；EII rho=0.733(p=0.023)/Simpson rho=0.673(p=0.038)；
rho=0.36(属)/0.08(科)；36.0%/31.8%；0.860 vs 0.657（FDR=0.00132）；R2 0.081–0.085；
beta=0.017（95% CI −0.024~0.058, p=0.426）
（机器核对：`投稿素材包\05_QC与溯源\01_ABSTRACT_NUMERIC_LOCK_QC.csv` 全 TRUE）

## 6. 已知缺口（投稿前必须处理）

1. **96F2–96F11 图美化迭代仍无独立源脚本**；97A 已补 `V2_97A_export_submission_bundle.ps1` 用于最终文件导出、压缩和 SHA256 清单，但不能反向重建缺失的图形美化源代码。
2. 60–87A 无 RUN_*.ps1 与 README_STEP（仅 87B 起有运行器、92A 起有说明）
3. 编号空洞：62/V2_02、V2_04 无产物；93M/93R_fix 等仅存在于 zip
4. 环境快照待汇总（`03_ENVIRONMENT.md`）
