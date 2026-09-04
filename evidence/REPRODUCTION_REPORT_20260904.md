# Sepsis V2 末端复现报告

## 结论

本次复现通过。基于冻结输入，11 个末端分析与投稿整合步骤全部正常退出，最终生成结果与 `E:/sepsis_project/results` 中的冻结结果在科学内容上保持一致。

- 正式运行时间：2026-09-04 12:19:00 至 12:51:14
- 11 个步骤退出码：全部为 0
- 正式运行累计耗时：1933.828 秒（约 32.2 分钟）
- 新生成结果目录：11 个
- 新生成文件：192 个
- 原始项目结果：未覆盖；本次结果写入独立运行目录

## 结果核对

共核对 192 个生成文件：

| 文件类型 | 数量 | 核对结果 |
|---|---:|---|
| CSV | 132 | 112 个原始内容完全一致；20 个仅运行根目录不同，根目录归一化后全部一致 |
| PDF | 10 | 全部重新渲染后图像 SHA256 一致 |
| PNG | 8 | 全部 SHA256 一致 |
| TXT | 31 | 路径与时间戳归一化后全部一致 |
| 完成标记 `.ok` | 11 | 路径与时间戳归一化后全部一致 |

不存在未解释的数值差异、表格维度差异、图像差异或缺失文件。

核对明细：

- `normalized_reproduction_comparison.csv`：CSV 与 PDF 的主核对表
- `additional_file_comparison.csv`：PNG、TXT、`.ok` 核对表
- `pdf_render_comparison.csv`：PDF 渲染图像核对表
- `shadow_vs_frozen.csv`：未做路径归一化的原始比较表

## 本轮发现并修复的问题

### 1. Windows R 的 UTF-8 环境变量不兼容

当前会话向 Windows 版 R 传入了无效的 `C.UTF-8` 区域设置，使 Step93X4 中的范围连接符“–”被错误读为字母 `b`。数值没有改变，但输出范围字符串损坏。

修复：运行器会在调用 R 前临时移除 `LANG`、`LC_ALL`、`LC_CTYPE`，并显式使用 `--encoding=UTF-8`。随后重跑 Step93X4 及所有下游投稿整合步骤，损坏字符串已恢复，比较通过。

### 2. Step95A2 候选扫描会被下游投稿目录反向污染

原脚本递归扫描整个 `results` 树。随着后续稿件包与复现结果增加，已复制的 Step93U 文本会再次进入候选集合，导致候选表行数随运行时点变化。

修复：将 Step93U 定位器限制在 Step95A2 定义时的上游结果范围，排除 Step95A2 自身及后续稿件、复现和排版目录。修复后候选目录表与原冻结结果维度一致，且所选权威来源未改变。

## 运行步骤

1. `93U_trajectory_ecology_robustness.R`
2. `93W_longitudinal_evidence_freeze.R`
3. `93X4_external_validation_evidence_freeze.R`
4. `94B2D_infection_source_evidence_freeze.R`
5. `V2_36I3_cross_cohort_taxonomic_reproducibility_FINAL.R`
6. `95A2_global_evidence_synthesis_source_anchor_fix.R`
7. `95B_manuscript_source_pack_assembly.R`
8. `95C_manuscript_architecture_figure1_and_main_tables.R`
9. `95E2_table1_analysis_population_fix_and_results_v2.R`
10. `95F3_results_and_table1_final_freeze_sentence_audit.R`
11. `95G_supplementary_tables_and_discussion_v1.R`

逐步退出码和日志路径见 `reproduction_step_summary.csv`；编码修复后的重跑记录见 `corrective_rerun_step_summary.csv`。

## 非阻断警告

R 运行时为 4.4.0，部分已安装包由 R 4.4.3 构建。所有步骤均正常完成，且最终文件比较全部通过，因此本轮未观察到由此导致的结果偏差。正式公开代码时仍建议使用 `renv` 锁定 R 与包版本。

## 范围限制

本次验证的是“冻结输入到投稿末端结果”的复现链，不是从约 253 GB 原始测序数据重新执行下载、DADA2/VSEARCH、去嵌合体和分类注释的全流程。运行目录中的 `data`、`metadata`、`processed_data` 以及 189 个上游结果目录是指向原项目的 Windows 目录连接；因此本次运行记录可在当前电脑复核，但不能脱离 `E:/sepsis_project` 单独搬到另一台电脑直接运行。
