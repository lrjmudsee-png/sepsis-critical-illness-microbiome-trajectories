# Sepsis V2 投稿收尾与修复清单

更新时间：2026-09-04

## A. 已修复并验证

- [x] 重建损坏的主表工作簿。新工作簿只保留最终 Table 1（10 队列特征）和 Table 2（核心纵向证据）。
- [x] 重建损坏的补充表工作簿，包含 S1-S11，删除原文件中缺失的 drawing/VML 关系和错误工作表维度。
- [x] 两个新工作簿均可重新导入，工作表范围正确，关键数值为数值类型，未发现公式错误。
- [x] 清理图注 DOCX：删除旧 Figure 3 图注、替换提示和 S2 DRAFT 说明，只保留一个 Figure 3 A-D 最终图注和最终 S1/S2 图注。
- [x] 修复 4 个 R 文件的语法错误：Step 68、80、89A、91B。
- [x] 修复 Step 61 和 Step 63 Python 文件中的无效转义警告。
- [x] 完成代码静态语法验证：182 个 R、5 个 Python、139 个 PowerShell 文件通过。
- [x] 将修复后的 6 个源代码文件同步到 V2 投稿代码快照。
- [x] 补齐复现入口、环境说明、数据访问说明和最终打包脚本。
- [x] 原损坏文件已备份至 `V2-SUBMISSION/99_ORIGINAL_BACKUP_20260904/`，不得用于投稿。

## B. 必须由作者补充，当前不能自动完成

- [ ] 确定最终作者顺序。
- [ ] 指定通讯作者并确认邮箱、电话和通讯地址。
- [ ] 补齐每位作者 ORCID。
- [ ] 完成 CRediT author-contribution statement。
- [ ] 填写 Funding、Acknowledgements 和 Competing interests。
- [ ] 确认 Ethics statement/public-data secondary-analysis statement 的正式措辞。
- [ ] 决定目标期刊，并确认是否需要单独 Title Page、Highlights、Graphical Abstract、Reporting Checklist 和推荐审稿人。
- [ ] 填完 Cover Letter 中的期刊名、编辑称谓、通讯作者和声明占位符。
- [ ] 将最终作者区和 Declarations 插入正文 DOCX。
- [ ] 建立公开代码/处理数据仓库，生成正式 URL、版本号和 DOI。
- [ ] 将 DOI/URL 写入正文、Cover Letter 和 Data/Code Availability statement。

## C. 建议在投稿前补强的方法报告

- [ ] 补充 26 个候选项目筛选为 10 个分析队列的检索日期、数据库、关键词、纳排标准和排除原因。
- [ ] 在 Supplementary Methods 中列出各队列 DADA2 的准确过滤/截短/maxEE/拼接参数。
- [ ] 说明 CLR 零值替换或伪计数、丰度/流行率过滤、分类聚合和未分类特征处理。
- [ ] 补充 PERMANOVA 置换次数与限制、随机种子、缺失值处理、dz 定义及模型诊断。
- [ ] 明确 3 个 progressive cohorts 用于跨队列分类学一致性分析的选择规则；若并非预注册，应标为探索性分析。
- [ ] 对参考文献 DOI、期刊名、卷页和 2026 年新文献做一次 PubMed/Crossref 最终核对。

## D. 仍存在的复现性限制

- [ ] 325 个历史脚本仍依赖硬编码的 `E:/sepsis_project` 项目路径；完整公共发布前应参数化并复测。
- [ ] 67 个脚本含自动安装依赖逻辑；建议用 `renv.lock` 或容器替代。
- [ ] 96F2-96F11 图形美化迭代的独立源脚本未恢复，无法完全重现最终版式生成过程。
- [ ] 尚未在全新机器上从 253 GB 原始数据端到端重跑 60-97 全链条；当前验证覆盖静态语法、冻结输入后的终端链和最终数值一致性。
- [ ] 历史 `results/` 中保留了失败或被后续 FIXED 步骤取代的目录；这些是溯源记录，不应打包给期刊。

## E. 投稿前最终放行检查

- [ ] 在 Microsoft Word 中打开正文、图注和 Cover Letter，检查字体、分页、行号、页码及特殊字符。
- [ ] 在 Microsoft Excel 中打开两个最终工作簿，确认无“文件已修复”提示。
- [ ] 检查正文、表格、图注和摘要中的全部样本量、效应量、p 值和 FDR 一致。
- [ ] 按目标期刊要求确认字数、摘要结构、主图/主表数量、TIFF 色彩模式和文件命名。
- [ ] 运行 `V2_97A_export_submission_bundle.ps1` 生成最终 ZIP 和 SHA256 清单。
- [ ] 只上传 `01_投稿用_最终文件` 内的现行文件，不上传 `99_ORIGINAL_BACKUP_20260904` 或历史失败结果目录。

