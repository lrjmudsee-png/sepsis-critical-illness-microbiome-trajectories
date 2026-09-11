"""Generate the requested evidence report from completed metadata/statistical audits."""
import csv
import hashlib
import json
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path
from collections import defaultdict

ROOT=Path('E:/sepsis_project')
OUT=ROOT/'results/V2_UPGRADE_20260907/98G_PRJNA1125274_PATIENT_RECONCILIATION_20260909'
UP=ROOT/'results/V2_UPGRADE_20260907'

def read(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(name,rows):
    rows=list(rows)
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        w.writeheader();w.writerows(rows)
def fmt(x):
    return f'{float(x):.8g}'
def md(rows,cols):
    header='| '+' | '.join(cols)+' |\n| '+' | '.join('---' for _ in cols)+' |\n'
    return header+''.join('| '+' | '.join(str(r.get(c,'')) for c in cols)+' |\n' for r in rows)

patients=read(OUT/'PRJNA1125274_patient_reconciliation_table.csv')
runs=read(OUT/'PRJNA1125274_run_level_reconciliation.csv')
centers=read(OUT/'PRJNA1125274_hospital_patient_sample_summary.csv')
comp=read(OUT/'PRJNA1125274_primary_results_before_after_correction.csv')
senspath=OUT/'PRJNA1125274_orphan_exclusion_sensitivity_comparison.csv'
sens=read(senspath) if senspath.exists() else []
assert sens, '98E orphan exclusion sensitivity must finish before reporting'
obj=read(OUT/'analysis_object_run_accounting.csv')
disp=read(UP/'98E_PRJNA1125274_EXTERNAL_VALIDATION/PRJNA1125274_external_validation_sample_displacement.csv')
primary_ids={r['patient_id'] for r in patients if r['in_current_30_primary']=='True'}
stages=[]
def stage(name,rs):
    group=defaultdict(set)
    for r in rs:group[r['patient_id']].add(r['timepoint'])
    return {'stage':name,'samples':len(rs),'patient_alias_groups':len(group),
            'patients_with_T0':sum('T0' in t for t in group.values()),
            'complete_T0_T1_T2_recomputed':sum(t=={'T0','T1','T2'} for t in group.values())}
stages.append(stage('PUBLIC_289_RUN_METADATA',runs))
stages.append(stage('EXISTING_98D_OBJECT',obj))
stages.append(stage('DEPTH_GE_2000', [r for r in obj if r['passes_depth']=='TRUE']))
stages.append(stage('HAS_DEPTH_PASSING_T0_DISPLACEMENT',disp))
stages.append(stage('CURRENT_PRIMARY_COMPLETE_CASES',[r for r in disp if r['patient_id'] in primary_ids]))
write('PRJNA1125274_reconciled_population_accounting.csv',stages)
losses=[r for r in runs if r['complete_T0_T1_T2']=='TRUE' and r['in_current_30_primary']=='False']
write('PRJNA1125274_complete_case_attrition_details.csv',losses)
source_facts=[
 {'source':'Paper abstract / Results 3.1','claim':'Total patients','reported_value':'132','audit':'Public alias groups 134; baseline-observed IDs 132','status':'DENOMINATOR_UNRESOLVED'},
 {'source':'Paper Results 3.1 and figure legends 1/2/4/5','claim':'Age strata','reported_value':'85 >=65; 49 <65','audit':'85+49=134, not 132. These are age strata, NOT hospital counts.','status':'PUBLICATION_INTERNAL_INCONSISTENCY'},
 {'source':'Paper Results 3.1','claim':'Colonization strata','reported_value':'94 not colonized;38 colonized','audit':'94+38=132; no per-patient clinical crosswalk','status':'AGGREGATE_ONLY'},
 {'source':'Paper figure legends','claim':'Sequence samples','reported_value':'289','audit':'289 unique BioSamples and runs locally and in current ENA','status':'MATCH'},
 {'source':'Paper Methods 2.1','claim':'Recruitment dates and centers','reported_value':'June 2022-August 2023; Alessandria and Novara','audit':'Center-specific n not reported in accessible article','status':'NO_CENTER_DENOMINATOR'},
 {'source':'Paper Methods 2.1','claim':'Exclusion','reported_value':'No signed informed consent','audit':'No numerical screened/excluded flow or named excluded IDs reported','status':'NO_PATIENT_LEVEL_EXCLUSION_LIST'},
 {'source':'Current Europe PMC JATS + Wiley full-text sections','claim':'Supplementary materials','reported_value':'pmc-prop-has-supplement=no; 0 supplementary-material nodes','audit':'No patient enrollment table / crosswalk found; not proof none held privately by authors','status':'NOT_PUBLICLY_AVAILABLE'},
 {'source':'ClinicalTrials.gov NCT07183722 retrieved 2026-09-09','claim':'Registry enrollment','reported_value':'86 ACTUAL; Alessandria location only','audit':'Cannot equate this single-site registry denominator with multicenter 132/134; no Novara crosswalk','status':'NOT_COMPARABLE_FOR_RESOLUTION'},
]
write('PRJNA1125274_publication_denominator_evidence.csv',source_facts)

main=[]
for r in comp:
    if r['sensitivity_set']=='ALL_COMPLETE_T0T1T2' and r['metric'] in ['bray_from_T0','aitchison_from_T0_CZM']:
        main.append({'指标':r['metric'],'完整患者数':r['n_complete_before']+' → '+r['n_complete_after'],
                     '均值Δ（核对前=核对后）':fmt(r['mean_contrast_after']),
                     '均值Δ 95% CI':f"[{fmt(r['mean_contrast_ci_low_after'])}, {fmt(r['mean_contrast_ci_high_after'])}]",
                     'Cohen dz':fmt(r['cohen_dz_after']),
                     'Hedges gz（95% CI）':f"{fmt(r['hedges_gz_after'])} [{fmt(r['gz_ci_low_after'])}, {fmt(r['gz_ci_high_after'])}]",
                     'paired t p':fmt(r['p_value_after'])})
sensitivity=[]
for r in sens:
    if r['sensitivity_set']=='ALL_COMPLETE_T0T1T2' and r['metric'] in ['bray_from_T0','aitchison_from_T0_CZM']:
        sensitivity.append({'指标':r['metric'],'n':r['n_patients_orphan_exclusion'],
                            '均值Δ':fmt(r['mean_paired_diff_orphan_exclusion']),
                            'Cohen dz':fmt(r['cohen_dz_orphan_exclusion']),
                            'Hedges gz（95% CI）':f"{fmt(r['hedges_gz_orphan_exclusion'])} [{fmt(r['gz_ci_low_orphan_exclusion'])}, {fmt(r['gz_ci_high_orphan_exclusion'])}]",
                            'paired t p':fmt(r['paired_t_p_orphan_exclusion'])})
center_table=[{'中心':r['hospital'],'编号组':r['n_patient_alias_groups'],'样本':r['n_samples'],'T0':r['samples_T0'],'T1':r['samples_T1'],'T2':r['samples_T2'],'完整ABC（深度前）':r['complete_before_depth'],'主分析（深度后）':r['complete_primary_after_depth']} for r in centers]
extra_table=[{'patient_id':r['patient_id'],'时间点':r['timepoint'],'alias':r['sample_alias'],'BioSample':r['biosample'],'Run':r['run_id'],'采样日期':r['collection_date'],'nonchim reads':r['nonchim_reads'],'在30人主分析':r['in_current_30_primary']} for r in runs if r['no_T0_alias_group']=='True']
single_ids=[r['patient_id'] for r in patients if r['single_sample_patient']=='True']
date_ids=[r['patient_id'] for r in patients if r['date_order_anomaly']=='True']
lost_ids=sorted({r['patient_id'] for r in losses})
title_conflicts=[r for r in runs if r['title_mismatch_category']=='HOSPITAL_CODE_CONFLICT']
title_table=[{'patient_id':r['patient_id'],'Run':r['run_id'],'正式alias':r['sample_alias'],'TITLE':r['xml_title'],'地理属性':r['hospital'],'在30人主分析':r['in_current_30_primary']} for r in title_conflicts]
maxerr=max(float(r['max_abs_difference_from_saved_98E']) for r in comp)
report=f'''# PRJNA1125274：132 vs 134 患者编号核对结论

核对日期：2026-09-09。仅使用公共元数据、已有98D计数对象与98E统计结果；没有重新执行DADA2、FASTQ读取、过滤、去噪、拼接或分类。

## 最终裁定

**最终分类：`unresolved`。最可能解释：`metadata bookkeeping issue`（基线人数与全部上传编号组的计数口径不同，或原论文人数报告不一致）。**

已经精确定位差额：全部289个Run映射为134个医院内患者编号组，其中132个编号具有T0。另两个编号为Novara的`SNO_4`和`SNO_34`，共3个晚期样本，均没有T0。132是可直接验证的“T0样本对应的独立编号数”；没有公开临床名册可以进一步证明它就是原文132人的确切名单。

**没有证据支持把两个编号并入其他患者，也没有证据认定其为两个确实额外招募的人。保留原始patient-time map。未生成所谓corrected_patient_time_map，以免把未经证实的身份关系写成更正。**

## 1. 可复核证据链

1. 从指定FINAL map及289-run manifest独立重解析`数字+A/B/C+医院后缀`。289条alias全部满足严格格式，与现有patient_id/timepoint解析结果逐条一致。Run、BioSample、SRS、alias以及patient-timepoint键均没有重复。中心后缀与每个BioSample的地理属性一致。
2. 逐条读取289份本地BioSample XML，并于2026-09-09重新获取ENA全部289条记录。当前ENA的run集合及alias与冻结版本一致；三个关键BioSample的在线XML仍保留相同编号、日期和Novara中心属性。
3. 原始alias中共有134个医院限定编号，但A/T0编号恰好132个：Alessandria 85个、Novara 47个。剩余`SNO_4`只有B/C，`SNO_34`只有C。旧版2026年8月的later-only排除表也已记录相同3条Run，证明差额不是9月7日98C新解析出来的偶发错误。
4. 正式论文摘要和Results 3.1写132人，定植分组94+38也为132；同一Results段落和多处图注的年龄分组却写85名≥65岁和49名＜65岁，合计134。**这是论文自身的数量不一致；85/49的年龄数字不可解释成论文报告的中心分布。**
5. 论文Methods 2.1给出了成人、入院科室、知情同意等纳入标准及无知情同意的排除标准，但未给逐中心人数、筛选/排除人数流程或个体ID清单。当前PMC/Europe PMC XML标记`pmc-prop-has-supplement=no`，无supplementary-material或media节点；Wiley全文也没有可供核对患者名单的Supporting Information入口。没有可公开取得的补充患者名册，不能声称已核对不存在的名册。
6. 论文所引临床注册NCT07183722当前记录为86名ACTUAL enrollment，地点仅Alessandria。这一登记不能直接用于确认多中心132或134人，也没有解决Novara两个编号的对应关系。

证据快照及下载日期见`source_evidence/`和`source_retrieval_log.csv`。原论文 DOI：[10.1002/mbo3.70301](https://onlinelibrary.wiley.com/doi/full/10.1002/mbo3.70301)。公共序列记录：[NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1125274)。注册记录：[NCT07183722](https://clinicaltrials.gov/study/NCT07183722)。

## 2. 医院、时间点及完整病例统计

以下全部是由公开alias和BioSample地理属性重建的统计，不是未经出处确认的论文中心人数。

{md(center_table,list(center_table[0]))}

{md(stages,list(stages[0]))}

37名在原始metadata中具备ABC，应用已有分析对象和2000 reads门槛后30名可进入当前完整病例主分析。失去完整性的7个编号：{', '.join(lost_ids)}。逐Run深度及缺失时间点见`PRJNA1125274_complete_case_attrition_details.csv`。

## 3. 两个额外编号的逐样本证据

{md(extra_table,list(extra_table[0]))}

三个样本均存在于已有98D对象并通过2000 reads门槛。但由于缺少自身T0，98E在建立个人基线后将其排除在displacement表之外。因此它们均不属于当前30名完整病例，也不在T0→T1配对次要分析中。

这些记录只有公开测序alias和采样日期等属性，没有可独立链接到临床纳入名单的patient/subject编号、年龄、性别、定植或脓毒症状态表。其他132个T0编号也并非已经得到逐人临床名册验证。

## 4. 重复、解析和异常检查

- Run、BioSample、正式alias、patient-timepoint键重复：均为0。
- 严格alias解析失败、患者编号解析不一致、医院后缀与地理属性不一致：均为0。
- 49个数字编号在两个医院复用，这是正常的医院内编号空间；`SAL_n`和`SNO_n`已分开，不能按裸数字合并。现有完整patient_id没有跨中心碰撞。
- 单样本编号共16个：{', '.join(single_ids)}。其中仅`SNO_34`无T0。单样本本身不是误标或应删除的证据。
- 日期倒序共12个编号：{', '.join(date_ids)}。原有mapping_status已标记这些异常；日期不具备自动改变患者身份或交换A/B/C的证据。
- 标题字段有4处与正式alias不完全一致：SNO_36的两处只是缺连接号，规范化后身份一致；另两处为下表的跨中心标题冲突。

{md(title_table,list(title_table[0]))}

`SAL_8`/`SNO_8`的正式alias、SUBMITTER_ID和geo相互支持，但TITLE医院后缀相反。这更像提交标题错误，也不能排除源样本标注问题。两名患者本来都已存在，因此交换T2也不能解释134→132；其中SAL_8进入30人主分析，需另行向数据提交者核实。未据此交换样本。

候选合并诊断表按“同中心、没有时间点冲突、日期兼容”列出可能配对。`SNO_4`可以与不止一个基线编号在日期上兼容，`SNO_34`也有多个候选。数字相近或日期连续不能证明同一患者，因此所有候选都标为NOT_ADJUDICATED，合并数为0。

## 5. 核对前后主结果

患者映射没有更正。复用原98E的`primary_contrast`函数，从冻结的患者位移表重新计算全部8组主/敏感性结果；与保存的98E数值最大绝对差为{maxerr:.3g}。对照文件中的“after”明确表示“核对后保留原映射”，不是伪造的132人更正版本。

主contrast定义为每人`D(T0,T2)-D(T0,T1)`。均值Δ的95% CI为配对差t区间；Hedges gz的CI沿用原98E的小样本近似计算，两个区间不要混淆。

{md(main,list(main[0]))}

## 6. 额外的98E统计敏感性复算

“没有进入配对检验”不等于完全不影响上游统计。原98E先在全部268个深度合格样本上筛选ASV，再剔除没有个人T0的记录。这三个样本参与了ASV总计数≥20且至少3样本出现的筛选。

在单独目录`98E_ORPHAN_EXCLUSION_SENSITIVITY/`中，从原98D对象重跑98E，假设剔除这三个来源未决样本；没有重跑98D，也没有更新原对象或原主结果。保留ASV从9373降为9262（111个差异），完整病例仍30人。此为**探索性剔除敏感性分析，不是患者身份更正**。

{md(sensitivity,list(sensitivity[0]))}

四种距离/零值处理及日期异常排除版本的前后效应、CI、p值见`PRJNA1125274_orphan_exclusion_sensitivity_comparison.csv`。是否维持原方向与统计结论应按表内结果表述，不能凭“两个编号未入主分析”直接断言所有数值完全不变。

本次实际复算中，两种主距离的Δ均为正、95% CI均未跨0、paired t p均＜0.05；30人的患者集合保持不变。对这三个样本的探索性剔除没有改变主结论，但不等于已解决源患者身份。

## 7. 裁定类别与待补证据

| 类别 | 本次判断 |
| --- | --- |
| metadata bookkeeping issue | 最可能：132基线编号与134全部上传编号组口径不同；论文自身也有132 vs 85+49的冲突。但未获得临床名册作最终确认。 |
| true extra uploaded patients | 公共库确有两个额外alias组及3个Run；不能据此证明确为另外两名独立临床患者。 |
| mapping error | 未发现本地patient_id解析错误。原数据库TITLE存在独立冲突，但不能据此证明SNO_4/SNO_34身份误标。 |
| unresolved | 最终分类。数值差额已完全定位，真实患者身份及论文人数口径仍未解决。 |

向原作者/数据提交者所需的最小证据：132人去标识纳入名单及中心；`4B-SNO`、`4C-SNO`、`34C-SNO`的临床patient ID与T0对应关系；这些患者是否被原文排除及原因；132总数与85/49年龄分组的校正说明；两个8C样本的alias/TITLE对应确认。本次未联系作者。

## 8. 同时发现的统计报告字段问题

- 原98E accounting的`patients_complete_T0T1T2=37`沿用了QC前complete标记，并非深度后真实完整人数。新`reconciled_population_accounting.csv`按实际timepoint重新计算，深度后为30。
- 原primary表的`n_excluded_anomaly=28`是位移表中异常样本行数，并非排除患者数；30→24表示该敏感性分析实际减少6名完整患者。
- 冻结SAP第2节写日期异常排除于primary paired analysis，但原98E将含异常的30人标作主分析、24人标作敏感性。此处存在SAP与结果命名不一致，须在投稿前说明。此次保持你指定的“当前30人主分析”作为核对基准，没有事后改变终点或主分析人群。

## 9. 文件和复现

- `PRJNA1125274_patient_reconciliation_table.csv`：134行，逐人汇总T0/T1/T2、Run、BioSample、alias、日期、中心、深度和主分析资格以及异常标记。
- `PRJNA1125274_unique_134_patient_ids.csv`：全部134个唯一编号。
- `PRJNA1125274_run_level_reconciliation.csv`：289行逐Run证据。
- `PRJNA1125274_primary_results_before_after_correction.csv`：8组原映射统计复核结果，mapping_changed=FALSE。
- `PRJNA1125274_orphan_exclusion_sensitivity_comparison.csv`：独立探索性剔除敏感性结果。
- 无`PRJNA1125274_corrected_patient_time_map.csv`，因为没有足以更正身份的证据。
- 原metadata、98E脚本/结果及98D分析对象SHA256前后相同。序列处理未执行。
- 源结果目录、data目录和9月8日交付包中的FINAL map及manifest逐文件SHA256一致。

## 附录：全部134个patient_id

Alessandria（85个）：

{', '.join(r['patient_id'] for r in patients if r['hospital']=='ALESSANDRIA')}

Novara（49个）：

{', '.join(r['patient_id'] for r in patients if r['hospital']=='NOVARA')}
'''
(OUT/'PRJNA1125274_132_vs_134_resolution.md').write_text(report,encoding='utf-8')

# Machine-readable field definitions for the wide patient and result tables.
codebook=[
 {'field':'patient_id','meaning':'Hospital-scoped alias group. Not clinically adjudicated real-person ID.'},
 {'field':'clinical_crosswalk_status','meaning':'No public enrollment-to-BioSample clinical crosswalk; NOT_PUBLICLY_AVAILABLE is not evidence of ineligibility.'},
 {'field':'n_samples / n_runs / n_unique_biosamples','meaning':'Public unique records before sequence-depth exclusion; each patient-timepoint has one sample.'},
 {'field':'has_T0 / has_T1 / has_T2','meaning':'Observed alias A/B/C under existing 98C mapping, not an invented missing-timepoint imputation.'},
 {'field':'complete_T0_T1_T2_before_depth','meaning':'All three aliases before depth filtering; total 37.'},
 {'field':'in_current_30_primary','meaning':'All three depth-passing timepoints in saved 98E displacement table; total 30.'},
 {'field':'T0_* / T1_* / T2_*','meaning':'Timepoint-specific run, BioSample, alias and date fields. Blank means absent, not zero or imputed.'},
 {'field':'numeric_id_shared_across_centers','meaning':'Same number reused in different hospitals. Namespaced patient IDs remain distinct; no merge is indicated.'},
 {'field':'strict_alias_and_xml_checks_pass','meaning':'Alias, accession, date, geography and normalized title agree. Title hospital conflicts are retained as flags.'},
 {'field':'mean_contrast_before/after','meaning':'Mean within-patient D(T0,T2)-D(T0,T1), separately per metric; after is unchanged-map verification.'},
 {'field':'mean_contrast_ci_*','meaning':'95% t interval for mean paired contrast.'},
 {'field':'gz_ci_*','meaning':'Original 98E approximate 95% interval for Hedges gz, distinct from mean contrast interval.'},
 {'field':'mapping_changed','meaning':'Always FALSE for this audit. No clinically unsupported remapping was performed.'},
]
write('PRJNA1125274_reconciliation_data_dictionary.csv',codebook)

# Preserve audit scripts alongside their outputs.
scripts=OUT/'scripts'
scripts.mkdir(exist_ok=True)
for name in ['98G_PRJNA1125274_patient_reconciliation.py','98G2_PRJNA1125274_verify_primary_statistics.R','98G3_PRJNA1125274_orphan_exclusion_sensitivity.R','98G4_PRJNA1125274_write_reconciliation_report.py']:
    shutil.copy2(ROOT/'code/03_data_processing'/name,scripts/name)

# Verify original hashes again after all analyses.
hashrows=read(OUT/'input_sha256_before.csv')
for r in hashrows:
    assert hashlib.file_digest(Path(r['path']).open('rb'),'sha256').hexdigest()==r['sha256'],r['path']
assert len(patients)==134 and len({r['patient_id'] for r in patients})==134
assert sum(int(r['n_samples']) for r in patients)==289
assert len(primary_ids)==30 and not primary_ids & {'SNO_4','SNO_34'}
assert all(float(r['max_abs_difference_from_saved_98E'])<1e-10 for r in comp)
assert all(int(r['n_patients_orphan_exclusion'])==int(r['n_patients_original']) for r in sens)
print(json.dumps({'report':str(OUT/'PRJNA1125274_132_vs_134_resolution.md'),'stages':stages,'primary':main,'exploratory_sensitivity':sensitivity,'original_input_hashes_unchanged':len(hashrows)},ensure_ascii=False,indent=2))
