"""Patient-level provenance audit only. Never opens FASTQ or invokes DADA2.

Writes an independent reconciliation folder; original maps/results are immutable.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import re
import shutil
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path('E:/sepsis_project')
UP = ROOT / 'results/V2_UPGRADE_20260907'
OUT = UP / '98G_PRJNA1125274_PATIENT_RECONCILIATION_20260909'
SOURCE = OUT / 'source_evidence'
DATA = ROOT / 'data/PRJNA1125274'
M = UP / '98C_PRJNA1125274_METADATA_FREEZE'
E = UP / '98E_PRJNA1125274_EXTERNAL_VALIDATION'
MAP = M / 'PRJNA1125274_FINAL_patient_time_map.csv'
MANIFEST = M / 'PRJNA1125274_FINAL_289_run_manifest.csv'
QC = M / 'PRJNA1125274_metadata_and_run_QC_summary.csv'
XML_DIR = DATA / '00_metadata/step76_sample_xml'
ORPHANS = {'SNO_4', 'SNO_34'}

def read_csv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def save_csv(name, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row))
    with (OUT / name).open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def joined(values):
    return ';'.join(dict.fromkeys(str(x) for x in values if str(x)))

def natural(pid):
    a, n = pid.split('_')
    return (a, int(n))

def get_source(item):
    name, url = item
    dest = SOURCE / name
    rec = {'file': name, 'url': url,
           'retrieved_utc': datetime.now(timezone.utc).isoformat()}
    try:
        if dest.exists():
            data = dest.read_bytes()
            rec['status'] = 'CACHED_THIS_AUDIT'
        else:
            req = urllib.request.Request(url, headers={'User-Agent': 'ResearchMetadataAudit/1.0'})
            with urllib.request.urlopen(req, timeout=35) as response:
                data = response.read()
                rec['http_status'] = response.status
            dest.write_bytes(data)
            rec['status'] = 'DOWNLOADED'
        rec.update(bytes=len(data), sha256=sha(dest))
    except Exception as ex:
        rec.update(status='UNAVAILABLE', error=str(ex))
    return rec

def acquire():
    SOURCE.mkdir(parents=True, exist_ok=True)
    urls = [
        ('article_europepmc.xml', 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13132800/fullTextXML'),
        ('clinicaltrials_NCT07183722.json', 'https://clinicaltrials.gov/api/v2/studies/NCT07183722'),
        ('current_ENA_run_metadata.tsv', 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1125274&result=read_run&fields=run_accession,sample_accession,secondary_sample_accession,sample_alias,sample_title,study_accession,library_strategy,library_layout&format=tsv&download=true'),
    ]
    for bs in ['SAMN41885716', 'SAMN41885724', 'SAMN41885844']:
        urls.append((bs + '_current.xml', 'https://www.ebi.ac.uk/ena/browser/api/xml/' + bs))
    with ThreadPoolExecutor(max_workers=4) as pool:
        log = list(pool.map(get_source, urls))
    save_csv('source_retrieval_log.csv', log)
    print(json.dumps(log, ensure_ascii=False, indent=2))

def audit():
    OUT.mkdir(parents=True, exist_ok=True)
    SOURCE.mkdir(exist_ok=True)
    rows, manifest = read_csv(MAP), read_csv(MANIFEST)
    assert len(rows) == len(manifest) == 289
    assert {x['run_id'] for x in rows} == {x['run_accession'] for x in manifest}
    mr = {x['run_accession']: x for x in manifest}
    disp = read_csv(E / 'PRJNA1125274_external_validation_sample_displacement.csv')
    dr = {r['run_id']: r for r in disp}
    object_csv = OUT / 'analysis_object_run_accounting.csv'
    obj = {r['run_id']: r for r in read_csv(object_csv)} if object_csv.exists() else {}
    current_path = SOURCE / 'current_ENA_run_metadata.tsv'
    current = {}
    if current_path.exists():
        with current_path.open(encoding='utf-8-sig') as f:
            current = {r['run_accession']: r for r in csv.DictReader(f, delimiter='\t')}
    live_xml_checks=[]
    for r in rows:
        p=SOURCE/(r['biosample']+'_current.xml')
        if not p.exists():continue
        s=ET.parse(p).getroot().find('.//SAMPLE')
        at={a.findtext('TAG'):a.findtext('VALUE','') for a in s.findall('./SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE')}
        live_xml_checks.append({'run_id':r['run_id'],'biosample':r['biosample'],
                                'current_alias':s.attrib.get('alias',''),
                                'alias_unchanged':s.attrib.get('alias','')==r['sample_alias'],
                                'date_unchanged':at.get('collection_date')==r['collection_date'],
                                'geo_unchanged':at.get('geo_loc_name')==r['geo'],
                                'source_file':str(p)})
    assert all(r['alias_unchanged'] and r['date_unchanged'] and r['geo_unchanged'] for r in live_xml_checks)
    save_csv('PRJNA1125274_key_BioSample_current_verification.csv',live_xml_checks)
    gp = defaultdict(list)
    for row in rows:
        gp[row['patient_id']].append(row)
    assert len(gp) == 134
    primary = {p for p in gp if {r['timepoint'] for r in disp if r['patient_id'] == p} == {'T0','T1','T2'}}
    assert len(primary) == 30
    assert primary.isdisjoint(ORPHANS)

    duplicate_fields = {f: Counter(r[f] for r in rows) for f in ['run_id','biosample','sample_alias','secondary_sample_accession']}
    pt = Counter((r['patient_id'],r['timepoint']) for r in rows)
    numeric = defaultdict(set)
    for r in rows:
        numeric[int(r['alias_number'])].add(r['hospital'])
    run_audit, attrs_records, xml_inputs = [], [], []
    for r in rows:
        xfile = XML_DIR / (r['biosample'] + '.xml')
        xml_inputs.append(xfile)
        sample = ET.parse(xfile).getroot().find('.//SAMPLE')
        attrs = defaultdict(list)
        for a in sample.findall('./SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE'):
            attrs[a.findtext('TAG')].append(a.findtext('VALUE',''))
        attrs_records.extend({'biosample':r['biosample'],'tag':k,'value':v} for k,vs in attrs.items() for v in vs)
        match = re.fullmatch(r'(\d+)([ABC])-(SAL|SNO)', r['sample_alias'])
        expected_pid = f'{match[3]}_{int(match[1])}' if match else ''
        title = sample.findtext('TITLE','')
        title_alias = re.sub(r'_S\d+$', '', title)
        normalized_title = re.sub(r'^(\d+[ABC])(SAL|SNO)$', r'\1-\2', title_alias)
        title_match_normalized = normalized_title == r['sample_alias']
        date = attrs.get('collection_date', [''])[0]
        geo = attrs.get('geo_loc_name',[''])[0]
        clinical = [k for k in attrs if re.search(r'patient|subject|host_age|host_sex|age|gender|sepsis|disease|outcome|enroll|consent', k, re.I) and k != 'NCBI submission package']
        ra = dict(r)
        ra.update(
            source_xml=str(xfile), xml_primary=sample.findtext('./IDENTIFIERS/PRIMARY_ID',''),
            xml_alias=sample.attrib.get('alias',''), xml_title=title, xml_collection_date=date,
            xml_geo=geo, strict_alias_parse_ok=bool(match), parsed_patient_id=expected_pid,
            patient_id_parse_match=(expected_pid==r['patient_id']),
            timepoint_parse_match=bool(match) and {'A':'T0','B':'T1','C':'T2'}[match[2]]==r['timepoint'],
            title_alias_match=(title_alias==r['sample_alias']),
            title_alias_match_after_hyphen_normalization=title_match_normalized,
            title_mismatch_category=('MATCH' if title_alias==r['sample_alias'] else 'MISSING_HYPHEN_ONLY' if title_match_normalized else 'HOSPITAL_CODE_CONFLICT'),
            xml_alias_match=(sample.attrib.get('alias','')==r['sample_alias']),
            xml_date_match=(date==r['collection_date']),
            xml_biosample_match=(sample.findtext('./IDENTIFIERS/PRIMARY_ID','')==r['biosample']),
            xml_secondary_match=(sample.findtext('./IDENTIFIERS/SECONDARY_ID','')==r['secondary_sample_accession']),
            geo_hospital_match=r['hospital'].lower() in geo.lower(),
            duplicate_patient_timepoint=pt[(r['patient_id'],r['timepoint'])]>1,
            duplicate_run=duplicate_fields['run_id'][r['run_id']]>1,
            duplicate_biosample=duplicate_fields['biosample'][r['biosample']]>1,
            duplicate_alias=duplicate_fields['sample_alias'][r['sample_alias']]>1,
            numeric_id_shared_across_hospitals=len(numeric[int(r['alias_number'])])>1,
            manifest_mapping_match=all(r[f]==mr[r['run_id']][f] for f in ['patient_id','sample_alias','biosample','timepoint','hospital','collection_date']),
            current_ena_run_present=r['run_id'] in current if current else 'NOT_CHECKED',
            current_ena_alias_match=current.get(r['run_id'],{}).get('sample_alias')==r['sample_alias'] if current else 'NOT_CHECKED',
            current_ena_biosample_match=current.get(r['run_id'],{}).get('sample_accession')==r['biosample'] if current else 'NOT_CHECKED',
            current_ena_secondary_match=current.get(r['run_id'],{}).get('secondary_sample_accession')==r['secondary_sample_accession'] if current else 'NOT_CHECKED',
            explicit_clinical_attribute_tags=joined(clinical),
            clinical_patient_crosswalk_status='NOT_PUBLICLY_AVAILABLE',
            in_analysis_object=r['run_id'] in obj if obj else 'NOT_CHECKED',
            nonchim_reads=obj.get(r['run_id'],{}).get('depth',''),
            passes_2000_reads=obj.get(r['run_id'],{}).get('passes_depth',''),
            in_displacement_table=r['run_id'] in dr,
            in_current_30_primary=r['patient_id'] in primary,
            no_T0_alias_group=r['patient_id'] in ORPHANS)
        run_audit.append(ra)
    save_csv('PRJNA1125274_run_level_reconciliation.csv', run_audit)
    save_csv('PRJNA1125274_BioSample_attributes_audit.csv',attrs_records)
    rmap={r['run_id']:r for r in run_audit}
    patients=[]
    anomaly_pairs=[]
    for pid in sorted(gp,key=natural):
        rs=sorted(gp[pid],key=lambda x:x['timepoint'])
        ars=[rmap[r['run_id']] for r in rs]
        tps={r['timepoint'] for r in rs}
        dates={r['timepoint']:r['collection_date'] for r in rs}
        reversed_dates=[a+'>'+b for a,b in [('T0','T1'),('T0','T2'),('T1','T2')] if a in dates and b in dates and dates[a]>dates[b]]
        same_dates=[a+'='+b for a,b in [('T0','T1'),('T0','T2'),('T1','T2')] if a in dates and b in dates and dates[a]==dates[b]]
        allsame=len(set(dates.values()))==1 and len(dates)>1
        flags=[]
        if 'T0' not in tps: flags.append('NO_T0_NO_CLINICAL_CROSSWALK')
        if len(rs)==1: flags.append('SINGLE_SAMPLE')
        if reversed_dates: flags.append('DATE_ORDER_ANOMALY')
        if same_dates: flags.append('SAME_DAY_TIMEPOINTS')
        if any(not a['patient_id_parse_match'] or not a['title_alias_match_after_hyphen_normalization'] or not a['xml_alias_match'] for a in ars): flags.append('SOURCE_TITLE_HOSPITAL_CONFLICT')
        if any(a['title_mismatch_category']=='MISSING_HYPHEN_ONLY' for a in ars): flags.append('TITLE_MISSING_HYPHEN_ONLY')
        if any(a['duplicate_patient_timepoint'] for a in ars): flags.append('DUPLICATE_PATIENT_TIMEPOINT')
        p={'patient_id':pid,'hospital':rs[0]['hospital'],'alias_number':int(rs[0]['alias_number']),
           'n_samples':len(rs),'n_unique_biosamples':len({r['biosample'] for r in rs}),
           'n_runs':len({r['run_id'] for r in rs}),'timepoint_pattern':joined(sorted(tps)),
           'has_T0':'T0' in tps,'has_T1':'T1' in tps,'has_T2':'T2' in tps,
           'complete_T0_T1_T2_before_depth':len(tps)==3,'in_current_30_primary':pid in primary,
           'n_depth_passing_samples':sum(a['passes_2000_reads']=='TRUE' for a in ars),
           'n_displacement_rows':sum(a['in_displacement_table'] for a in ars),
           'sample_aliases':joined(r['sample_alias'] for r in rs),
           'biosamples':joined(r['biosample'] for r in rs),'run_accessions':joined(r['run_id'] for r in rs),
           'collection_dates':joined(r['timepoint']+':'+r['collection_date'] for r in rs),
           'date_order_anomaly':bool(reversed_dates),'date_order_pairs':joined(reversed_dates),
           'same_day_pairs':joined(same_dates),'all_timepoints_same_date':allsame,
           'single_sample_patient':len(rs)==1,
           'numeric_id_shared_across_centers':any(a['numeric_id_shared_across_hospitals'] for a in ars),
           'cross_center_patient_key_conflict':len({r['hospital'] for r in rs})>1,
           'strict_alias_and_xml_checks_pass':all(a['strict_alias_parse_ok'] and a['patient_id_parse_match'] and a['title_alias_match_after_hyphen_normalization'] and a['xml_alias_match'] and a['xml_biosample_match'] and a['xml_date_match'] and a['geo_hospital_match'] for a in ars),
           'clinical_crosswalk_status':'NOT_PUBLICLY_AVAILABLE',
           'relationship_to_paper_132':'T0_OBSERVED_ID_NUMERICALLY_MATCHES_BASELINE_COUNT' if 'T0' in tps else 'EXTRA_ALIAS_GROUP_RELATIVE_TO_132_T0_IDS_UNRESOLVED',
           'flags':joined(flags),'mapping_correction_applied':False,
           'adjudication':'UNRESOLVED_CLINICAL_IDENTITY_KEEP_ORIGINAL' if pid in ORPHANS else 'NO_IDENTITY_CORRECTION_SUPPORTED'}
        for t in ['T0','T1','T2']:
            sub=[r for r in rs if r['timepoint']==t]
            p.update({t+'_n_samples':len(sub),t+'_biosamples':joined(r['biosample'] for r in sub),
                      t+'_runs':joined(r['run_id'] for r in sub),t+'_aliases':joined(r['sample_alias'] for r in sub),
                      t+'_dates':joined(r['collection_date'] for r in sub)})
        patients.append(p)
    save_csv('PRJNA1125274_patient_reconciliation_table.csv',patients)
    save_csv('PRJNA1125274_unique_134_patient_ids.csv',({'patient_id':p['patient_id']} for p in patients))
    save_csv('PRJNA1125274_flagged_patients.csv',(p for p in patients if p['flags']))
    save_csv('PRJNA1125274_extra_alias_groups.csv',(r for r in run_audit if r['no_T0_alias_group']))
    save_csv('PRJNA1125274_title_alias_discrepancies.csv',(r for r in run_audit if not r['title_alias_match']))
    centers=[]
    for h in ['ALESSANDRIA','NOVARA','TOTAL']:
        ps=[p for p in patients if h=='TOTAL' or p['hospital']==h]
        rs=[r for r in run_audit if h=='TOTAL' or r['hospital']==h]
        centers.append({'hospital':h,'n_patient_alias_groups':len(ps),'n_samples':len(rs),
                        'n_T0_patients':sum(p['has_T0'] for p in ps),
                        'n_no_T0_patients':sum(not p['has_T0'] for p in ps),
                        'samples_T0':sum(r['timepoint']=='T0' for r in rs),
                        'samples_T1':sum(r['timepoint']=='T1' for r in rs),
                        'samples_T2':sum(r['timepoint']=='T2' for r in rs),
                        'complete_before_depth':sum(p['complete_T0_T1_T2_before_depth'] for p in ps),
                        'complete_primary_after_depth':sum(p['in_current_30_primary'] for p in ps),
                        'single_sample_patients':sum(p['single_sample_patient'] for p in ps),
                        'paper_hospital_patient_count':'NOT_REPORTED_IN_ACCESSIBLE_ARTICLE',
                        'count_source':'RECONSTRUCTED_FROM_PUBLIC_ALIASES_AND_GEO'})
    save_csv('PRJNA1125274_hospital_patient_sample_summary.csv',centers)
    cross=[{'alias_number':n,'hospitals':joined(sorted(hs)),
            'patient_ids':joined(sorted({r['patient_id'] for r in rows if int(r['alias_number'])==n},key=natural)),
            'decision':'KEEP_SEPARATE_HOSPITAL_NAMESPACE'} for n,hs in sorted(numeric.items()) if len(hs)>1]
    save_csv('PRJNA1125274_cross_center_numeric_ids.csv',cross)

    # Flag candidate typographical merges conservatively; no candidate is adopted.
    candidates=[]
    for orphan in patients:
        if orphan['patient_id'] not in ORPHANS: continue
        src=gp[orphan['patient_id']]
        st={r['timepoint'] for r in src}
        for dest in patients:
            if dest['patient_id'] in ORPHANS or dest['hospital']!=orphan['hospital']:continue
            drs=gp[dest['patient_id']]
            dt={r['timepoint'] for r in drs}
            if st & dt: continue
            dates=sorted([(r['timepoint'],r['collection_date']) for r in src+drs])
            chron=all(dates[i][1]<=dates[i+1][1] for i in range(len(dates)-1))
            a,b=str(orphan['alias_number']),str(dest['alias_number'])
            edit_one=(len(a)==len(b) and sum(x!=y for x,y in zip(a,b))==1) or (len(a)==len(b)+1 and any(a[:i]+a[i+1:]==b for i in range(len(a)))) or (len(b)==len(a)+1 and any(b[:i]+b[i+1:]==a for i in range(len(b))))
            if chron:
                candidates.append({'orphan_patient_id':orphan['patient_id'],'candidate_patient_id':dest['patient_id'],
                                   'candidate_aliases':dest['sample_aliases'],'candidate_dates':dest['collection_dates'],
                                   'no_timepoint_collision':True,'chronology_compatible':chron,'numeric_id_edit_distance_one':edit_one,
                                   'identity_evidence':'NONE_DATE_AND_NUMBER_SIMILARITY_IS_NOT_IDENTITY',
                                   'merge_applied':False})
    save_csv('PRJNA1125274_potential_merge_candidates_NOT_ADJUDICATED.csv',candidates,
             ['orphan_patient_id','candidate_patient_id','candidate_aliases','candidate_dates','no_timepoint_collision','chronology_compatible','numeric_id_edit_distance_one','identity_evidence','merge_applied'])
    # Freeze source hashes, including existing analysis object. Never open raw sequences.
    fixed=[MAP,MANIFEST,QC,E/'PRJNA1125274_external_validation_primary_contrast.csv',
           E/'PRJNA1125274_external_validation_sample_displacement.csv',
           ROOT/'code/03_data_processing/98E_PRJNA1125274_external_validation.R',
           DATA/'03_dada2/PRJNA1125274_analysis_object_external_validation.rds',
           ROOT/'PMC13132800_fulltext.xml']
    hashes=[{'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)} for p in fixed+xml_inputs]
    hashfile=OUT/'input_sha256_before.csv'
    if not hashfile.exists():save_csv('input_sha256_before.csv',hashes)
    old={r['path']:r for r in read_csv(hashfile)}
    assert all(old[r['path']]['sha256']==r['sha256'] for r in hashes)
    save_csv('input_sha256_after.csv',hashes)
    snap=[MAP,MANIFEST,QC,E/'PRJNA1125274_external_validation_primary_contrast.csv']
    for p in snap: shutil.copy2(p,SOURCE/p.name)
    consistency=[]
    for name in [MAP.name,MANIFEST.name]:
        authoritative=M/name
        copy_candidates=[DATA/name,ROOT/'V2_UPGRADE_20260908_DELIVERABLE/02_Results/98C_PRJNA1125274_METADATA_FREEZE'/name]
        for other in copy_candidates:
            consistency.append({'authoritative_path':str(authoritative),'comparison_path':str(other),
                                'exists':other.exists(),'sha256_identical':other.exists() and sha(authoritative)==sha(other)})
    save_csv('PRJNA1125274_input_copy_consistency.csv',consistency)
    paper=ET.parse(ROOT/'PMC13132800_fulltext.xml').getroot()
    paperfacts={
        'title':paper.findtext('./front/article-meta/title-group/article-title'),
        'supplementary_material_nodes':len(paper.findall('.//supplementary-material')),
        'media_nodes':len(paper.findall('.//media')),
        'supplement_property':[ ''.join(e.itertext()) for e in paper.findall('.//custom-meta') if e.findtext('meta-name')=='pmc-prop-has-supplement'],
        'relevant_sections':[]}
    for sec in paper.findall('.//body//sec'):
        for para in sec.findall('./p'):
            txt=''.join(para.itertext())
            if re.search(r'132 patients|132 hospitalized|85 patients|Participants were recruited|Samples were collected',txt):
                paperfacts['relevant_sections'].append({'section_id':sec.attrib.get('id'), 'section_title':''.join(sec.find('title').itertext()) if sec.find('title') is not None else '', 'text':txt})
    (OUT/'article_audit_evidence.json').write_text(json.dumps(paperfacts,ensure_ascii=False,indent=2),encoding='utf-8')
    flag_counts=Counter(f for p in patients for f in p['flags'].split(';') if f)
    checkcols=['strict_alias_parse_ok','patient_id_parse_match','timepoint_parse_match','title_alias_match','xml_alias_match','xml_date_match','xml_biosample_match','xml_secondary_match','geo_hospital_match','manifest_mapping_match']
    errors={k:sum(not r[k] for r in run_audit) for k in checkcols}
    summary={'patients':len(patients),'runs':len(rows),'centers':centers,'orphan_ids':sorted(ORPHANS,key=natural),
             'primary_ids':sorted(primary,key=natural),'flag_counts':dict(flag_counts),'source_check_errors':errors,
             'cross_center_reused_numbers':len(cross),'merge_candidates':candidates,
             'clinical_attribute_tags':sorted({r['tag'] for r in attrs_records if re.search(r'host_age|host_sex|patient|subject|disease|sepsis',r['tag'],re.I)}),
             'current_ena_runs':len(current),
             'current_ena_run_set_identical':set(current)=={r['run_id'] for r in rows} if current else None,
             'current_ena_alias_mismatches':sum(r['current_ena_alias_match'] is False for r in run_audit),
             'current_ena_biosample_mismatches':sum(r['current_ena_biosample_match'] is False for r in run_audit),
             'current_ena_secondary_mismatches':sum(r['current_ena_secondary_match'] is False for r in run_audit),
             'duplicate_run_count':sum(r['duplicate_run'] for r in run_audit),
             'duplicate_biosample_count':sum(r['duplicate_biosample'] for r in run_audit),
             'duplicate_patient_timepoint_count':sum(r['duplicate_patient_timepoint'] for r in run_audit),
             'applied_patient_identity_corrections':0,
             'final_classification':'unresolved',
             'most_likely_explanation':'metadata bookkeeping issue (baseline denominator / publication inconsistency); individual clinical identities unresolved'}
    (OUT/'reconciliation_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['primary_ids','merge_candidates']},ensure_ascii=False,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--acquire', action='store_true')
    args=ap.parse_args()
    if args.acquire: acquire()
    audit()
