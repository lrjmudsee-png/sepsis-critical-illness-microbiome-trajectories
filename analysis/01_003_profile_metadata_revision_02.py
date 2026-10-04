# -*- coding: utf-8 -*-
from pathlib import Path
from datetime import datetime
import re
import pandas as pd

INV = Path(r"E:\sepsis_project\results\V2_00_inventory")
OUT = Path(r"E:\sepsis_project\results\V2_01_metadata_profile")
OUT.mkdir(parents=True, exist_ok=True)

ROLE_PATTERNS = {
    "patient_id":[r"patient",r"subject",r"participant",r"individual",r"host.?id",r"donor"],
    "sample_id":[r"sample.?id",r"sample_accession",r"biosample",r"specimen.?id"],
    "run_id":[r"run.?id",r"run_accession",r"^accession$",r"^run$",r"sra.?run"],
    "time":[r"time",r"timepoint",r"time_point",r"^day$",r"_day$",r"visit",r"collection.?date",r"sampling.?date",r"icu.?day",r"hospital.?day"],
    "sepsis_group":[r"sepsis",r"diagnos",r"disease",r"group",r"case.?control",r"phenotype"],
    "infection_source":[r"infection.?source",r"infection.?site",r"site.?of.?infection",r"focus.?of.?infection"],
    "outcome":[r"outcome",r"mortality",r"death",r"dead",r"survival",r"survivor",r"deceased"],
    "age":[r"^age$",r"age.?year",r"host.?age"],
    "sex":[r"^sex$",r"gender",r"host.?sex"],
    "antibiotics":[r"antibiotic",r"antimicrobial",r"abx"],
    "severity":[r"sofa",r"apache",r"severity",r"septic.?shock",r"shock"],
}

FILE_COLS=["Project","FileName","FullPath","Exists","LoadStatus","Rows","Columns","DetectedRoles","SuggestedPrimary","Error"]
COL_COLS=["Project","FileName","FullPath","Column","NormalizedColumn","InferredRoles","NonNull","Unique","ExampleValues"]

def latest(pattern):
    fs=sorted(INV.glob(pattern), key=lambda p:p.stat().st_mtime, reverse=True)
    if not fs: raise FileNotFoundError(pattern)
    return fs[0]

def norm(x):
    x=str(x).strip().lower()
    x=re.sub(r"[\s\-/\\\.\(\)\[\]:]+","_",x)
    return re.sub(r"_+","_",x).strip("_")

def infer(c):
    n=norm(c); out=[]
    for role,pats in ROLE_PATTERNS.items():
        if any(re.search(p,n,re.I) for p in pats): out.append(role)
    return out

def read_table(p):
    ext=p.suffix.lower()
    if ext==".csv":
        return pd.read_csv(p,dtype=str,low_memory=False,encoding_errors="replace")
    if ext in {".tsv",".txt"}:
        try:
            df=pd.read_csv(p,sep="\t",dtype=str,low_memory=False,encoding_errors="replace")
            if df.shape[1]>1: return df
        except: pass
        return pd.read_csv(p,sep=None,engine="python",dtype=str,on_bad_lines="skip",encoding_errors="replace")
    if ext in {".xlsx",".xls"}:
        xl=pd.ExcelFile(p); frames=[]
        for s in xl.sheet_names:
            d=pd.read_excel(p,sheet_name=s,dtype=str); d.insert(0,"__sheet__",s); frames.append(d)
        return pd.concat(frames,ignore_index=True,sort=False) if frames else pd.DataFrame()
    raise ValueError(f"Unsupported format: {ext}")

def exvals(s,n=5):
    try:
        v=s.dropna().astype(str).map(str.strip)
        return " | ".join(v[v!=""].drop_duplicates().head(n).tolist())
    except: return ""

def role_columns(cp, project, role):
    if cp.empty or not {"Project","Column","InferredRoles"}.issubset(cp.columns): return ""
    sub=cp[cp["Project"]==project]
    if sub.empty: return ""
    m=sub["InferredRoles"].fillna("").apply(lambda x: role in str(x).split(";"))
    return ";".join(sub.loc[m,"Column"].dropna().astype(str).drop_duplicates().tolist())

def main():
    invp=latest("V2_project_inventory_*.csv")
    candp=latest("V2_metadata_candidates_*.csv")
    inv=pd.read_csv(invp,dtype=str,encoding_errors="replace")
    cand=pd.read_csv(candp,dtype=str,encoding_errors="replace")

    file_rows=[]; col_rows=[]

    for _,r in cand.iterrows():
        project=str(r.get("Project","")).strip()
        raw=str(r.get("FullPath","")).strip()
        if not project or not raw or raw.lower()=="nan": continue
        p=Path(raw)
        base={"Project":project,"FileName":p.name,"FullPath":raw,"Exists":p.exists(),"LoadStatus":"","Rows":None,"Columns":None,"DetectedRoles":"","SuggestedPrimary":"NO","Error":""}
        if not p.exists():
            base["LoadStatus"]="MISSING_FILE"; base["Error"]="Path does not exist"; file_rows.append(base); continue
        try:
            df=read_table(p); df.columns=[str(c) for c in df.columns]
            roles=set()
            for c in df.columns:
                rr=infer(c); roles.update(rr); s=df[c]
                col_rows.append({"Project":project,"FileName":p.name,"FullPath":raw,"Column":c,"NormalizedColumn":norm(c),"InferredRoles":";".join(rr),"NonNull":int(s.notna().sum()),"Unique":int(s.dropna().astype(str).nunique()),"ExampleValues":exvals(s)})
            base.update({"LoadStatus":"OK","Rows":len(df),"Columns":len(df.columns),"DetectedRoles":";".join(sorted(roles))})
        except Exception as e:
            base["LoadStatus"]="ERROR"; base["Error"]=f"{type(e).__name__}: {e}"
        file_rows.append(base)

    fp=pd.DataFrame(file_rows,columns=FILE_COLS)
    cp=pd.DataFrame(col_rows,columns=COL_COLS)

    # choose primary: readable file with most useful roles, then most rows
    if not fp.empty:
        for project in fp["Project"].dropna().unique():
            sub=fp[(fp["Project"]==project)&(fp["LoadStatus"]=="OK")].copy()
            if sub.empty: continue
            sub["score"]=sub["DetectedRoles"].fillna("").apply(lambda x: len([z for z in str(x).split(";") if z]))*100 + pd.to_numeric(sub["Rows"],errors="coerce").fillna(0)
            fp.loc[sub["score"].idxmax(),"SuggestedPrimary"]="YES"

    projects=[]
    for project in inv["Project"].dropna().astype(str).drop_duplicates():
        readable=fp[(fp["Project"]==project)&(fp["LoadStatus"]=="OK")]
        failed=fp[(fp["Project"]==project)&(fp["LoadStatus"]!="OK")]
        primary=readable[readable["SuggestedPrimary"]=="YES"]
        pc={role:role_columns(cp,project,role) for role in ROLE_PATTERNS}
        if readable.empty:
            ready="NO_READABLE_METADATA"
        elif pc["patient_id"] and pc["time"] and (pc["sample_id"] or pc["run_id"]):
            ready="READY_FOR_MANUAL_MAPPING"
        elif pc["time"] and (pc["sample_id"] or pc["run_id"]):
            ready="MISSING_PATIENT_ID"
        elif pc["patient_id"] and (pc["sample_id"] or pc["run_id"]):
            ready="MISSING_TIME"
        else:
            ready="PARTIAL_METADATA"
        projects.append({
            "Project":project,
            "MetadataFilesReadable":len(readable),
            "MetadataFilesFailed":len(failed),
            "SuggestedPrimaryFile":"" if primary.empty else str(primary.iloc[0]["FileName"]),
            "SuggestedPrimaryPath":"" if primary.empty else str(primary.iloc[0]["FullPath"]),
            "PatientColumns":pc["patient_id"],"SampleColumns":pc["sample_id"],"RunColumns":pc["run_id"],"TimeColumns":pc["time"],
            "SepsisGroupColumns":pc["sepsis_group"],"InfectionSourceColumns":pc["infection_source"],"OutcomeColumns":pc["outcome"],
            "AgeColumns":pc["age"],"SexColumns":pc["sex"],"AntibioticColumns":pc["antibiotics"],"SeverityColumns":pc["severity"],
            "Readiness":ready
        })

    pr=pd.DataFrame(projects)
    ts=datetime.now().strftime("%Y%m%d_%H%M%S")
    f1=OUT/f"V2_metadata_file_profile_{ts}.csv"
    f2=OUT/f"V2_metadata_column_profile_{ts}.csv"
    f3=OUT/f"V2_project_metadata_readiness_{ts}.csv"
    f4=OUT/f"V2_metadata_load_errors_{ts}.csv"
    fp.to_csv(f1,index=False,encoding="utf-8-sig")
    cp.to_csv(f2,index=False,encoding="utf-8-sig")
    pr.to_csv(f3,index=False,encoding="utf-8-sig")
    fp[fp["LoadStatus"]!="OK"].to_csv(f4,index=False,encoding="utf-8-sig")
    print("STEP 02 COMPLETE")
    print(pr[["Project","MetadataFilesReadable","MetadataFilesFailed","SuggestedPrimaryFile","Readiness"]].to_string(index=False))
    print(f1); print(f2); print(f3); print(f4)

if __name__=="__main__":
    main()
