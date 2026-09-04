# -*- coding: utf-8 -*-
from __future__ import annotations

import re
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
from datetime import datetime

import pandas as pd

DATA_ROOT = Path(r"E:\sepsis_project\data")
OUT_ROOT = Path(r"E:\sepsis_project\results\V2_05_metadata_coverage")
OUT_ROOT.mkdir(parents=True, exist_ok=True)

TARGET_PROJECTS = [
    "PRJEB37289",
    "PRJEB67798",
    "PRJEB82425",
    "PRJNA516701",
    "PRJNA578267",
    "PRJNA595346",
    "PRJNA884103",
]

DOMAINS = {
    "patient_id": [r"\bpatient\b", r"\bsubject\b", r"participant", r"patient.?id", r"subject.?id", r"individual.?id"],
    "sample_id": [r"sample.?id", r"sample_accession", r"biosample", r"sequence.?name", r"sample.?name", r"specimen"],
    "run_id": [r"run_accession", r"sra.?run", r"\brun.?id\b", r"\bcrr\d+", r"\bsrr\d+", r"\berr\d+"],
    "time": [r"time.?point", r"timepoint", r"\bday\b", r"\bdol\b", r"collection.?date", r"sampling.?date", r"visit", r"baseline", r"discharge", r"intubation", r"post.?antibiotic", r"post.?admission"],
    "infection_group": [r"\bsepsis\b", r"septic.?shock", r"\bvap\b", r"pneumonia", r"infected", r"infection.?group", r"bacteremia", r"bacteraemia", r"bloodstream.?infection"],
    "infection_source": [r"infection.?source", r"infection.?site", r"site.?of.?infection", r"abdominal.?origin", r"pulmonary.?infection", r"urinary.?infection", r"bloodstream.?infection"],
    "outcome": [r"\boutcome\b", r"mortality", r"\bdeath\b", r"\bdied\b", r"survival", r"survivor", r"deceased", r"icu.?length.?of.?stay", r"hospital.?length.?of.?stay", r"ventilator.?free", r"adverse.?outcome"],
    "antibiotics": [r"antibiotic", r"antimicrobial", r"\babx\b", r"cefazolin", r"clindamycin", r"tobramycin", r"colistin", r"cefotaxime", r"antibiotic.?score"],
    "severity": [r"\bsofa\b", r"\bapache\b", r"severity", r"shock", r"vasopressor", r"mechanical.?ventilation", r"ventilation.?duration"],
    "demographics": [r"\bage\b", r"\bsex\b", r"gender", r"birth.?weight", r"gestational.?age", r"\bbmi\b"],
    "body_site": [r"body.?site", r"sample.?type", r"rectal", r"stool", r"faecal", r"fecal", r"tracheal", r"endotracheal", r"gut", r"lung"],
}

TABULAR_EXTS = {".csv",".tsv",".txt",".xlsx",".xls"}
TEXT_EXTS = {".md",".rmd",".r",".py",".json",".xml"}
UNPARSED_EXTS = {".pdf",".rdata",".rds",".biom",".zip",".gz"}
MAX_TEXT_CHARS = 2_000_000
MAX_TABULAR_ROWS = 3000

def detect_domains(text: str):
    t = str(text).lower()
    out = {}
    for domain, pats in DOMAINS.items():
        out[domain] = any(re.search(p, t, re.I) for p in pats)
    return out

def read_tabular(path: Path):
    ext = path.suffix.lower()

    if ext == ".csv":
        try:
            df = pd.read_csv(path, dtype=str, nrows=MAX_TABULAR_ROWS, low_memory=False, encoding_errors="replace")
            return ({"sheet":"","df":df},)
        except Exception:
            return tuple()

    if ext in {".tsv",".txt"}:
        try:
            df = pd.read_csv(path, sep="\t", dtype=str, nrows=MAX_TABULAR_ROWS, low_memory=False, encoding_errors="replace")
            if df.shape[1] > 1:
                return ({"sheet":"","df":df},)
        except Exception:
            pass
        try:
            df = pd.read_csv(path, sep=None, engine="python", dtype=str, nrows=MAX_TABULAR_ROWS, on_bad_lines="skip", encoding_errors="replace")
            return ({"sheet":"","df":df},)
        except Exception:
            return tuple()

    if ext in {".xlsx",".xls"}:
        try:
            xl = pd.ExcelFile(path)
        except Exception:
            return tuple()
        out = []
        for sheet in xl.sheet_names:
            try:
                df = pd.read_excel(path, sheet_name=sheet, dtype=str, nrows=MAX_TABULAR_ROWS)
                out.append({"sheet":sheet,"df":df})
            except Exception:
                pass
        return tuple(out)

    return tuple()

def docx_text(path: Path):
    try:
        with zipfile.ZipFile(path) as z:
            parts = []
            for member in z.namelist():
                if member.startswith("word/") and member.endswith(".xml"):
                    try:
                        root = ET.fromstring(z.read(member))
                        parts.extend(x.text for x in root.iter() if x.text)
                    except Exception:
                        pass
            return " ".join(parts)[:MAX_TEXT_CHARS]
    except Exception:
        return ""

def text_file(path: Path):
    try:
        return path.read_text(encoding="utf-8", errors="replace")[:MAX_TEXT_CHARS]
    except Exception:
        return ""

def project_files(project: str):
    base = DATA_ROOT / project / "00_metadata"
    if not base.exists():
        return []
    return [p for p in base.rglob("*") if p.is_file()]

def main():
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    file_rows = []
    project_domain_sources = {p:{d:[] for d in DOMAINS} for p in TARGET_PROJECTS}
    project_structured_sources = {p:{d:[] for d in DOMAINS} for p in TARGET_PROJECTS}
    project_unparsed = {p:[] for p in TARGET_PROJECTS}

    for project in TARGET_PROJECTS:
        files = project_files(project)
        print(f"Scanning {project}: {len(files)} files")

        for path in files:
            ext = path.suffix.lower()
            rel = str(path.relative_to(DATA_ROOT / project))
            detected = {d:False for d in DOMAINS}
            structured = {d:False for d in DOMAINS}
            evidence = []
            parse_status = "UNPARSED"

            if ext in TABULAR_EXTS:
                tables = read_tabular(path)
                if tables:
                    parse_status = "TABULAR_PARSED"

                for item in tables:
                    df = item["df"]
                    cols = [str(c) for c in df.columns]
                    hdet = detect_domains(" | ".join(cols))

                    vals = []
                    try:
                        for c in df.columns[:100]:
                            vals.extend(
                                df[c].dropna().astype(str).drop_duplicates().head(20).tolist()
                            )
                    except Exception:
                        pass
                    vdet = detect_domains(" | ".join(vals)[:MAX_TEXT_CHARS])

                    for d in DOMAINS:
                        if hdet[d] or vdet[d]:
                            detected[d] = True
                        if hdet[d]:
                            structured[d] = True

                    prefix = f"sheet={item['sheet']}; " if item["sheet"] else ""
                    evidence.append(prefix + "columns=" + ";".join(cols[:40]))

            elif ext in TEXT_EXTS:
                parse_status = "TEXT_PARSED"
                ddet = detect_domains(text_file(path))
                for d in DOMAINS:
                    detected[d] = ddet[d]

            elif ext == ".docx":
                parse_status = "DOCX_TEXT_PARSED"
                ddet = detect_domains(docx_text(path))
                for d in DOMAINS:
                    detected[d] = ddet[d]

            elif ext in UNPARSED_EXTS:
                project_unparsed[project].append(rel)
                parse_status = "PRESENT_NOT_PARSED"

            filename_det = detect_domains(path.name)
            for d in DOMAINS:
                if filename_det[d]:
                    detected[d] = True
                if detected[d]:
                    project_domain_sources[project][d].append(rel)
                if structured[d]:
                    project_structured_sources[project][d].append(rel)

            file_rows.append({
                "Project": project,
                "RelativePath": rel,
                "FileName": path.name,
                "Extension": ext,
                "SizeKB": round(path.stat().st_size/1024,2),
                "ParseStatus": parse_status,
                "DetectedDomains": ";".join(d for d,v in detected.items() if v),
                "StructuredDomains": ";".join(d for d,v in structured.items() if v),
                "EvidencePreview": " || ".join(evidence)[:3000],
            })

    coverage_rows = []
    manual_rows = []
    priority_rows = []
    CORE = ["patient_id","sample_id","time"]
    CLINICAL = ["infection_group","infection_source","outcome","antibiotics","severity"]

    for project in TARGET_PROJECTS:
        row = {"Project":project}
        structured_count = 0
        detected_count = 0

        for d in DOMAINS:
            detected = bool(project_domain_sources[project][d])
            structured = bool(project_structured_sources[project][d])
            if structured:
                status = "STRUCTURED"
                structured_count += 1
            elif detected:
                status = "PRESENT_UNSTRUCTURED"
            else:
                status = "NOT_FOUND"
            if detected:
                detected_count += 1
            row[f"{d}_status"] = status
            row[f"{d}_sources"] = " | ".join(list(dict.fromkeys(project_domain_sources[project][d]))[:10])

        core_structured = all(row[f"{d}_status"]=="STRUCTURED" for d in CORE)
        if core_structured:
            mapping_ready = "YES"
        elif all(row[f"{d}_status"]!="NOT_FOUND" for d in CORE):
            mapping_ready = "MANUAL_LINKAGE_REVIEW"
        else:
            mapping_ready = "NO"

        row["CoreMappingReady"] = mapping_ready
        row["StructuredDomainCount"] = structured_count
        row["DetectedDomainCount"] = detected_count
        row["UnparsedPotentialFiles"] = " | ".join(project_unparsed[project][:20])
        coverage_rows.append(row)

        missing_core = [d for d in CORE if row[f"{d}_status"]=="NOT_FOUND"]
        missing_clin = [d for d in CLINICAL if row[f"{d}_status"]=="NOT_FOUND"]
        unstructured = [d for d in CORE+CLINICAL if row[f"{d}_status"]=="PRESENT_UNSTRUCTURED"]

        if missing_core:
            level = "HIGH"
            action = "MANUAL_COLLECTION_REQUIRED: obtain patient/sample/time linkage from supplement, repository, or authors."
        elif unstructured:
            level = "MEDIUM"
            action = "NO NEW DATA COLLECTION YET: parse/link existing supplement or repository resources first."
        elif missing_clin:
            level = "LOW"
            action = "CORE MAPPING AVAILABLE; collect missing clinical covariates only if required by the planned model."
        else:
            level = "NONE"
            action = "NO MANUAL COLLECTION CURRENTLY REQUIRED."

        if project == "PRJEB67798" and row["patient_id_status"]!="STRUCTURED":
            level = "HIGH"
            action = "AFTER PARSING EXISTING SUPPLEMENT: if patient-level linkage is absent, author contact is likely required because public source data exclude patient data."

        manual_rows.append({
            "Project":project,
            "Priority":level,
            "MissingCore":";".join(missing_core),
            "MissingClinical":";".join(missing_clin),
            "PresentButUnstructured":";".join(unstructured),
            "RecommendedAction":action,
        })

        score = (100 if core_structured else 0) + structured_count*10 + detected_count
        priority_rows.append({
            "Project":project,
            "MappingReadiness":mapping_ready,
            "CoverageScore":score,
            "RecommendedNextStep":(
                "BUILD_SAMPLE_PATIENT_TIME_MAP" if core_structured
                else "PARSE_EXISTING_RESOURCES" if not missing_core
                else "FILL_CORE_MAPPING_GAPS"
            ),
        })

    coverage = pd.DataFrame(coverage_rows)
    files = pd.DataFrame(file_rows)
    manual = pd.DataFrame(manual_rows)
    priority = pd.DataFrame(priority_rows).sort_values(["CoverageScore","Project"], ascending=[False,True])

    p1 = OUT_ROOT / f"V2_metadata_coverage_audit_{stamp}.csv"
    p2 = OUT_ROOT / f"V2_metadata_file_inventory_{stamp}.csv"
    p3 = OUT_ROOT / f"V2_manual_collection_required_{stamp}.csv"
    p4 = OUT_ROOT / f"V2_mapping_priority_{stamp}.csv"

    coverage.to_csv(p1,index=False,encoding="utf-8-sig")
    files.to_csv(p2,index=False,encoding="utf-8-sig")
    manual.to_csv(p3,index=False,encoding="utf-8-sig")
    priority.to_csv(p4,index=False,encoding="utf-8-sig")

    print()
    print("================ STEP 05 COMPLETE ================")
    print(manual[["Project","Priority","MissingCore","MissingClinical","RecommendedAction"]].to_string(index=False))
    print()
    print(p1)
    print(p2)
    print(p3)
    print(p4)

if __name__ == "__main__":
    main()
