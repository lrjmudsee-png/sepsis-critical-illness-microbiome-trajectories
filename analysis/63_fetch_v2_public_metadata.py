# -*- coding: utf-8 -*-
r"""
Sepsis V2 - Step 03 Metadata Collector
Automatically fetch public repository metadata and public supplementary resources.

Targets:
PRJEB37289
PRJEB67798
PRJEB82425
PRJNA516701
PRJNA578267
PRJNA595346
PRJNA884103

Default data root:
E:\sepsis_project\data

For each project, outputs are saved under:
<project>\00_metadata\auto_fetched\

Also writes a collection report under:
E:\sepsis_project\results\V2_03_metadata_collection\

No raw sequencing files are modified.

Run:
python "E:\sepsis_project\code\03_data_processing\63_fetch_v2_public_metadata.py"
"""

from __future__ import annotations

import csv
import io
import json
import os
import re
import shutil
import tarfile
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path
from datetime import datetime


DATA_ROOT = Path(r"E:\sepsis_project\data")
RESULT_ROOT = Path(r"E:\sepsis_project\results\V2_03_metadata_collection")
RESULT_ROOT.mkdir(parents=True, exist_ok=True)

USER_AGENT = "SepsisV2MetadataCollector/1.0"
NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "")
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")

PROJECTS = {
    "PRJEB37289": {
        "paper_title": "Integrative Transkingdom Analysis of the Gut Microbiome in Antibiotic Perturbation and Critical Illness",
        "doi": "10.1128/msystems.01148-20",
        "pmcid": "PMC8546997",
        "github": "Bwhaak/MOFA_microbiome",
        "notes": "Processed data links are provided through the public analysis repository."
    },
    "PRJEB67798": {
        "paper_title": "Intestinal dysbiosis as an intraoperative predictor of septic complications: evidence from human surgical cohorts and preclinical models of peritoneal sepsis",
        "doi": "10.1038/s41598-023-49034-z",
        "pmcid": "PMC10739899",
        "github": "",
        "notes": "BORIS public source data exclude patient data; supplementary XLSX is downloaded separately when possible."
    },
    "PRJEB82425": {
        "paper_title": "Lung and gut microbiota profiling in intensive care unit patients: a prospective pilot study",
        "doi": "10.1186/s12879-025-10825-6",
        "pmcid": "PMC11972518",
        "github": "",
        "notes": "Supplement contains infectious episode and antibiotic/sampling information."
    },
    "PRJNA516701": {
        "paper_title": "Rectal Swabs from Critically Ill Patients Provide Discordant Representations of the Gut Microbiome Compared to Stool Samples",
        "doi": "10.1128/msphere.00358-19",
        "pmcid": "PMC6656869",
        "github": "",
        "notes": "Study used baseline, middle and late ICU sampling intervals."
    },
    "PRJNA578267": {
        "paper_title": "Marked Changes in Gut Microbiota in Cardio-Surgical Intensive Care Unit Patients: A Longitudinal Cohort Study",
        "doi": "",
        "pmcid": "PMC6974539",
        "github": "",
        "notes": "Public supplementary material includes antibiotic-use details."
    },
    "PRJNA595346": {
        "paper_title": "Longitudinal multicompartment characterization of host-microbiota interactions in patients with acute respiratory failure",
        "doi": "10.1038/s41467-024-48819-8",
        "pmcid": "PMC11148165",
        "github": "MicrobiomeALIR/MultiCompartmentMicrobiome",
        "notes": "De-identified clinical and processed microbiome data are public in the authors' GitHub repository."
    },
    "PRJNA884103": {
        "paper_title": "Gut pathogen colonization precedes bloodstream infection in the neonatal intensive care unit",
        "doi": "10.1126/scitranslmed.adg5562",
        "pmcid": "PMC10259202",
        "github": "DJSchwartzLab/NICUBSI",
        "notes": "The paper states that clinical metadata are available under this NCBI BioProject."
    },
}

# Known direct public resources not reliably exposed by PMC package alone.
DIRECT_FILES = {
    "PRJEB67798": [
        (
            "07_known_public_supplement",
            "PRJEB67798_Supplementary_Information_2.xlsx",
            "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41598-023-49034-z/MediaObjects/41598_2023_49034_MOESM2_ESM.xlsx"
        ),
        (
            "07_known_public_supplement",
            "PRJEB67798_Supplementary_Information_1.pdf",
            "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41598-023-49034-z/MediaObjects/41598_2023_49034_MOESM1_ESM.pdf"
        ),
    ],
}


ENA_BASE = "https://www.ebi.ac.uk/ena/portal/api"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PMC_OA_API = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi"


def http_get(url: str, timeout: int = 120, binary: bool = False):
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": USER_AGENT,
            "Accept": "*/*",
        }
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        data = resp.read()
    return data if binary else data.decode("utf-8", errors="replace")


def save_text(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8-sig")


def save_binary(path: Path, data: bytes):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)


def download_file(url: str, out: Path, timeout: int = 300):
    out.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as r, open(out, "wb") as f:
        shutil.copyfileobj(r, f)
    return out


def ncbi_params(extra: dict):
    p = dict(extra)
    p["tool"] = "SepsisV2MetadataCollector"
    if NCBI_EMAIL:
        p["email"] = NCBI_EMAIL
    if NCBI_API_KEY:
        p["api_key"] = NCBI_API_KEY
    return p


def ncbi_url(endpoint: str, params: dict):
    return f"{EUTILS_BASE}/{endpoint}?" + urllib.parse.urlencode(ncbi_params(params))


def sleep_ncbi():
    # Conservative rate limit without API key.
    time.sleep(0.38 if not NCBI_API_KEY else 0.12)


def ena_return_fields(result: str):
    url = f"{ENA_BASE}/returnFields?" + urllib.parse.urlencode({
        "dataPortal": "ena",
        "format": "tsv",
        "result": result
    })
    text = http_get(url)
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    return [r["columnId"] for r in rows if r.get("columnId")]


def ena_read_run(project: str, outdir: Path):
    wanted = [
        "study_accession","secondary_study_accession",
        "sample_accession","secondary_sample_accession",
        "sample_alias","sample_title","sample_description",
        "run_accession","run_alias",
        "experiment_accession","experiment_alias","experiment_title",
        "collection_date","age","sex","host","host_scientific_name",
        "host_sex","host_status","host_body_site","host_phenotype",
        "disease","isolation_source","tissue_type","country",
        "library_strategy","library_source","library_selection",
        "library_layout","instrument_platform","instrument_model",
        "fastq_ftp","fastq_md5","fastq_bytes"
    ]
    available = set(ena_return_fields("read_run"))
    fields = [x for x in wanted if x in available]
    url = f"{ENA_BASE}/filereport?" + urllib.parse.urlencode({
        "accession": project,
        "result": "read_run",
        "fields": ",".join(fields),
        "format": "tsv",
        "download": "true"
    })
    text = http_get(url, timeout=180)
    save_text(outdir / "01_ENA_read_run_metadata.tsv", text)
    return max(0, len(text.splitlines()) - 1)


def ena_sample_metadata(project: str, outdir: Path):
    wanted = [
        "accession","secondary_sample_accession","sample_alias",
        "sample_title","sample_description","collection_date",
        "age","sex","host","host_scientific_name","host_sex",
        "host_status","host_body_site","host_phenotype",
        "disease","isolation_source","tissue_type","country"
    ]
    available = set(ena_return_fields("sample"))
    fields = [x for x in wanted if x in available]

    # ENA studies can be addressed by primary or secondary project accession.
    queries = [
        f'study_accession="{project}"',
        f'secondary_study_accession="{project}"'
    ]
    best = ""
    for query in queries:
        url = f"{ENA_BASE}/search?" + urllib.parse.urlencode({
            "result": "sample",
            "query": query,
            "fields": ",".join(fields),
            "format": "tsv",
            "limit": "0"
        })
        try:
            text = http_get(url, timeout=180)
            if len(text.splitlines()) > len(best.splitlines()):
                best = text
        except Exception:
            pass

    if best:
        save_text(outdir / "02_ENA_sample_metadata.tsv", best)
        return max(0, len(best.splitlines()) - 1)
    return 0


def ncbi_find_bioproject_uid(project: str):
    url = ncbi_url("esearch.fcgi", {
        "db": "bioproject",
        "term": f"{project}[Project Accession]",
        "retmode": "xml",
        "retmax": "10",
    })
    xml = http_get(url)
    sleep_ncbi()
    root = ET.fromstring(xml)
    ids = [x.text for x in root.findall(".//IdList/Id") if x.text]
    if ids:
        return ids[0]

    # Secondary accession fallback.
    url = ncbi_url("esearch.fcgi", {
        "db": "bioproject",
        "term": project,
        "retmode": "xml",
        "retmax": "10",
    })
    xml = http_get(url)
    sleep_ncbi()
    root = ET.fromstring(xml)
    ids = [x.text for x in root.findall(".//IdList/Id") if x.text]
    return ids[0] if ids else None


def ncbi_bioproject_and_biosamples(project: str, outdir: Path):
    uid = ncbi_find_bioproject_uid(project)
    if not uid:
        return {"bioproject_uid": "", "biosample_count": 0}

    url = ncbi_url("efetch.fcgi", {
        "db": "bioproject",
        "id": uid,
        "retmode": "xml"
    })
    xml = http_get(url, timeout=180)
    sleep_ncbi()
    save_text(outdir / "03_NCBI_BioProject.xml", xml)

    elink = ncbi_url("elink.fcgi", {
        "dbfrom": "bioproject",
        "db": "biosample",
        "id": uid,
        "retmode": "xml"
    })
    link_xml = http_get(elink, timeout=180)
    sleep_ncbi()
    save_text(outdir / "04_NCBI_BioProject_to_BioSample_links.xml", link_xml)

    root = ET.fromstring(link_xml)
    ids = []
    for linksetdb in root.findall(".//LinkSetDb"):
        dbto = linksetdb.findtext("DbTo")
        if dbto == "biosample":
            ids.extend(
                x.text for x in linksetdb.findall("./Link/Id")
                if x.text
            )

    # Keep order, remove duplicates.
    ids = list(dict.fromkeys(ids))
    if not ids:
        return {"bioproject_uid": uid, "biosample_count": 0}

    raw_xml_parts = []
    rows = []
    all_attr_names = set()

    for start in range(0, len(ids), 200):
        batch = ids[start:start+200]
        url = ncbi_url("efetch.fcgi", {
            "db": "biosample",
            "id": ",".join(batch),
            "retmode": "xml"
        })
        part = http_get(url, timeout=180)
        sleep_ncbi()
        raw_xml_parts.append(part)

        try:
            rroot = ET.fromstring(part)
        except ET.ParseError:
            continue

        for bs in rroot.findall(".//BioSample"):
            row = {
                "BioSample_accession": bs.attrib.get("accession", ""),
                "BioSample_id": bs.attrib.get("id", ""),
                "publication_date": bs.attrib.get("publication_date", ""),
                "last_update": bs.attrib.get("last_update", ""),
            }

            desc = bs.find("./Description")
            if desc is not None:
                row["Title"] = desc.findtext("Title", default="")
                org = desc.find("./Organism")
                if org is not None:
                    row["Organism"] = org.attrib.get("taxonomy_name", "")
                    row["Taxonomy_ID"] = org.attrib.get("taxonomy_id", "")

            for attr in bs.findall(".//Attributes/Attribute"):
                name = (
                    attr.attrib.get("harmonized_name")
                    or attr.attrib.get("attribute_name")
                    or "attribute"
                )
                name = re.sub(r"[^A-Za-z0-9_]+", "_", name).strip("_")
                if not name:
                    continue
                value = (attr.text or "").strip()
                if name in row and row[name] and value and row[name] != value:
                    row[name] = f"{row[name]} | {value}"
                else:
                    row[name] = value
                all_attr_names.add(name)

            rows.append(row)

    save_text(
        outdir / "05_NCBI_BioSample_raw_batches.xml.txt",
        "\n\n<!-- NEXT_BATCH -->\n\n".join(raw_xml_parts)
    )

    base_cols = [
        "BioSample_accession","BioSample_id","Title","Organism",
        "Taxonomy_ID","publication_date","last_update"
    ]
    cols = base_cols + sorted(all_attr_names - set(base_cols))
    csv_path = outdir / "06_NCBI_BioSample_attributes.csv"
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    return {"bioproject_uid": uid, "biosample_count": len(rows)}


def pmc_oa_package(pmcid: str, outdir: Path):
    if not pmcid:
        return "NO_PMCID"

    api = PMC_OA_API + "?" + urllib.parse.urlencode({"id": pmcid})
    xml = http_get(api, timeout=120)
    save_text(outdir / "PMC_OA_api_response.xml", xml)

    root = ET.fromstring(xml)
    links = root.findall(".//link")
    if not links:
        return "NO_OA_PACKAGE"

    preferred = None
    for link in links:
        if link.attrib.get("format") == "tgz":
            preferred = link.attrib.get("href")
            break
    if not preferred:
        for link in links:
            href = link.attrib.get("href")
            if href:
                preferred = href
                break

    if not preferred:
        return "NO_DOWNLOAD_LINK"

    if preferred.startswith("ftp://ftp.ncbi.nlm.nih.gov"):
        preferred = preferred.replace(
            "ftp://ftp.ncbi.nlm.nih.gov",
            "https://ftp.ncbi.nlm.nih.gov",
            1
        )

    archive_dir = outdir / "08_PMC_open_access_package"
    archive_dir.mkdir(parents=True, exist_ok=True)
    suffix = ".tar.gz" if (".tar.gz" in preferred or ".tgz" in preferred) else Path(urllib.parse.urlparse(preferred).path).suffix
    archive = archive_dir / f"{pmcid}_oa_package{suffix or '.bin'}"

    download_file(preferred, archive, timeout=600)

    extract_dir = archive_dir / "extracted"
    extract_dir.mkdir(parents=True, exist_ok=True)

    try:
        if tarfile.is_tarfile(archive):
            with tarfile.open(archive, "r:*") as tf:
                tf.extractall(extract_dir)
            return "DOWNLOADED_AND_EXTRACTED"
        if zipfile.is_zipfile(archive):
            with zipfile.ZipFile(archive) as zf:
                zf.extractall(extract_dir)
            return "DOWNLOADED_AND_EXTRACTED"
    except Exception as e:
        return f"DOWNLOADED_EXTRACTION_FAILED: {e}"

    return "DOWNLOADED_NOT_ARCHIVE"


def github_snapshot(repo: str, outdir: Path):
    if not repo:
        return "NO_GITHUB"

    api = f"https://api.github.com/repos/{repo}"
    try:
        info = json.loads(http_get(api))
        branch = info.get("default_branch") or "main"
    except Exception:
        branch = "main"

    url = f"https://codeload.github.com/{repo}/zip/refs/heads/{branch}"
    dest = outdir / "09_public_analysis_repository"
    dest.mkdir(parents=True, exist_ok=True)
    zpath = dest / (repo.replace("/", "__") + ".zip")

    try:
        download_file(url, zpath, timeout=600)
    except Exception:
        if branch != "master":
            branch = "master"
            url = f"https://codeload.github.com/{repo}/zip/refs/heads/{branch}"
            download_file(url, zpath, timeout=600)
        else:
            raise

    extract = dest / "extracted"
    extract.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zpath) as zf:
        zf.extractall(extract)
    return f"DOWNLOADED:{branch}"


def download_direct_files(project: str, outdir: Path):
    statuses = []
    for folder, name, url in DIRECT_FILES.get(project, []):
        d = outdir / folder
        d.mkdir(parents=True, exist_ok=True)
        try:
            download_file(url, d / name, timeout=600)
            statuses.append(f"OK:{name}")
        except Exception as e:
            statuses.append(f"FAILED:{name}:{e}")
    return " ; ".join(statuses)


def write_source_registry(project: str, cfg: dict, outdir: Path):
    lines = [
        f"Project: {project}",
        f"Paper title: {cfg.get('paper_title','')}",
        f"DOI: {cfg.get('doi','')}",
        f"PMCID: {cfg.get('pmcid','')}",
        f"GitHub: {cfg.get('github','')}",
        f"Notes: {cfg.get('notes','')}",
        "",
        "This file is a source registry only. Do not treat inferred sample naming rules as verified clinical metadata."
    ]
    save_text(outdir / "00_source_registry.txt", "\n".join(lines))


def main():
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_rows = []

    print("Sepsis V2 metadata collector")
    print("Data root:", DATA_ROOT)
    print()

    for project, cfg in PROJECTS.items():
        project_dir = DATA_ROOT / project
        meta_dir = project_dir / "00_metadata" / "auto_fetched"
        meta_dir.mkdir(parents=True, exist_ok=True)

        print("=" * 70)
        print("PROJECT:", project)
        print("SAVE TO:", meta_dir)

        write_source_registry(project, cfg, meta_dir)

        row = {
            "Project": project,
            "ProjectDir": str(project_dir),
            "MetadataDir": str(meta_dir),
            "PaperTitle": cfg.get("paper_title",""),
            "DOI": cfg.get("doi",""),
            "PMCID": cfg.get("pmcid",""),
            "GitHub": cfg.get("github",""),
            "ENA_Run_Rows": "",
            "ENA_Sample_Rows": "",
            "NCBI_BioProject_UID": "",
            "NCBI_BioSample_Rows": "",
            "PMC_Package": "",
            "GitHub_Snapshot": "",
            "Direct_Public_Files": "",
            "Errors": "",
        }

        errors = []

        try:
            row["ENA_Run_Rows"] = ena_read_run(project, meta_dir)
            print("  ENA read_run:", row["ENA_Run_Rows"])
        except Exception as e:
            errors.append(f"ENA read_run: {type(e).__name__}: {e}")

        try:
            row["ENA_Sample_Rows"] = ena_sample_metadata(project, meta_dir)
            print("  ENA samples:", row["ENA_Sample_Rows"])
        except Exception as e:
            errors.append(f"ENA sample: {type(e).__name__}: {e}")

        try:
            ncbi = ncbi_bioproject_and_biosamples(project, meta_dir)
            row["NCBI_BioProject_UID"] = ncbi["bioproject_uid"]
            row["NCBI_BioSample_Rows"] = ncbi["biosample_count"]
            print("  NCBI BioSamples:", row["NCBI_BioSample_Rows"])
        except Exception as e:
            errors.append(f"NCBI: {type(e).__name__}: {e}")

        try:
            row["PMC_Package"] = pmc_oa_package(cfg.get("pmcid",""), meta_dir)
            print("  PMC:", row["PMC_Package"])
        except Exception as e:
            errors.append(f"PMC: {type(e).__name__}: {e}")

        try:
            row["GitHub_Snapshot"] = github_snapshot(cfg.get("github",""), meta_dir)
            print("  GitHub:", row["GitHub_Snapshot"])
        except Exception as e:
            errors.append(f"GitHub: {type(e).__name__}: {e}")

        try:
            row["Direct_Public_Files"] = download_direct_files(project, meta_dir)
            if row["Direct_Public_Files"]:
                print("  Direct public files:", row["Direct_Public_Files"])
        except Exception as e:
            errors.append(f"Direct files: {type(e).__name__}: {e}")

        row["Errors"] = " || ".join(errors)
        report_rows.append(row)

    report_path = RESULT_ROOT / f"V2_public_metadata_collection_report_{stamp}.csv"
    fields = list(report_rows[0].keys())
    with open(report_path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(report_rows)

    print()
    print("=" * 70)
    print("DONE")
    print("Collection report:")
    print(report_path)
    print("=" * 70)


if __name__ == "__main__":
    main()
