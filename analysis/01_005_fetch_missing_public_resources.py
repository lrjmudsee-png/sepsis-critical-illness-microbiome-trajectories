# -*- coding: utf-8 -*-
from __future__ import annotations

import argparse
import html
import json
import re
import shutil
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
from datetime import datetime
import csv

DATA_ROOT = Path(r"E:\sepsis_project\data")
RESULT_ROOT = Path(r"E:\sepsis_project\results\V2_04_metadata_enrichment")
RESULT_ROOT.mkdir(parents=True, exist_ok=True)
UA = "Mozilla/5.0 SepsisV2MetadataEnrichment/1.0"

SPRINGER_PRJEB82425 = (
    "https://media.springernature.com/original/"
    "springer-static/esm/art%3A10.1186%2Fs12879-025-10825-6/"
    "MediaObjects/12879_2025_10825_MOESM1_ESM.docx"
)
GITHUB_PRJNA516701_REPO = "MicrobiomeALIR/ICUgutMbioMethods"
FRONTIERS_PRJNA578267_HTML = (
    "https://www.frontiersin.org/journals/"
    "cellular-and-infection-microbiology/articles/"
    "10.3389/fcimb.2019.00467/full"
)
FRONTIERS_PRJNA578267_XML = (
    "https://www.frontiersin.org/journals/"
    "cellular-and-infection-microbiology/articles/"
    "10.3389/fcimb.2019.00467/xml"
)

def request_bytes(url: str, timeout: int = 300) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "*/*"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()

def request_text(url: str, timeout: int = 180) -> str:
    return request_bytes(url, timeout).decode("utf-8", errors="replace")

def download(url: str, path: Path, force: bool = False, timeout: int = 600):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0 and not force:
        return "EXISTS"
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as r, open(path, "wb") as f:
        shutil.copyfileobj(r, f)
    return "DOWNLOADED"

def extract_zip(path: Path, outdir: Path, force: bool = False):
    if not path.exists():
        return "ZIP_MISSING"
    outdir.mkdir(parents=True, exist_ok=True)
    marker = outdir / ".extract_complete"
    if marker.exists() and not force:
        return "EXISTS"
    with zipfile.ZipFile(path) as zf:
        zf.extractall(outdir)
    marker.write_text("ok", encoding="utf-8")
    return "EXTRACTED"

def github_snapshot(repo: str, outdir: Path, force: bool = False):
    api = f"https://api.github.com/repos/{repo}"
    try:
        info = json.loads(request_text(api))
        branch = info.get("default_branch") or "master"
    except Exception:
        branch = "master"

    target = outdir / "09_public_analysis_repository"
    target.mkdir(parents=True, exist_ok=True)
    zip_path = target / (repo.replace("/", "__") + ".zip")
    url = f"https://codeload.github.com/{repo}/zip/refs/heads/{branch}"

    try:
        status = download(url, zip_path, force=force)
    except Exception:
        branch = "master" if branch != "master" else "main"
        url = f"https://codeload.github.com/{repo}/zip/refs/heads/{branch}"
        status = download(url, zip_path, force=force)

    extract_status = extract_zip(zip_path, target / "extracted", force=force)
    return f"{status};{extract_status};branch={branch}"

def normalize_found_url(url: str, base: str) -> str:
    url = html.unescape(url).replace("\\/", "/")
    return urllib.parse.urljoin(base, url)

def find_supplement_urls(text: str, base_url: str) -> list[str]:
    candidates = []

    pattern1 = r'(?:href|src)\s*=\s*["\']([^"\']+)["\']'
    for match in re.findall(pattern1, text, flags=re.I):
        low = match.lower()
        if (
            "supplement" in low
            or "supplementary" in low
            or "mediaobjects" in low
            or re.search(r"\.(xlsx?|docx?|pdf|zip|csv|tsv)(?:\?|$)", low)
        ):
            candidates.append(normalize_found_url(match, base_url))

    pattern2 = r'(?:xlink:href|href)\s*=\s*["\']([^"\']+)["\']'
    for match in re.findall(pattern2, text, flags=re.I):
        low = match.lower()
        if (
            "supplement" in low
            or "supplementary" in low
            or re.search(r"\.(xlsx?|docx?|pdf|zip|csv|tsv)(?:\?|$)", low)
        ):
            candidates.append(normalize_found_url(match, base_url))

    for match in re.findall(r'https?://[^\s"<>]+', text):
        low = match.lower()
        if (
            "supplement" in low
            or "supplementary" in low
            or re.search(r"\.(xlsx?|docx?|pdf|zip|csv|tsv)(?:\?|$)", low)
        ):
            candidates.append(html.unescape(match.rstrip(").,;'")))

    out = []
    for u in candidates:
        if "#supplementary-material" in u and u.split("#")[0] == base_url.split("#")[0]:
            continue
        if u not in out:
            out.append(u)
    return out

def safe_filename_from_url(url: str, idx: int) -> str:
    parsed = urllib.parse.urlparse(url)
    name = Path(urllib.parse.unquote(parsed.path)).name
    if not name or "." not in name:
        name = f"supplement_candidate_{idx:02d}.bin"
    return re.sub(r'[<>:"/\\|?*]+', "_", name)

def fetch_frontiers_supplements(outdir: Path, force: bool = False):
    folder = outdir / "07_known_public_supplement"
    folder.mkdir(parents=True, exist_ok=True)

    sources = []
    errors = []

    for kind, url in [("html", FRONTIERS_PRJNA578267_HTML), ("xml", FRONTIERS_PRJNA578267_XML)]:
        try:
            text = request_text(url)
            (folder / f"PRJNA578267_article_source.{kind}.txt").write_text(text, encoding="utf-8")
            sources.extend(find_supplement_urls(text, url))
        except Exception as e:
            errors.append(f"{kind}:{type(e).__name__}:{e}")

    sources = list(dict.fromkeys(sources))
    (folder / "PRJNA578267_discovered_supplement_urls.txt").write_text(
        "\n".join(sources), encoding="utf-8"
    )

    statuses = []
    for idx, url in enumerate(sources, start=1):
        low = url.lower()
        if not (
            "supplement" in low
            or "supplementary" in low
            or re.search(r"\.(xlsx?|docx?|zip|csv|tsv)(?:\?|$)", low)
        ):
            continue

        name = safe_filename_from_url(url, idx)
        target = folder / name
        try:
            status = download(url, target, force=force)
            head = target.read_bytes()[:300].lower()
            if b"<html" in head or b"<!doctype html" in head:
                target.unlink(missing_ok=True)
                statuses.append(f"HTML_PAGE_SKIPPED:{url}")
            else:
                statuses.append(f"{status}:{target.name}")
        except Exception as e:
            statuses.append(f"FAILED:{url}:{type(e).__name__}:{e}")

    return " | ".join(statuses + errors)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    rows = []

    p = "PRJEB82425"
    out = DATA_ROOT / p / "00_metadata" / "auto_fetched"
    supp = out / "07_known_public_supplement"
    supp.mkdir(parents=True, exist_ok=True)
    target = supp / "PRJEB82425_Supplementary_Material_1.docx"
    try:
        status = download(SPRINGER_PRJEB82425, target, force=args.force)
    except Exception as e:
        status = f"FAILED:{type(e).__name__}:{e}"
    rows.append({
        "Project": p,
        "Resource": "Springer Supplementary Material 1",
        "Status": status,
        "LocalPath": str(target),
        "Source": SPRINGER_PRJEB82425,
    })

    p = "PRJNA516701"
    out = DATA_ROOT / p / "00_metadata" / "auto_fetched"
    try:
        status = github_snapshot(GITHUB_PRJNA516701_REPO, out, force=args.force)
    except Exception as e:
        status = f"FAILED:{type(e).__name__}:{e}"
    rows.append({
        "Project": p,
        "Resource": "ICUgutMbioMethods GitHub snapshot",
        "Status": status,
        "LocalPath": str(out / "09_public_analysis_repository"),
        "Source": f"https://github.com/{GITHUB_PRJNA516701_REPO}",
    })

    p = "PRJNA578267"
    out = DATA_ROOT / p / "00_metadata" / "auto_fetched"
    try:
        status = fetch_frontiers_supplements(out, force=args.force)
        if not status:
            status = "NO_DIRECT_FILE_DISCOVERED"
    except Exception as e:
        status = f"FAILED:{type(e).__name__}:{e}"
    rows.append({
        "Project": p,
        "Resource": "Frontiers supplementary discovery/download",
        "Status": status,
        "LocalPath": str(out / "07_known_public_supplement"),
        "Source": FRONTIERS_PRJNA578267_HTML,
    })

    report = RESULT_ROOT / f"V2_metadata_enrichment_report_{stamp}.csv"
    with open(report, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["Project","Resource","Status","LocalPath","Source"])
        w.writeheader()
        w.writerows(rows)

    print("STEP 04 COMPLETE")
    for r in rows:
        print(r["Project"], "|", r["Resource"], "|", r["Status"])
    print("Report:", report)

if __name__ == "__main__":
    main()
