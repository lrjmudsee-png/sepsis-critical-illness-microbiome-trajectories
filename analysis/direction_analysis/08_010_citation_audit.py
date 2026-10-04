#!/usr/bin/env python3
"""Read-only Crossref and in-text/reference consistency audit for the draft manuscript."""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
from difflib import SequenceMatcher
from pathlib import Path

from docx import Document


DOI_RE = re.compile(r"https?://doi\.org/([^\s]+?)(?:[.,;)]?$)", re.I)
REF_RE = re.compile(r"^(\d+)\.\s+(.*)$", re.S)
CITATION_RE = re.compile(r"\[([0-9,;\-–\s]+)\]")


def normalize(value: str) -> str:
    value = unicodedata.normalize("NFKD", value).lower()
    value = re.sub(r"[^a-z0-9]+", " ", value)
    return " ".join(value.split())


def title_from_reference(reference: str) -> str:
    without_doi = re.sub(r"\s*https?://doi\.org/\S+\.?$", "", reference, flags=re.I)
    parts = [part.strip() for part in without_doi.split(". ")]
    return parts[1] if len(parts) >= 3 else ""


def expand_numbers(content: str) -> set[int]:
    numbers: set[int] = set()
    for part in re.split(r"[,;]", content.replace("–", "-")):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            bounds = [item.strip() for item in part.split("-", 1)]
            if all(item.isdigit() for item in bounds):
                low, high = map(int, bounds)
                numbers.update(range(min(low, high), max(low, high) + 1))
        elif part.isdigit():
            numbers.add(int(part))
    return numbers


def crossref(doi: str) -> tuple[dict, str]:
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "sepsis-manuscript-citation-audit/1.0 (mailto:metadata-audit@example.invalid)"},
    )
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            payload = json.load(response)
        return payload.get("message", {}), ""
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, json.JSONDecodeError) as error:
        return {}, f"{type(error).__name__}: {error}"


def first(value):
    return value[0] if isinstance(value, list) and value else ""


def published_year(message: dict) -> str:
    for key in ("published-print", "published-online", "published", "issued"):
        parts = message.get(key, {}).get("date-parts", [])
        if parts and parts[0]:
            return str(parts[0][0])
    return ""


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if path.exists():
        raise RuntimeError(f"Refusing to overwrite: {path}")
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manuscript", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    manuscript = Path(args.manuscript).resolve()
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "citation_crossref_audit.csv"
    report_path = output / "citation_audit_report.md"
    if csv_path.exists() or report_path.exists():
        raise RuntimeError("Refusing to overwrite an existing citation audit")

    document = Document(manuscript)
    paragraphs = [" ".join(paragraph.text.split()) for paragraph in document.paragraphs]
    try:
        reference_start = next(index for index, text in enumerate(paragraphs) if text.strip().lower() == "references")
    except StopIteration as error:
        raise RuntimeError("References heading not found") from error

    cited: set[int] = set()
    for text in paragraphs[:reference_start]:
        for match in CITATION_RE.finditer(text):
            cited.update(expand_numbers(match.group(1)))

    references: list[tuple[int, str]] = []
    for text in paragraphs[reference_start + 1 :]:
        match = REF_RE.match(text)
        if match:
            references.append((int(match.group(1)), match.group(2)))
    if not references:
        raise RuntimeError("No numbered references found")

    rows: list[dict[str, object]] = []
    for number, reference in references:
        doi_match = DOI_RE.search(reference)
        doi = doi_match.group(1).rstrip(".,;") if doi_match else ""
        message, error = crossref(doi) if doi else ({}, "DOI absent")
        crossref_title = first(message.get("title", []))
        draft_title = title_from_reference(reference)
        title_similarity = SequenceMatcher(None, normalize(draft_title), normalize(crossref_title)).ratio() if crossref_title else 0.0
        year = published_year(message)
        year_match = bool(year and re.search(rf"\b{re.escape(year)}\b", reference))
        container = first(message.get("container-title", []))
        container_tokens = set(normalize(container).split())
        reference_tokens = set(normalize(reference).split())
        container_overlap = len(container_tokens & reference_tokens) / len(container_tokens) if container_tokens else 0.0
        update_relations = message.get("relation", {})
        relation_flag = "|".join(sorted(update_relations)) if update_relations else "NONE_REPORTED_BY_CROSSREF"
        if not doi or error or not crossref_title:
            status = "FAIL"
        elif title_similarity >= 0.90 and year_match and container_overlap >= 0.50:
            status = "PASS"
        elif title_similarity >= 0.75 and year_match:
            status = "REVIEW"
        else:
            status = "FAIL"
        rows.append(
            {
                "reference_number": number,
                "cited_in_text": number in cited,
                "doi": doi,
                "doi_url": f"https://doi.org/{doi}" if doi else "",
                "crossref_status": "FOUND" if message else "NOT_FOUND",
                "draft_title": draft_title,
                "crossref_title": crossref_title,
                "title_similarity": f"{title_similarity:.6f}",
                "crossref_container": container,
                "container_token_overlap": f"{container_overlap:.6f}",
                "crossref_year": year,
                "year_matches_draft": year_match,
                "crossref_relation_flag": relation_flag,
                "automated_status": status,
                "lookup_error": error,
                "reference_text": reference,
            }
        )
        time.sleep(0.10)

    write_csv(csv_path, rows)
    numbers = {number for number, _ in references}
    orphan_in_text = sorted(cited - numbers)
    orphan_reference = sorted(numbers - cited)
    counts = {status: sum(row["automated_status"] == status for row in rows) for status in ("PASS", "REVIEW", "FAIL")}
    manual = [row for row in rows if row["automated_status"] != "PASS"]
    relation_rows = [row for row in rows if row["crossref_relation_flag"] != "NONE_REPORTED_BY_CROSSREF"]
    manual_lines = "\n".join(
        f"- Reference {row['reference_number']}: {row['automated_status']}; {row['lookup_error'] or 'metadata mismatch requires inspection'}"
        for row in manual
    ) or "- None from DOI/title/year/container matching."
    report = f"""# Citation audit report

## Material Passport

- Origin Skill: academic-research-suite / citation-compliance
- Origin Mode: citation audit only
- Verification Status: MACHINE_METADATA_CHECK_COMPLETE; CLAIM_SUPPORT_AND_RETRACTION_REVIEW_REMAIN_HUMAN
- Version: citation-audit-v1

## Summary

- Numbered reference entries: {len(rows)}
- Distinct reference numbers cited in text: {len(cited)}
- Orphan in-text numbers: {', '.join(map(str, orphan_in_text)) if orphan_in_text else 'none'}
- Uncited reference entries: {', '.join(map(str, orphan_reference)) if orphan_reference else 'none'}
- Crossref metadata PASS: {counts['PASS']}
- Crossref metadata REVIEW: {counts['REVIEW']}
- Crossref metadata FAIL: {counts['FAIL']}
- Crossref relation/update fields present: {len(relation_rows)}

## Items requiring manual review after machine checking

{manual_lines}

## Interpretation boundary

The machine audit checks DOI registry resolution plus title, journal-token and publication-year consistency. It does not establish that each cited article supports the exact manuscript sentence, and absence of a Crossref relation flag is not proof that an article has never been retracted or corrected. Final claim-to-source reading and a dedicated retraction/expression-of-concern check remain author responsibilities.

No manuscript file was modified.
"""
    report_path.write_text(report, encoding="utf-8")
    print(f"Citation audit complete: {counts}; orphan_in_text={orphan_in_text}; orphan_reference={orphan_reference}")
    return 2 if counts["FAIL"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
