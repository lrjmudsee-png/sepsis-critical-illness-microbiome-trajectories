# -*- coding: utf-8 -*-
r"""
Sepsis V2 - Step 02
Profile metadata files and infer candidate columns for longitudinal analysis.

This script:
1. Finds the latest Step 01 inventory outputs automatically.
2. Reads all metadata candidate files.
3. Profiles columns and example values.
4. Infers likely roles:
   patient_id, sample_id, run_id, time, sepsis_group,
   infection_source, outcome, age, sex, antibiotics, severity.
5. Chooses a suggested primary metadata source per project.
6. Does NOT modify any raw or processed data.

Run on Windows:
python "E:\sepsis_project\code\03_data_processing\61_profile_v2_metadata.py"
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from datetime import datetime

import pandas as pd


INVENTORY_DIR = Path(r"E:\sepsis_project\results\V2_00_inventory")
OUTPUT_DIR = Path(r"E:\sepsis_project\results\V2_01_metadata_profile")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


ROLE_PATTERNS = {
    "patient_id": [
        r"patient", r"subject", r"participant", r"individual",
        r"host.?id", r"donor", r"person.?id"
    ],
    "sample_id": [
        r"sample.?id", r"sample_accession", r"biosample",
        r"specimen.?id"
    ],
    "run_id": [
        r"run.?id", r"run_accession", r"accession",
        r"^run$", r"sra.?run"
    ],
    "time": [
        r"time", r"timepoint", r"time_point", r"day",
        r"visit", r"collection.?date", r"sampling.?date",
        r"collection.?time", r"sampling.?time",
        r"icu.?day", r"hospital.?day"
    ],
    "sepsis_group": [
        r"sepsis", r"diagnos", r"disease", r"group",
        r"case.?control", r"clinical.?group", r"phenotype"
    ],
    "infection_source": [
        r"infection.?source", r"infection.?site",
        r"site.?of.?infection", r"primary.?infection",
        r"focus.?of.?infection"
    ],
    "outcome": [
        r"outcome", r"mortality", r"death", r"dead",
        r"survival", r"survivor", r"deceased",
        r"discharge.?status"
    ],
    "age": [r"^age$", r"age.?year", r"host.?age"],
    "sex": [r"^sex$", r"gender", r"host.?sex"],
    "antibiotics": [
        r"antibiotic", r"antimicrobial", r"abx"
    ],
    "severity": [
        r"sofa", r"apache", r"severity", r"septic.?shock",
        r"shock"
    ],
}


def latest_file(pattern: str) -> Path:
    files = sorted(
        INVENTORY_DIR.glob(pattern),
        key=lambda p: p.stat().st_mtime,
        reverse=True
    )
    if not files:
        raise FileNotFoundError(
            f"No file matching {pattern} under {INVENTORY_DIR}"
        )
    return files[0]


def normalize_col(x: str) -> str:
    x = str(x).strip().lower()
    x = re.sub(r"[\s\-/\\\.\(\)\[\]:]+", "_", x)
    x = re.sub(r"_+", "_", x)
    return x.strip("_")


def infer_roles(column_name: str) -> list[str]:
    n = normalize_col(column_name)
    roles = []

    for role, patterns in ROLE_PATTERNS.items():
        for pat in patterns:
            if re.search(pat, n, flags=re.I):
                roles.append(role)
                break

    return roles


def read_table(path: Path) -> tuple[pd.DataFrame, str]:
    ext = path.suffix.lower()

    if ext == ".csv":
        return pd.read_csv(path, dtype=str, low_memory=False), "csv"

    if ext in {".tsv", ".txt"}:
        # First try tab, then separator inference.
        try:
            df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
            if df.shape[1] > 1:
                return df, "tsv"
        except Exception:
            pass

        return (
            pd.read_csv(
                path,
                sep=None,
                engine="python",
                dtype=str,
                on_bad_lines="skip"
            ),
            "text_auto"
        )

    if ext in {".xlsx", ".xls"}:
        xl = pd.ExcelFile(path)
        frames = []

        for sheet in xl.sheet_names:
            temp = pd.read_excel(path, sheet_name=sheet, dtype=str)
            temp.insert(0, "__sheet__", sheet)
            frames.append(temp)

        if not frames:
            return pd.DataFrame(), "excel"

        # Keep all sheets for profiling. Columns are unioned automatically.
        return pd.concat(frames, ignore_index=True, sort=False), "excel"

    raise ValueError(f"Unsupported metadata format: {ext}")


def sample_values(series: pd.Series, n: int = 5) -> str:
    vals = (
        series.dropna()
        .astype(str)
        .map(str.strip)
    )
    vals = vals[vals != ""].drop_duplicates().head(n).tolist()
    return " | ".join(vals)


def filename_score(filename: str) -> int:
    n = filename.lower()
    score = 0

    if "metadata_final" in n:
        score += 120
    elif "metadata_after_taxonomy_cleaning" in n:
        score += 100
    elif "clinical" in n:
        score += 90
    elif "metadata" in n:
        score += 70

    if "run_manifest" in n or "ena_run_manifest" in n:
        score += 35

    if "human_only" in n or "gut_clear" in n or "all_gut" in n:
        score += 10

    if "sample_read_depth" in n:
        score -= 80

    return score


def role_bonus(role_set: set[str]) -> int:
    score = 0
    weights = {
        "patient_id": 50,
        "time": 50,
        "sample_id": 25,
        "run_id": 20,
        "sepsis_group": 25,
        "infection_source": 20,
        "outcome": 20,
        "age": 5,
        "sex": 5,
        "antibiotics": 10,
        "severity": 10,
    }
    for role, w in weights.items():
        if role in role_set:
            score += w
    return score


def main():
    inventory_path = latest_file("V2_project_inventory_*.csv")
    candidates_path = latest_file("V2_metadata_candidates_*.csv")

    print("Using inventory:")
    print(inventory_path)
    print("Using metadata candidates:")
    print(candidates_path)
    print()

    inventory = pd.read_csv(inventory_path, dtype=str)
    candidates = pd.read_csv(candidates_path, dtype=str)

    all_projects = (
        inventory["Project"]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .tolist()
    )

    file_rows = []
    column_rows = []

    for _, rec in candidates.iterrows():
        project = str(rec["Project"])
        raw_path = str(rec["FullPath"])
        path = Path(raw_path)

        print(f"Profiling: {project} | {path.name}")

        base_file_row = {
            "Project": project,
            "FileName": path.name,
            "FullPath": raw_path,
            "Exists": path.exists(),
            "LoadStatus": "",
            "Rows": None,
            "Columns": None,
            "DetectedRoles": "",
            "FilenameScore": filename_score(path.name),
            "RoleBonus": 0,
            "TotalScore": None,
            "SuggestedPrimary": "NO",
            "Error": "",
        }

        if not path.exists():
            base_file_row["LoadStatus"] = "MISSING_FILE"
            base_file_row["Error"] = "Path does not exist"
            file_rows.append(base_file_row)
            continue

        try:
            df, file_type = read_table(path)
            df.columns = [str(c) for c in df.columns]

            detected_roles = set()

            for col in df.columns:
                roles = infer_roles(col)
                detected_roles.update(roles)

                s = df[col]
                nonnull = int(s.notna().sum())
                unique = int(s.dropna().astype(str).nunique())

                column_rows.append({
                    "Project": project,
                    "FileName": path.name,
                    "FullPath": raw_path,
                    "Column": col,
                    "NormalizedColumn": normalize_col(col),
                    "InferredRoles": ";".join(roles),
                    "NonNull": nonnull,
                    "Unique": unique,
                    "ExampleValues": sample_values(s, 5),
                })

            bonus = role_bonus(detected_roles)
            total = filename_score(path.name) + bonus

            base_file_row.update({
                "LoadStatus": "OK",
                "Rows": len(df),
                "Columns": len(df.columns),
                "DetectedRoles": ";".join(sorted(detected_roles)),
                "RoleBonus": bonus,
                "TotalScore": total,
                "Error": "",
            })

            file_rows.append(base_file_row)

        except Exception as e:
            base_file_row["LoadStatus"] = "ERROR"
            base_file_row["Error"] = repr(e)
            file_rows.append(base_file_row)

    file_profile = pd.DataFrame(file_rows)
    column_profile = pd.DataFrame(column_rows)

    # Choose one suggested primary metadata file per project.
    if not file_profile.empty:
        for project in file_profile["Project"].dropna().unique():
            mask = (
                (file_profile["Project"] == project)
                & (file_profile["LoadStatus"] == "OK")
            )
            sub = file_profile.loc[mask].copy()

            if len(sub):
                sub["TotalScoreNumeric"] = pd.to_numeric(
                    sub["TotalScore"],
                    errors="coerce"
                ).fillna(-999999)

                best_idx = sub.sort_values(
                    ["TotalScoreNumeric", "Rows"],
                    ascending=[False, False]
                ).index[0]

                file_profile.loc[best_idx, "SuggestedPrimary"] = "YES"

    # Build project-level readiness.
    project_rows = []

    for project in all_projects:
        files = file_profile[
            (file_profile["Project"] == project)
            & (file_profile["LoadStatus"] == "OK")
        ].copy()

        primary = files[files["SuggestedPrimary"] == "YES"]

        if len(primary):
            primary_row = primary.iloc[0]
            primary_path = primary_row["FullPath"]
            primary_name = primary_row["FileName"]
        else:
            primary_path = ""
            primary_name = ""

        cols = column_profile[column_profile["Project"] == project]

        role_to_columns = {}
        for role in ROLE_PATTERNS:
            matched = cols[
                cols["InferredRoles"]
                .fillna("")
                .str.split(";")
                .apply(lambda xs: role in xs)
            ]
            role_to_columns[role] = ";".join(
                matched["Column"].drop_duplicates().astype(str).tolist()
            )

        has_patient = bool(role_to_columns["patient_id"])
        has_time = bool(role_to_columns["time"])
        has_sample_or_run = bool(
            role_to_columns["sample_id"] or role_to_columns["run_id"]
        )

        if len(files) == 0:
            readiness = "NO_METADATA"
        elif has_patient and has_time and has_sample_or_run:
            readiness = "READY_FOR_MANUAL_MAPPING"
        else:
            readiness = "PARTIAL_METADATA"

        project_rows.append({
            "Project": project,
            "MetadataFilesReadable": len(files),
            "SuggestedPrimaryFile": primary_name,
            "SuggestedPrimaryPath": primary_path,
            "PatientColumns": role_to_columns["patient_id"],
            "SampleColumns": role_to_columns["sample_id"],
            "RunColumns": role_to_columns["run_id"],
            "TimeColumns": role_to_columns["time"],
            "SepsisGroupColumns": role_to_columns["sepsis_group"],
            "InfectionSourceColumns": role_to_columns["infection_source"],
            "OutcomeColumns": role_to_columns["outcome"],
            "AgeColumns": role_to_columns["age"],
            "SexColumns": role_to_columns["sex"],
            "AntibioticColumns": role_to_columns["antibiotics"],
            "SeverityColumns": role_to_columns["severity"],
            "Readiness": readiness,
        })

    project_readiness = pd.DataFrame(project_rows)

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")

    file_out = OUTPUT_DIR / f"V2_metadata_file_profile_{ts}.csv"
    col_out = OUTPUT_DIR / f"V2_metadata_column_profile_{ts}.csv"
    project_out = OUTPUT_DIR / f"V2_project_metadata_readiness_{ts}.csv"

    file_profile.to_csv(file_out, index=False, encoding="utf-8-sig")
    column_profile.to_csv(col_out, index=False, encoding="utf-8-sig")
    project_readiness.to_csv(
        project_out,
        index=False,
        encoding="utf-8-sig"
    )

    print()
    print("================ STEP 02 COMPLETE ================")
    print(project_readiness[
        [
            "Project",
            "MetadataFilesReadable",
            "SuggestedPrimaryFile",
            "Readiness",
        ]
    ].to_string(index=False))

    print()
    print("File profile:")
    print(file_out)
    print("Column profile:")
    print(col_out)
    print("Project readiness:")
    print(project_out)
    print("==================================================")


if __name__ == "__main__":
    main()
