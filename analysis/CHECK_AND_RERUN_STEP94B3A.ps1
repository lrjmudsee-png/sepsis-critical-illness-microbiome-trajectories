$ErrorActionPreference = "Stop"

$Root = "E:\sepsis_project"
$CodeDir = Join-Path $Root "code\03_data_processing"
$Script = Join-Path $CodeDir "94B3A_CRA002354_exploratory_genus_source_analysis.R"
$Rscript = "C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Out = Join-Path $Root "results\V2_34B3A_CRA002354_EXPLORATORY_GENUS_SOURCE_ANALYSIS"
$Diag = Join-Path $Root "results\V2_34B3A_DIAGNOSTIC_AND_RETRY"

New-Item -ItemType Directory -Force -Path $Diag | Out-Null
New-Item -ItemType Directory -Force -Path $Out | Out-Null

$complete = Join-Path $Out "_STEP94B3A_COMPLETE.ok"

# 1. Pre-run status
$filesBefore = Get-ChildItem $Out -Force -ErrorAction SilentlyContinue |
    Select-Object Name,Length,LastWriteTime

$filesBefore |
    Export-Csv -NoTypeInformation -Encoding UTF8 (Join-Path $Diag "01_B3A_FILES_BEFORE.csv")

[pscustomobject]@{
    script_exists = Test-Path $Script
    rscript_exists = Test-Path $Rscript
    output_dir_exists = Test-Path $Out
    complete_flag_exists = Test-Path $complete
} | Export-Csv -NoTypeInformation -Encoding UTF8 (Join-Path $Diag "02_B3A_PRECHECK.csv")

if (Test-Path $complete) {
    "B3A already complete. No rerun performed." |
        Out-File -Encoding UTF8 (Join-Path $Diag "03_B3A_DECISION.txt")
    Write-Host "B3A already complete; no rerun."
    exit 0
}

if (-not (Test-Path $Script)) {
    "ERROR: B3A R script missing: $Script" |
        Out-File -Encoding UTF8 (Join-Path $Diag "03_B3A_DECISION.txt")
    throw "B3A R script is missing."
}

if (-not (Test-Path $Rscript)) {
    "ERROR: Rscript missing: $Rscript" |
        Out-File -Encoding UTF8 (Join-Path $Diag "03_B3A_DECISION.txt")
    throw "Rscript executable is missing."
}

"Rerunning B3A with stdout/stderr capture." |
    Out-File -Encoding UTF8 (Join-Path $Diag "03_B3A_DECISION.txt")

$stdout = Join-Path $Diag "04_B3A_STDOUT.txt"
$stderr = Join-Path $Diag "05_B3A_STDERR.txt"

$p = Start-Process `
    -FilePath $Rscript `
    -ArgumentList "`"$Script`"" `
    -NoNewWindow `
    -Wait `
    -PassThru `
    -RedirectStandardOutput $stdout `
    -RedirectStandardError $stderr

[pscustomobject]@{
    exit_code = $p.ExitCode
    completed_flag_after = Test-Path $complete
} | Export-Csv -NoTypeInformation -Encoding UTF8 (Join-Path $Diag "06_B3A_RUN_STATUS.csv")

# 2. Post-run output inventory
$filesAfter = Get-ChildItem $Out -Force -ErrorAction SilentlyContinue |
    Select-Object Name,Length,LastWriteTime

$filesAfter |
    Export-Csv -NoTypeInformation -Encoding UTF8 (Join-Path $Diag "07_B3A_FILES_AFTER.csv")

# 3. Copy compact copies of result summaries if they exist
$targets = @(
    "02_RAREFIED_GENUS_COVERAGE_SUMMARY.csv",
    "04_EXPLORATORY_BASELINE_GENUS_SOURCE_CONTRAST.csv",
    "06_EXPLORATORY_LONGITUDINAL_GENUS_SOURCE_TIME_INTERACTIONS.csv",
    "07_SENSITIVITY_EXCLUDE_OTHER_UNKNOWN_GENUS_INTERACTIONS.csv",
    "08_PRIMARY_LONGITUDINAL_HIT_ROBUSTNESS.csv",
    "09_DESCRIPTIVE_PULMONARY_VS_ABDOMINAL_GI_BASELINE_GENUS.csv",
    "12_GENUS_SOURCE_EVIDENCE_SUMMARY.csv",
    "13_STEP94B3A_INTERPRETATION.txt",
    "_STEP94B3A_COMPLETE.ok"
)

foreach ($f in $targets) {
    $src = Join-Path $Out $f
    if (Test-Path $src) {
        Copy-Item $src (Join-Path $Diag ("RESULT_" + $f)) -Force
    }
}

# 4. Final decision
if ($p.ExitCode -eq 0 -and (Test-Path $complete)) {
    "B3A_COMPLETE_SUCCESSFULLY" |
        Out-File -Encoding UTF8 (Join-Path $Diag "08_FINAL_DECISION.txt")
    Write-Host "B3A COMPLETE SUCCESSFULLY"
} else {
    "B3A_FAILED_OR_INCOMPLETE; inspect 05_B3A_STDERR.txt" |
        Out-File -Encoding UTF8 (Join-Path $Diag "08_FINAL_DECISION.txt")
    Write-Host "B3A FAILED OR INCOMPLETE"
    Write-Host "Inspect:" $stderr
    exit 1
}
