$ErrorActionPreference = "Stop"

$Root = "E:\sepsis_project"
$Script = Join-Path $Root "code\03_data_processing\91B4_infection_source_freeze_and_analysis_ATTRITION_FIXED.R"

if (-not (Test-Path $Script)) {
    Write-Host "ERROR: Step91B4 R script not found:"
    Write-Host $Script
    exit 1
}

$cmd = Get-Command Rscript.exe -ErrorAction SilentlyContinue

if ($null -ne $cmd) {
    Write-Host "Using Rscript:"
    Write-Host $cmd.Source
    & $cmd.Source $Script
    exit $LASTEXITCODE
}

$candidates = @()

if (Test-Path "C:\Program Files\R") {
    $candidates += Get-ChildItem "C:\Program Files\R" -Recurse -Filter Rscript.exe -ErrorAction SilentlyContinue |
        Where-Object { $_.FullName -notmatch "\\i386\\" }
}

$userRoots = @(
    "$env:LOCALAPPDATA\Programs\R",
    "$env:LOCALAPPDATA\R",
    "$env:USERPROFILE\AppData\Local\Programs\R"
)

foreach ($r in $userRoots) {
    if (Test-Path $r) {
        $candidates += Get-ChildItem $r -Recurse -Filter Rscript.exe -ErrorAction SilentlyContinue |
            Where-Object { $_.FullName -notmatch "\\i386\\" }
    }
}

$selected = $candidates |
    Sort-Object FullName -Descending |
    Select-Object -First 1

if ($null -eq $selected) {
    Write-Host "ERROR: Rscript.exe could not be located automatically."
    Write-Host "Run in RStudio:"
    Write-Host 'source("E:/sepsis_project/code/03_data_processing/91B4_infection_source_freeze_and_analysis_ATTRITION_FIXED.R", echo=TRUE)'
    exit 1
}

Write-Host "Using Rscript:"
Write-Host $selected.FullName

& $selected.FullName $Script
exit $LASTEXITCODE
