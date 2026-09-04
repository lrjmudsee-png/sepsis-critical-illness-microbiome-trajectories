$ErrorActionPreference = "Stop"

$Root = "E:\sepsis_project"
$Script = Join-Path $Root "code\03_data_processing\88A_longitudinal_diversity_and_displacement.R"

if (-not (Test-Path $Script)) {
    Write-Host "ERROR: Step88A R script not found:"
    Write-Host $Script
    exit 1
}

# 1) Try PATH first
$cmd = Get-Command Rscript.exe -ErrorAction SilentlyContinue

if ($null -ne $cmd) {
    Write-Host "Using Rscript:"
    Write-Host $cmd.Source
    & $cmd.Source $Script
    exit $LASTEXITCODE
}

# 2) Search common machine-wide installation folder
$candidates = @()

if (Test-Path "C:\Program Files\R") {
    $candidates += Get-ChildItem "C:\Program Files\R" -Recurse -Filter Rscript.exe -ErrorAction SilentlyContinue |
        Where-Object { $_.FullName -notmatch "\\i386\\" }
}

# 3) Search common per-user installation folders
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
    Write-Host ""
    Write-Host "If RStudio works, run this inside RStudio:"
    Write-Host 'source("E:/sepsis_project/code/03_data_processing/88A_longitudinal_diversity_and_displacement.R", echo=TRUE)'
    exit 1
}

Write-Host "Using Rscript:"
Write-Host $selected.FullName

& $selected.FullName $Script
exit $LASTEXITCODE
