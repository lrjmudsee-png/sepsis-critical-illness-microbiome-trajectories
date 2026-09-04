$ErrorActionPreference = "Stop"

$Root = "E:\sepsis_project"
$Script = Join-Path $Root "code\87B_freeze_replicate_rule_and_build_analysis_objects.R"

if (-not (Test-Path $Script)) {
    Write-Host "ERROR: R script not found:"
    Write-Host $Script
    Write-Host ""
    Write-Host "Copy 87B_freeze_replicate_rule_and_build_analysis_objects.R into E:\sepsis_project\code\ first."
    exit 1
}

$Rscript = Get-Command Rscript.exe -ErrorAction SilentlyContinue

if ($null -ne $Rscript) {
    & $Rscript.Source $Script
    exit $LASTEXITCODE
}

$CommonR = @(
    "C:\Program Files\R\R-4.4.0\bin\Rscript.exe",
    "C:\Program Files\R\R-4.4.1\bin\Rscript.exe",
    "C:\Program Files\R\R-4.4.2\bin\Rscript.exe",
    "C:\Program Files\R\R-4.4.3\bin\Rscript.exe"
)

$Found = $CommonR | Where-Object { Test-Path $_ } | Select-Object -First 1

if ($null -eq $Found) {
    Write-Host "ERROR: Rscript.exe was not found."
    Write-Host "Run the R script from RStudio with:"
    Write-Host 'source("E:/sepsis_project/code/87B_freeze_replicate_rule_and_build_analysis_objects.R", echo=TRUE)'
    exit 1
}

& $Found $Script
exit $LASTEXITCODE
