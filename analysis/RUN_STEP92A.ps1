$ErrorActionPreference="Stop"

$Root="E:\sepsis_project"

$Script=Join-Path `
$Root `
"code\03_data_processing\92A_EII_feasibility_audit.R"

if(-not(Test-Path $Script)){
    Write-Host "Missing script:"
    Write-Host $Script
    exit 1
}

$r=Get-Command Rscript.exe -ErrorAction SilentlyContinue

if($null -ne $r){
    & $r.Source $Script
    exit $LASTEXITCODE
}

Write-Host "Rscript not found."
Write-Host "Use RStudio:"
Write-Host 'source("E:/sepsis_project/code/03_data_processing/92A_EII_feasibility_audit.R")'
