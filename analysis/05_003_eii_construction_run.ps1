$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_003_eii_construction.R')

if(!(Test-Path $Rscript)){
    Write-Host "Rscript missing:"
    Write-Host $Rscript
    exit 1
}

& $Rscript $Script
exit $LASTEXITCODE
