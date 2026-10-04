$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_005_eii_stability_comparison.R')

& $Rscript $Script

exit $LASTEXITCODE
