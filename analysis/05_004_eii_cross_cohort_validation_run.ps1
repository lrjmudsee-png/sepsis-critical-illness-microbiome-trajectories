$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_004_eii_cross_cohort_validation.R')

& $Rscript $Script

exit $LASTEXITCODE
