$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_002_eii_baseline_harmonization_corrected.R')

& $Rscript $Script

exit $LASTEXITCODE
