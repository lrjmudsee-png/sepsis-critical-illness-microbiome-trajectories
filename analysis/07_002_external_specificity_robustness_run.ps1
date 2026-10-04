$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '07_002_external_specificity_robustness.R')

& $Rscript $Script

exit $LASTEXITCODE
