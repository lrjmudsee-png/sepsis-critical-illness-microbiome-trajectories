$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '06_018_cra002354_genus_coverage_robustness.R')

& $Rscript $Script
exit $LASTEXITCODE
