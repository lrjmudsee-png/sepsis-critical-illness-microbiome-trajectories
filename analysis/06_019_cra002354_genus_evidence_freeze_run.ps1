$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '06_019_cra002354_genus_evidence_freeze.R')

& $Rscript $Script
exit $LASTEXITCODE
