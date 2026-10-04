$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '05_050_longitudinal_evidence_freeze.R')

& $Rscript $Script

exit $LASTEXITCODE
