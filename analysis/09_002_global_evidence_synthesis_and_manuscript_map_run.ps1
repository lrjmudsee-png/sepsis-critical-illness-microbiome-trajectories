$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '09_002_global_evidence_synthesis_and_manuscript_map.R')

& $Rscript $Script
exit $LASTEXITCODE
