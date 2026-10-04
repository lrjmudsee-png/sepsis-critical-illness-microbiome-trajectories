$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '09_001_global_evidence_synthesis_source_anchor_fix.R')

& $Rscript $Script
exit $LASTEXITCODE
