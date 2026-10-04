$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '06_012_94b2a3_source_model_spec_qc_fix.R')

& $Rscript $Script

exit $LASTEXITCODE
