$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '07_004_prjna1010969_external_state_validation.R')

& $Rscript $Script

exit $LASTEXITCODE
