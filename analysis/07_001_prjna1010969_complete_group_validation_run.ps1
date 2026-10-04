$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '07_001_prjna1010969_complete_group_validation.R')

& $Rscript $Script

exit $LASTEXITCODE
