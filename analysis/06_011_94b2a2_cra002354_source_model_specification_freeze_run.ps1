$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '06_011_94b2a2_cra002354_source_model_specification_freeze.R')

& $Rscript $Script

exit $LASTEXITCODE
