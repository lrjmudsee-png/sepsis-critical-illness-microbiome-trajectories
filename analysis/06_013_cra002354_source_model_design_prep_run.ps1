$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '06_013_cra002354_source_model_design_prep.R')

& $Rscript $Script

exit $LASTEXITCODE
