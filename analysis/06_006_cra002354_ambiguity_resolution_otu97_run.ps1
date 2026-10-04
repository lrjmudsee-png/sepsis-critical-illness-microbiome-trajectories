$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '06_006_cra002354_ambiguity_resolution_otu97.R')

& $Rscript $Script

exit $LASTEXITCODE
