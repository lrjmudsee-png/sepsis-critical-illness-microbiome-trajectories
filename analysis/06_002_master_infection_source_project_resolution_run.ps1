$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '06_002_master_infection_source_project_resolution.R')

& $Rscript $Script

exit $LASTEXITCODE
