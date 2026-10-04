$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '06_004_infection_source_feasibility_audit.R')

& $Rscript $Script

exit $LASTEXITCODE
