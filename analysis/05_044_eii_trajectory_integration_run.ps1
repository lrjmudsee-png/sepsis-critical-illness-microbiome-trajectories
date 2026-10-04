
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_044_eii_trajectory_integration.R')

& $Rscript $Script

exit $LASTEXITCODE
