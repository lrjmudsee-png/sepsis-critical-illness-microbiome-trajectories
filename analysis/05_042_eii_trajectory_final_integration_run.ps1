$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_042_eii_trajectory_final_integration.R')

& $Rscript $Script

exit $LASTEXITCODE
