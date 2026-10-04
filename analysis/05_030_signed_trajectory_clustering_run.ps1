
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_030_signed_trajectory_clustering.R')

& $Rscript $Script

exit $LASTEXITCODE
