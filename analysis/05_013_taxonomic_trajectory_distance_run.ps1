$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_013_taxonomic_trajectory_distance.R')

& $Rscript $Script

exit $LASTEXITCODE
