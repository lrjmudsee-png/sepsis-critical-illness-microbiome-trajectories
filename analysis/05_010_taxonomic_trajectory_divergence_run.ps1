$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_010_taxonomic_trajectory_divergence.R')

& $Rscript $Script

exit $LASTEXITCODE
