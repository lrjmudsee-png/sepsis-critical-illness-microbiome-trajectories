$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_007_taxonomic_trajectory_heterogeneity.R')

& $Rscript $Script

exit $LASTEXITCODE
