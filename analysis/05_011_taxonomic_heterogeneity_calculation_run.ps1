$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_011_taxonomic_heterogeneity_calculation.R')

& $Rscript $Script

exit $LASTEXITCODE
