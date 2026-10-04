$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_009_taxonomic_heterogeneity_index.R')

& $Rscript $Script

exit $LASTEXITCODE
