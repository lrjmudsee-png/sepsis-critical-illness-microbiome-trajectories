$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_008_taxonomic_divergence_vs_eii.R')

& $Rscript $Script

exit $LASTEXITCODE
