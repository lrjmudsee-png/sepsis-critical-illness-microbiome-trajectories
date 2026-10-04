
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_023_taxonomic_abundance_metadata_integration.R')

& $Rscript $Script

exit $LASTEXITCODE
