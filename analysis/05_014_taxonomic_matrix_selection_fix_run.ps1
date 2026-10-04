$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_014_taxonomic_matrix_selection_fix.R')

& $Rscript $Script

exit $LASTEXITCODE
