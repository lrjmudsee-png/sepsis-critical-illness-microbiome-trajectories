$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '09_007_publication_ready_figures_and_table_placement.R')

& $Rscript $Script
exit $LASTEXITCODE
