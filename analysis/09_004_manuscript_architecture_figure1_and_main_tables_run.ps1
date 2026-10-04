$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '09_004_manuscript_architecture_figure1_and_main_tables.R')

& $Rscript $Script
exit $LASTEXITCODE
