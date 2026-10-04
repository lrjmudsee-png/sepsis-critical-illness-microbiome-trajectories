$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '09_006_final_figure_formatting_fix.R')

& $Rscript $Script
exit $LASTEXITCODE
