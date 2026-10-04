$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '06_017_cra002354_exploratory_genus_source_analysis.R')

& $Rscript $Script
exit $LASTEXITCODE
