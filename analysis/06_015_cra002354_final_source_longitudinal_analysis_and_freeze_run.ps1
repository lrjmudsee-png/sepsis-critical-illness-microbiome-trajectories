$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '06_015_cra002354_final_source_longitudinal_analysis_and_freeze.R')

& $Rscript $Script
exit $LASTEXITCODE
