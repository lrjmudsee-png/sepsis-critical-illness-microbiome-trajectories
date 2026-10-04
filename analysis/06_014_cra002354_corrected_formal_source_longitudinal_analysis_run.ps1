$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '06_014_cra002354_corrected_formal_source_longitudinal_analysis.R')

& $Rscript $Script
exit $LASTEXITCODE
