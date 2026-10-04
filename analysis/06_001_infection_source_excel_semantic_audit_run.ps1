$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '06_001_infection_source_excel_semantic_audit.R')

& $Rscript $Script

exit $LASTEXITCODE
