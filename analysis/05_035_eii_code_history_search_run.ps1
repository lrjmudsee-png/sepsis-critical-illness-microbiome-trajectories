
$ErrorActionPreference="Stop"
$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '05_035_eii_code_history_search.R')
& $Rscript $Script
exit $LASTEXITCODE
