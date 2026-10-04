
$ErrorActionPreference="Stop"
$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '05_034_eii_all_format_discovery.R')
& $Rscript $Script
exit $LASTEXITCODE
