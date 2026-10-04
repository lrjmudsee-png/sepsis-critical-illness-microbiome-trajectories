
$ErrorActionPreference="Stop"
$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '05_032_eii_real_table_identification.R')
& $Rscript $Script
exit $LASTEXITCODE
