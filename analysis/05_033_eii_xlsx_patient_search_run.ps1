
$ErrorActionPreference="Stop"
$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '05_033_eii_xlsx_patient_search.R')
& $Rscript $Script
exit $LASTEXITCODE
