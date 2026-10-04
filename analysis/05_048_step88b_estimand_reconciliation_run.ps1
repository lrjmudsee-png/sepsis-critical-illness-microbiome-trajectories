$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '05_048_step88b_estimand_reconciliation.R')

& $Rscript $Script

exit $LASTEXITCODE
