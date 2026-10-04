$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '06_003_cra002354_rescue_and_readiness_audit.R')

& $Rscript $Script

exit $LASTEXITCODE
