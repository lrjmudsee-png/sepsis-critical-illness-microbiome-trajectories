$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '05_043_forced_run_bridge.R')

& $Rscript $Script

exit $LASTEXITCODE
