
$ErrorActionPreference="Stop"
$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '05_040_fix_always_output.R')
& $Rscript $Script
exit $LASTEXITCODE
