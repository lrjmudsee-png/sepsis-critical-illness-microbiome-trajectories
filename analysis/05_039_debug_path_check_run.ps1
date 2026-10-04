
$ErrorActionPreference="Stop"
$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '05_039_debug_path_check.R')
& $Rscript $Script
exit $LASTEXITCODE
