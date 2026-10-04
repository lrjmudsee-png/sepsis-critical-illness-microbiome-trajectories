
$ErrorActionPreference="Stop"
$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '05_027_fix_na_trajectory.R')
& $Rscript $Script
exit $LASTEXITCODE
