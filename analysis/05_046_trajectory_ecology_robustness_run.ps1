$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '05_046_trajectory_ecology_robustness.R')

& $Rscript $Script

exit $LASTEXITCODE
