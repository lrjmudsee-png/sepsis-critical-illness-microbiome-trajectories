$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '05_049_cross_cohort_trajectory_polarity.R')

& $Rscript $Script

exit $LASTEXITCODE
