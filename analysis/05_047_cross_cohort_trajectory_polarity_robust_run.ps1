$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '05_047_cross_cohort_trajectory_polarity_robust.R')

& $Rscript $Script

exit $LASTEXITCODE
