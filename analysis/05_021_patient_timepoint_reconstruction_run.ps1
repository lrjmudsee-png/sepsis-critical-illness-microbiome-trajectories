
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_021_patient_timepoint_reconstruction.R')

& $Rscript $Script

exit $LASTEXITCODE
