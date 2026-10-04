
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_025_trajectory_patient_id_reconstruction_fix.R')

& $Rscript $Script

exit $LASTEXITCODE
