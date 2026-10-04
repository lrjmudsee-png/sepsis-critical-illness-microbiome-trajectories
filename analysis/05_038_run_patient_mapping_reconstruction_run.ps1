
$ErrorActionPreference="Stop"
$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '05_038_run_patient_mapping_reconstruction.R')
& $Rscript $Script
exit $LASTEXITCODE
