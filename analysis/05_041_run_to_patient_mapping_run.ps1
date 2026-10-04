
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_041_run_to_patient_mapping.R')

& $Rscript $Script

exit $LASTEXITCODE
