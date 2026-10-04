
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_037_patient_id_mapping_audit.R')

& $Rscript $Script

exit $LASTEXITCODE
