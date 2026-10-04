
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_045_fix_patient_id_type.R')

& $Rscript $Script

exit $LASTEXITCODE
