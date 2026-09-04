
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script="E:\sepsis_project\code\03_data_processing\93N4_patient_timepoint_reconstruction.R"

& $Rscript $Script

exit $LASTEXITCODE
