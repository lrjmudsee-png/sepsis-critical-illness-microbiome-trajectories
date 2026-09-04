
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script="E:\sepsis_project\code\03_data_processing\93Q_trajectory_patient_id_reconstruction_fix.R"

& $Rscript $Script

exit $LASTEXITCODE
