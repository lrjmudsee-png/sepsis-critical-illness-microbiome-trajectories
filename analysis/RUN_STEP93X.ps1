$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\93X_PRJNA1010969_external_state_validation.R"

& $Rscript $Script

exit $LASTEXITCODE
