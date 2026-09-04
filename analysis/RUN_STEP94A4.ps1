$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\94A4_CRA002354_rescue_and_readiness_audit.R"

& $Rscript $Script

exit $LASTEXITCODE
