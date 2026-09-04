$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\94B2A3_source_model_spec_qc_fix.R"

& $Rscript $Script

exit $LASTEXITCODE
