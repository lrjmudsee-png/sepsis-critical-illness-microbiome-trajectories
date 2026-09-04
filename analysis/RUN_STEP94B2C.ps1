$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\94B2C_CRA002354_final_source_longitudinal_analysis_and_freeze.R"

& $Rscript $Script
exit $LASTEXITCODE
