$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\94A3_master_infection_source_project_resolution.R"

& $Rscript $Script

exit $LASTEXITCODE
