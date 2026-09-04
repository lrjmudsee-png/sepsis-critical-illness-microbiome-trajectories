$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\94A_infection_source_feasibility_audit.R"

& $Rscript $Script

exit $LASTEXITCODE
