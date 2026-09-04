$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\94B1E_CRA002354_taxonomy_id_repair.R"

& $Rscript $Script
exit $LASTEXITCODE
