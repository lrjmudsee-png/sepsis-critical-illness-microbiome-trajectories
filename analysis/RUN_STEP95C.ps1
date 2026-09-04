$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\95C_manuscript_architecture_figure1_and_main_tables.R"

& $Rscript $Script
exit $LASTEXITCODE
