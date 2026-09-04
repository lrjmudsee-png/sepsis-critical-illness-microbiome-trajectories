$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\95D_publication_ready_figures_and_table_placement.R"

& $Rscript $Script
exit $LASTEXITCODE
