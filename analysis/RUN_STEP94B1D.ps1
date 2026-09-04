$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\94B1D_CRA002354_nonchimeric_otu_table_and_rarefaction_fix.R"

& $Rscript $Script

exit $LASTEXITCODE
