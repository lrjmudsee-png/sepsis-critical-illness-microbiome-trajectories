$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script="E:\sepsis_project\code\03_data_processing\93K_genus_matrix_validation.R"

& $Rscript $Script

exit $LASTEXITCODE
