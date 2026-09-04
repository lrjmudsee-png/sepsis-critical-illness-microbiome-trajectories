
$ErrorActionPreference="Stop"
$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\93T3_EII_real_table_identification.R"
& $Rscript $Script
exit $LASTEXITCODE
