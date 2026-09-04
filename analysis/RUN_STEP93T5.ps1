
$ErrorActionPreference="Stop"
$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\93T5_EII_all_format_discovery.R"
& $Rscript $Script
exit $LASTEXITCODE
