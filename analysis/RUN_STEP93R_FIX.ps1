
$ErrorActionPreference="Stop"
$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\93R_fix_na_trajectory.R"
& $Rscript $Script
exit $LASTEXITCODE
