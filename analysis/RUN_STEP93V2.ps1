$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\93V2_cross_cohort_trajectory_polarity_robust.R"

& $Rscript $Script

exit $LASTEXITCODE
