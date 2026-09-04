
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\V2_36I2_comparison_selection_fix.R"

& $Rscript $Script

$Out="E:\sepsis_project\results\V2_36I2_CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY_FIXED"

if(!(Test-Path "$Out\01_COMPARISON_SELECTION_AUDIT.csv")){
throw "Missing Step36I2 output"
}

Write-Host "STEP36I2 COMPLETE"
