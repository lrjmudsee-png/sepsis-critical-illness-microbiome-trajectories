
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\V2_36D_value_level_clinical_feasibility.R"
$Out="E:\sepsis_project\results\V2_36D_VALUE_LEVEL_CLINICAL_FEASIBILITY"
$Zip="E:\sepsis_project\results\V2_36D_VALUE_LEVEL_CLINICAL_FEASIBILITY_FOR_UPLOAD.zip"

& $Rscript $Script

$Required=@(
  "01_value_level_clinical_variable_inventory.csv",
  "02_analysis_specific_feasibility.csv",
  "03_antimicrobial_detail_file_inventory.csv",
  "04_DECISION_SUMMARY.txt",
  "_STEP36D_COMPLETE.txt"
)

foreach ($f in $Required) {
  $p=Join-Path $Out $f
  if (!(Test-Path $p)) {
    throw "Missing Step36D output: $p"
  }
}

if (Test-Path $Zip) { Remove-Item $Zip -Force }

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP36D COMPLETE" -ForegroundColor Green
Write-Host "Upload this file:"
Write-Host $Zip -ForegroundColor Cyan
