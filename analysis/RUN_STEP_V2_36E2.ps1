
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\V2_36E2_raw_clinical_metadata_recovery_audit_fixed.R"
$Out="E:\sepsis_project\results\V2_36E2_RAW_CLINICAL_METADATA_RECOVERY_AUDIT_FIXED"
$Zip="E:\sepsis_project\results\V2_36E2_RAW_CLINICAL_METADATA_RECOVERY_AUDIT_FIXED_FOR_UPLOAD.zip"

& $Rscript $Script

$Required=@(
  "01_candidate_tabular_files.csv",
  "02_recovered_clinical_variables_value_level.csv",
  "03_cohort_level_recovery_summary.csv",
  "04_manual_mapping_shortlist.csv",
  "05_DECISION_SUMMARY.txt",
  "_STEP36E2_COMPLETE.txt"
)

foreach ($f in $Required) {
  $p=Join-Path $Out $f
  if (!(Test-Path $p)) {
    throw "Missing Step36E2 output: $p"
  }
}

Write-Host ""
Write-Host "Decision summary:" -ForegroundColor Cyan
Get-Content (Join-Path $Out "05_DECISION_SUMMARY.txt")

if (Test-Path $Zip) {
  Remove-Item $Zip -Force
}

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP36E2 COMPLETE" -ForegroundColor Green
Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
