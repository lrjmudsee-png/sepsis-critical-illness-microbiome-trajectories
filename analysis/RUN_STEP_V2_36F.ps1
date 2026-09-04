
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\V2_36F_PRJNA851469_clinical_recovery_linkage.R"
$Out="E:\sepsis_project\results\V2_36F_PRJNA851469_CLINICAL_RECOVERY_LINKAGE"
$Zip="E:\sepsis_project\results\V2_36F_PRJNA851469_CLINICAL_RECOVERY_LINKAGE_FOR_UPLOAD.zip"

& $Rscript $Script

$Required=@(
  "01_patient_clinical_mapping_curated.csv",
  "02_PRJNA851469_displacement_with_recovered_clinical.csv",
  "02B_unmatched_ICU_patient_ids.csv",
  "03_day3_landmark_candidate.csv",
  "04_analysis_readiness_summary.csv",
  "05_linkage_QC_summary.csv",
  "06_DECISION_SUMMARY.txt",
  "_STEP36F_COMPLETE.txt"
)

foreach ($f in $Required) {
  $p=Join-Path $Out $f
  if (!(Test-Path $p)) {
    throw "Missing Step36F output: $p"
  }
}

Write-Host ""
Write-Host "STEP36F DECISION:" -ForegroundColor Cyan
Get-Content (Join-Path $Out "06_DECISION_SUMMARY.txt")

if (Test-Path $Zip) { Remove-Item $Zip -Force }
Compress-Archive -Path "$Out\*" -DestinationPath $Zip -CompressionLevel Optimal

Write-Host ""
Write-Host "Upload:" -ForegroundColor Green
Write-Host $Zip -ForegroundColor Cyan
