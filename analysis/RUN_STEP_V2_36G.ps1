
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\V2_36G_PRJNA851469_day3_clinical_outcome_models.R"
$Out="E:\sepsis_project\results\V2_36G_PRJNA851469_DAY3_CLINICAL_OUTCOME_MODELS"
$Zip="E:\sepsis_project\results\V2_36G_PRJNA851469_DAY3_CLINICAL_OUTCOME_MODELS_FOR_UPLOAD.zip"

& $Rscript $Script

$Required=@(
  "01_day3_landmark_analysis_population.csv",
  "02_descriptive_by_subsequent_event.csv",
  "03_cox_landmark_models.csv",
  "04_cox_PH_diagnostics.csv",
  "05_logistic_horizon_models.csv",
  "06_bootstrap_primary_cox.csv",
  "07_permutation_primary_cox.csv",
  "08_EVIDENCE_DECISION.csv",
  "09_Figure_candidate_day3_displacement_event.pdf",
  "09_Figure_candidate_day3_displacement_event.png",
  "10_MANUSCRIPT_SAFE_SUMMARY.txt",
  "_STEP36G_COMPLETE.txt"
)

foreach ($f in $Required) {
  $p=Join-Path $Out $f
  if (!(Test-Path $p)) {
    throw "Missing Step36G output: $p"
  }
}

Write-Host ""
Write-Host "STEP36G RESULT:" -ForegroundColor Cyan
Get-Content (Join-Path $Out "10_MANUSCRIPT_SAFE_SUMMARY.txt")

if (Test-Path $Zip) { Remove-Item $Zip -Force }

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "Upload:" -ForegroundColor Green
Write-Host $Zip -ForegroundColor Cyan
