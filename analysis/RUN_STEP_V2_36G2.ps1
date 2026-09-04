
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\V2_36G2_PRJNA851469_day3_clinical_outcome_models_FIXED.R"

$Out="E:\sepsis_project\results\V2_36G2_PRJNA851469_DAY3_CLINICAL_OUTCOME_MODELS_FIXED"
$Zip="E:\sepsis_project\results\V2_36G2_PRJNA851469_DAY3_CLINICAL_OUTCOME_MODELS_FIXED_FOR_UPLOAD.zip"

New-Item -ItemType Directory -Force -Path $Out | Out-Null

$Stdout=Join-Path $Out "_RUN_STDOUT.txt"
$Stderr=Join-Path $Out "_RUN_STDERR.txt"

$p=Start-Process `
  -FilePath $Rscript `
  -ArgumentList "`"$Script`"" `
  -NoNewWindow `
  -Wait `
  -PassThru `
  -RedirectStandardOutput $Stdout `
  -RedirectStandardError $Stderr

Write-Host "R exit code:" $p.ExitCode

if ($p.ExitCode -ne 0) {

  Write-Host ""
  Write-Host "R STDERR:" -ForegroundColor Red

  if (Test-Path $Stderr) {
    Get-Content $Stderr
  }

  throw "Step36G2 R script failed."
}

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
  "_STEP36G2_COMPLETE.txt"
)

foreach ($f in $Required) {

  $fp=Join-Path $Out $f

  if (!(Test-Path $fp)) {
    throw "Missing Step36G2 output: $fp"
  }
}

Write-Host ""
Write-Host "STEP36G2 RESULT:" -ForegroundColor Cyan
Get-Content (Join-Path $Out "10_MANUSCRIPT_SAFE_SUMMARY.txt")

if (Test-Path $Zip) {
  Remove-Item $Zip -Force
}

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP36G2 COMPLETE" -ForegroundColor Green
Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
