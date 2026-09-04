
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\V2_36I_cross_cohort_taxonomic_reproducibility.R"

$Out="E:\sepsis_project\results\V2_36I_CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY"
$Zip="E:\sepsis_project\results\V2_36I_CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY_FOR_UPLOAD.zip"

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

if($p.ExitCode -ne 0) {
  Write-Host ""
  Write-Host "R STDERR:" -ForegroundColor Red
  if(Test-Path $Stderr) {
    Get-Content $Stderr
  }
  throw "Step36I R script failed."
}

$Required=@(
  "00_INPUT_REGISTRY.csv",
  "01_COMPARISON_DEFINITION_REGISTRY.csv",
  "02_HARMONIZED_WITHIN_COHORT_TAXON_EFFECTS.csv",
  "03_PRIMARY3_PAIRWISE_TAXONOMIC_CONCORDANCE.csv",
  "04_PRIMARY3_PAIRWISE_SENSITIVITY_ALL_TAXA.csv",
  "05_PRIMARY3_SHARED_TAXA_DIRECTION_PATTERNS.csv",
  "06_PRIMARY3_THREEWAY_DIRECTION_SUMMARY.csv",
  "07_CONTEXT_EXTERNAL_AND_CONTROL_CONCORDANCE.csv",
  "08_ECOLOGICAL_VS_TAXONOMIC_SYNTHESIS.csv",
  "09_EVIDENCE_DECISION.csv",
  "10_Figure_candidate_pairwise_effect_concordance.pdf",
  "10_Figure_candidate_pairwise_effect_concordance.png",
  "11_MANUSCRIPT_SAFE_SUMMARY.txt",
  "_STEP36I_COMPLETE.txt"
)

foreach($f in $Required) {
  $fp=Join-Path $Out $f
  if(!(Test-Path $fp)) {
    throw "Missing Step36I output: $fp"
  }
}

Write-Host ""
Write-Host "STEP36I RESULT:" -ForegroundColor Cyan
Get-Content (Join-Path $Out "11_MANUSCRIPT_SAFE_SUMMARY.txt")

if(Test-Path $Zip) {
  Remove-Item $Zip -Force
}

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP36I COMPLETE" -ForegroundColor Green
Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
