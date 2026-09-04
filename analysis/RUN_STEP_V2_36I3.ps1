
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script="E:\sepsis_project\code\03_data_processing\V2_36I3_cross_cohort_taxonomic_reproducibility_FINAL.R"

$Out="E:\sepsis_project\results\V2_36I3_CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY_FINAL"

$Zip="E:\sepsis_project\results\V2_36I3_CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY_FINAL_FOR_UPLOAD.zip"

New-Item `
  -ItemType Directory `
  -Force `
  -Path $Out | Out-Null

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

if($p.ExitCode -ne 0){

  Write-Host ""
  Write-Host "R STDERR:" -ForegroundColor Red

  if(Test-Path $Stderr){
    Get-Content $Stderr
  }

  throw "Step36I3 R script failed."
}

$Required=@(
  "00_INPUT_REGISTRY.csv",
  "01_COMPARISON_SELECTION_AUDIT.csv",
  "02_HARMONIZED_WITHIN_COHORT_EFFECTS.csv",
  "03_PRIMARY3_PAIRWISE_CONCORDANCE.csv",
  "04_ALL_TAXA_SENSITIVITY_PAIRWISE_CONCORDANCE.csv",
  "05_PRIMARY3_SHARED_TAXA_DIRECTION_PATTERNS.csv",
  "06_THREEWAY_DIRECTION_SUMMARY.csv",
  "07_CONTEXT_EXTERNAL_CONTROL_CONCORDANCE.csv",
  "08_ECOLOGICAL_VS_TAXONOMIC_SYNTHESIS.csv",
  "09_COMMON_ANCHOR_SENSITIVITY_SUMMARY.csv",
  "10_EVIDENCE_DECISION.csv",
  "11_Figure_candidate_taxonomic_concordance.pdf",
  "11_Figure_candidate_taxonomic_concordance.png",
  "12_MANUSCRIPT_SAFE_SUMMARY.txt",
  "_STEP36I3_COMPLETE.txt"
)

foreach($f in $Required){

  $fp=Join-Path $Out $f

  if(!(Test-Path $fp)){
    throw "Missing Step36I3 output: $fp"
  }
}

Write-Host ""
Write-Host "STEP36I3 RESULT:" -ForegroundColor Cyan

Get-Content (
  Join-Path `
    $Out `
    "12_MANUSCRIPT_SAFE_SUMMARY.txt"
)

if(Test-Path $Zip){
  Remove-Item $Zip -Force
}

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP36I3 COMPLETE" -ForegroundColor Green
Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
