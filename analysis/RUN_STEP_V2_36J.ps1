
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script="E:\sepsis_project\code\03_data_processing\V2_36J_microbiome_enhancement_evidence_freeze.R"

$Out="E:\sepsis_project\results\V2_36J_MICROBIOME_ENHANCEMENT_EVIDENCE_FREEZE"

$Zip="E:\sepsis_project\results\V2_36J_MICROBIOME_ENHANCEMENT_EVIDENCE_FREEZE_FOR_UPLOAD.zip"

& $Rscript $Script

$Required=@(
  "01_ENHANCEMENT_BRANCH_FREEZE.csv",
  "02_FROZEN_ENHANCEMENT_METRICS.csv",
  "03_FROZEN_MAIN_TEXT_WORDING.txt",
  "04_DISCUSSION_INTEGRATION_POINTS.txt",
  "05_FIGURE_ARCHITECTURE_RECOMMENDATION.txt",
  "SOURCE_primary3_pairwise_concordance.csv",
  "SOURCE_threeway_direction_summary.csv",
  "SOURCE_taxonomic_concordance_candidate.pdf",
  "SOURCE_taxonomic_concordance_candidate.png",
  "_STEP36J_COMPLETE.txt"
)

foreach($f in $Required){

  $p=Join-Path $Out $f

  if(!(Test-Path $p)){
    throw "Missing Step36J output: $p"
  }
}

if(Test-Path $Zip){
  Remove-Item $Zip -Force
}

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP36J COMPLETE" -ForegroundColor Green
Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
