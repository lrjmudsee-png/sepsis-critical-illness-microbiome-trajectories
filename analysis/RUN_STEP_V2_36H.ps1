
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\V2_36H_cross_cohort_taxonomic_reproducibility_audit.R"
$Out="E:\sepsis_project\results\V2_36H_CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY_AUDIT"
$Zip="E:\sepsis_project\results\V2_36H_CROSS_COHORT_TAXONOMIC_REPRODUCIBILITY_AUDIT_FOR_UPLOAD.zip"

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
  if(Test-Path $Stderr) {
    Get-Content $Stderr
  }
  throw "Step36H R script failed."
}

$Required=@(
  "01_taxonomic_table_candidate_inventory.csv",
  "02_scored_taxonomic_candidates.csv",
  "03_best_candidates_by_project.csv",
  "04_longitudinal_metadata_candidate_inventory.csv",
  "05_cohort_taxonomic_reproducibility_feasibility.csv",
  "06_DECISION_SUMMARY.txt",
  "_STEP36H_COMPLETE.txt"
)

foreach($f in $Required) {
  $p2=Join-Path $Out $f
  if(!(Test-Path $p2)) {
    throw "Missing Step36H output: $p2"
  }
}

Write-Host ""
Write-Host "STEP36H DECISION:" -ForegroundColor Cyan
Get-Content (Join-Path $Out "06_DECISION_SUMMARY.txt")

if(Test-Path $Zip) {
  Remove-Item $Zip -Force
}

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "Upload:" -ForegroundColor Green
Write-Host $Zip -ForegroundColor Cyan
