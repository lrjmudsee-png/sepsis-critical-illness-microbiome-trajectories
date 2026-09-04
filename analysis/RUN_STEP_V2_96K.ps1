
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script="E:\sepsis_project\code\03_data_processing\V2_96K_final_manuscript_text_and_audit.R"

$Out="E:\sepsis_project\results\V2_96K_FINAL_MANUSCRIPT_TEXT_AND_AUDIT"

$Zip="E:\sepsis_project\results\V2_96K_FINAL_MANUSCRIPT_TEXT_AND_AUDIT_FOR_UPLOAD.zip"

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

if($p.ExitCode -ne 0){

  Write-Host ""
  Write-Host "R STDERR:" -ForegroundColor Red

  if(Test-Path $Stderr){
    Get-Content $Stderr
  }

  throw "Step96K R script failed."
}

$Required=@(
  "Manuscript_v5_SUBMISSION_DRAFT.txt",
  "01_ABSTRACT_NUMERIC_LOCK_QC.csv",
  "02_FINAL_INTERNAL_LANGUAGE_QC.csv",
  "03_SCIENTIFIC_GUARDRAIL_QC.csv",
  "04_FINAL_REFERENCE_NUMBER_MAP.csv",
  "05_FINAL_MAIN_FIGURE_ASSET_REGISTRY.csv",
  "06_TEXT_LENGTH_METRICS.csv",
  "07_SUBMISSION_READINESS_AUDIT.csv",
  "08_STEP96K_SUMMARY.txt",
  "_STEP96K_COMPLETE.txt"
)

foreach($f in $Required){

  $fp=Join-Path $Out $f

  if(!(Test-Path $fp)){
    throw "Missing Step96K output: $fp"
  }
}

Write-Host ""
Write-Host "STEP96K SUMMARY:" -ForegroundColor Cyan
Get-Content (Join-Path $Out "08_STEP96K_SUMMARY.txt")

if(Test-Path $Zip){
  Remove-Item $Zip -Force
}

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP96K COMPLETE" -ForegroundColor Green
Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
