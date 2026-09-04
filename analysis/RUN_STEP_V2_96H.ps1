
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script="E:\sepsis_project\code\03_data_processing\V2_96H_introduction_integration.R"

$NewRefs="E:\sepsis_project\code\03_data_processing\V2_96H_NEW_VERIFIED_REFERENCES.csv"

$Out="E:\sepsis_project\results\V2_96H_INTRODUCTION_INTEGRATION"

$Zip="E:\sepsis_project\results\V2_96H_INTRODUCTION_INTEGRATION_FOR_UPLOAD.zip"

if(!(Test-Path $NewRefs)){
  throw "Missing verified Introduction reference file: $NewRefs"
}

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

  throw "Step96H R script failed."
}

$Required=@(
  "Manuscript_Core_v3_WITH_INTRODUCTION.txt",
  "01_VERIFIED_REFERENCE_POOL_EXPANDED.csv",
  "02_CITATION_KEYS_USED_AFTER_INTRODUCTION.csv",
  "03_MISSING_CITATION_KEYS.csv",
  "05_INTRODUCTION_AND_INTERNAL_LANGUAGE_QC.csv",
  "06_INTRODUCTION_QC.csv",
  "07_STEP96H_SUMMARY.txt",
  "_STEP96H_COMPLETE.txt"
)

foreach($f in $Required){

  $fp=Join-Path $Out $f

  if(!(Test-Path $fp)){
    throw "Missing Step96H output: $fp"
  }
}

Write-Host ""
Write-Host "STEP96H SUMMARY:" -ForegroundColor Cyan
Get-Content (Join-Path $Out "07_STEP96H_SUMMARY.txt")

if(Test-Path $Zip){
  Remove-Item $Zip -Force
}

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP96H COMPLETE" -ForegroundColor Green
Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
