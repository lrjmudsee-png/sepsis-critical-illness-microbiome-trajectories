
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script="E:\sepsis_project\code\03_data_processing\V2_96D2_methods_source_verification_audit.R"

$Out="E:\sepsis_project\results\V2_96D2_METHODS_SOURCE_VERIFICATION_AUDIT"

$Zip="E:\sepsis_project\results\V2_96D2_METHODS_SOURCE_VERIFICATION_AUDIT_FOR_UPLOAD.zip"

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

  throw "Step96D2 R script failed."
}

$Required=@(
  "01A_VERIFY1_PROCESSING_SCRIPT_CANDIDATES.csv",
  "01B_VERIFY1_UPSTREAM_PROCESSING_EVIDENCE.csv",
  "01C_VERIFY1_CRA002354_PROCESSING_EVIDENCE.csv",
  "02_VERIFY2_STEP88A2_MODEL_EVIDENCE.csv",
  "02B_STEP88A2_NUMBERED_SOURCE.txt",
  "03_VERIFY3_STEP94B2_MODEL_EVIDENCE.csv",
  "04A_PACKAGE_USAGE_IN_SOURCE_CODE.csv",
  "04_VERIFY4_SOFTWARE_VERSIONS.csv",
  "05_VERIFY_RESOLUTION_SUMMARY.txt",
  "_STEP96D2_COMPLETE.txt"
)

foreach($f in $Required){

  $fp=Join-Path $Out $f

  if(!(Test-Path $fp)){
    throw "Missing Step96D2 output: $fp"
  }
}

Write-Host ""
Write-Host "STEP96D2 SUMMARY:" -ForegroundColor Cyan
Get-Content (Join-Path $Out "05_VERIFY_RESOLUTION_SUMMARY.txt")

if(Test-Path $Zip){
  Remove-Item $Zip -Force
}

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP96D2 COMPLETE" -ForegroundColor Green
Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
