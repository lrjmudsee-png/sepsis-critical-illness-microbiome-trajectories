
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script="E:\sepsis_project\code\03_data_processing\V2_96C_discussion_integration.R"

$Out="E:\sepsis_project\results\V2_96C_DISCUSSION_INTEGRATION"

$Zip="E:\sepsis_project\results\V2_96C_DISCUSSION_INTEGRATION_FOR_UPLOAD.zip"

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

  throw "Step96C R script failed."
}

$Required=@(
  "Discussion_Draft_v2_WITH_CROSS_COHORT_REPRODUCIBILITY.txt",
  "01_DISCUSSION_INTEGRATION_AUDIT.csv",
  "03_STEP96C_SUMMARY.txt",
  "_STEP96C_COMPLETE.txt"
)

foreach($f in $Required){

  $fp=Join-Path $Out $f

  if(!(Test-Path $fp)){
    throw "Missing Step96C output: $fp"
  }
}

# Copy verified literature map from code package if present.
$Map="E:\sepsis_project\code\03_data_processing\V2_96C_VERIFIED_LITERATURE_MAP.csv"

if(Test-Path $Map){
  Copy-Item $Map (Join-Path $Out "02_VERIFIED_LITERATURE_MAP.csv") -Force
}else{
  throw "Missing V2_96C_VERIFIED_LITERATURE_MAP.csv in code directory."
}

Write-Host ""
Write-Host "STEP96C SUMMARY:" -ForegroundColor Cyan
Get-Content (Join-Path $Out "03_STEP96C_SUMMARY.txt")

if(Test-Path $Zip){
  Remove-Item $Zip -Force
}

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP96C COMPLETE" -ForegroundColor Green
Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
