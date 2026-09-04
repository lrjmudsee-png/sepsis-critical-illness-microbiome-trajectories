
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script="E:\sepsis_project\code\03_data_processing\V2_96B_results_integration.R"

$Out="E:\sepsis_project\results\V2_96B_RESULTS_INTEGRATION"

$Zip="E:\sepsis_project\results\V2_96B_RESULTS_INTEGRATION_FOR_UPLOAD.zip"

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

  throw "Step96B R script failed."
}

$Required=@(
  "Results_Draft_v6_WITH_CROSS_COHORT_REPRODUCIBILITY.txt",
  "01_INSERTION_AUDIT.csv",
  "02_SOURCE_METRICS_USED.csv",
  "03_STEP96B_SUMMARY.txt",
  "_STEP96B_COMPLETE.txt"
)

foreach($f in $Required){

  $fp=Join-Path $Out $f

  if(!(Test-Path $fp)){
    throw "Missing Step96B output: $fp"
  }
}

Write-Host ""
Write-Host "STEP96B SUMMARY:" -ForegroundColor Cyan
Get-Content (Join-Path $Out "03_STEP96B_SUMMARY.txt")

if(Test-Path $Zip){
  Remove-Item $Zip -Force
}

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP96B COMPLETE" -ForegroundColor Green
Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
