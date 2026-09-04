
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\V2_96D_methods_manuscriptization.R"

$Out="E:\sepsis_project\results\V2_96D_METHODS_MANUSCRIPTIZATION"
$Zip="E:\sepsis_project\results\V2_96D_METHODS_MANUSCRIPTIZATION_FOR_UPLOAD.zip"

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

  throw "Step96D R script failed."
}

$Required=@(
  "01_METHODS_PROVENANCE_DIRECTORY_REGISTRY.csv",
  "02_RELEVANT_ANALYSIS_CODE_REGISTRY.csv",
  "03_KEY_FROZEN_FILE_REGISTRY.csv",
  "Methods_Draft_v1_FROM_FROZEN_V2_ANALYSES.txt",
  "04_METHODS_ITEMS_REQUIRING_VERIFICATION.csv",
  "05_STEP96D_SUMMARY.txt",
  "_STEP96D_COMPLETE.txt"
)

foreach($f in $Required){

  $fp=Join-Path $Out $f

  if(!(Test-Path $fp)){
    throw "Missing Step96D output: $fp"
  }
}

Write-Host ""
Write-Host "STEP96D SUMMARY:" -ForegroundColor Cyan
Get-Content (Join-Path $Out "05_STEP96D_SUMMARY.txt")

if(Test-Path $Zip){
  Remove-Item $Zip -Force
}

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP96D COMPLETE" -ForegroundColor Green
Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
