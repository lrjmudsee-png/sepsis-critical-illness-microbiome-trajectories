
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\V2_96D3_methods_manuscriptization_verified.R"

$Out="E:\sepsis_project\results\V2_96D3_METHODS_MANUSCRIPTIZATION_VERIFIED"
$Zip="E:\sepsis_project\results\V2_96D3_METHODS_MANUSCRIPTIZATION_VERIFIED_FOR_UPLOAD.zip"

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

  throw "Step96D3 R script failed."
}

$Required=@(
  "Methods_Draft_v2_SOURCE_VERIFIED.txt",
  "01_VERIFY_RESOLUTION_REGISTRY.csv",
  "02_CORE_SOFTWARE_VERSIONS_USED_IN_METHODS.csv",
  "03_METHODS_V2_QC.csv",
  "04_STEP96D3_SUMMARY.txt",
  "_STEP96D3_COMPLETE.txt"
)

foreach($f in $Required){

  $fp=Join-Path $Out $f

  if(!(Test-Path $fp)){
    throw "Missing Step96D3 output: $fp"
  }
}

Write-Host ""
Write-Host "STEP96D3 SUMMARY:" -ForegroundColor Cyan
Get-Content (Join-Path $Out "04_STEP96D3_SUMMARY.txt")

if(Test-Path $Zip){
  Remove-Item $Zip -Force
}

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP96D3 COMPLETE" -ForegroundColor Green
Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
