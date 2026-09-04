
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\V2_96I2_reference_format_polish.R"

$Out="E:\sepsis_project\results\V2_96I2_REFERENCE_FORMAT_POLISH"
$Zip="E:\sepsis_project\results\V2_96I2_REFERENCE_FORMAT_POLISH_FOR_UPLOAD.zip"

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

  throw "Step96I2 R script failed."
}

$Required=@(
  "01_FINAL_NUMBERED_REFERENCE_LIST_POLISHED.txt",
  "02_FINAL_REFERENCE_NUMBER_MAP_UNCHANGED.csv",
  "03_REFERENCE_FORMAT_QC.csv",
  "04_STEP96I2_SUMMARY.txt",
  "Manuscript_Core_v4_1_REFERENCES_POLISHED.txt",
  "_STEP96I2_COMPLETE.txt"
)

foreach($f in $Required){

  $fp=Join-Path $Out $f

  if(!(Test-Path $fp)){
    throw "Missing Step96I2 output: $fp"
  }
}

Write-Host ""
Write-Host "STEP96I2 SUMMARY:" -ForegroundColor Cyan
Get-Content (Join-Path $Out "04_STEP96I2_SUMMARY.txt")

if(Test-Path $Zip){
  Remove-Item $Zip -Force
}

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP96I2 COMPLETE" -ForegroundColor Green
Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
