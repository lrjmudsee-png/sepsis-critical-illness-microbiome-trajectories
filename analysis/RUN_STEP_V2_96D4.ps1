
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script="E:\sepsis_project\code\03_data_processing\V2_96D4_methods_final_manuscript_polish.R"

$Out="E:\sepsis_project\results\V2_96D4_METHODS_FINAL_MANUSCRIPT_POLISH"

$Zip="E:\sepsis_project\results\V2_96D4_METHODS_FINAL_MANUSCRIPT_POLISH_FOR_UPLOAD.zip"

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

  throw "Step96D4 R script failed."
}

$Required=@(
  "Methods_v3_MANUSCRIPT_READY.txt",
  "01_INTERNAL_WORKFLOW_LANGUAGE_QC.csv",
  "02_METHODS_CONTENT_PRESERVATION_QC.csv",
  "03_METHODS_V3_FINAL_AUDIT.csv",
  "04_STEP96D4_SUMMARY.txt",
  "_STEP96D4_COMPLETE.txt"
)

foreach($f in $Required){

  $fp=Join-Path $Out $f

  if(!(Test-Path $fp)){
    throw "Missing Step96D4 output: $fp"
  }
}

Write-Host ""
Write-Host "STEP96D4 SUMMARY:" -ForegroundColor Cyan
Get-Content (Join-Path $Out "04_STEP96D4_SUMMARY.txt")

if(Test-Path $Zip){
  Remove-Item $Zip -Force
}

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP96D4 COMPLETE" -ForegroundColor Green
Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
