
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script="E:\sepsis_project\code\03_data_processing\V2_96I_final_reference_numbering.R"

$Out="E:\sepsis_project\results\V2_96I_FINAL_REFERENCE_NUMBERING"

$Zip="E:\sepsis_project\results\V2_96I_FINAL_REFERENCE_NUMBERING_FOR_UPLOAD.zip"

New-Item `
  -ItemType Directory `
  -Force `
  -Path $Out | Out-Null

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

  throw "Step96I R script failed."
}

$Required=@(
  "Manuscript_Core_v4_NUMBERED_REFERENCES.txt",
  "01_FINAL_REFERENCE_NUMBER_MAP.csv",
  "02_CITATION_BLOCK_TO_NUMBER_MAP.csv",
  "03_FINAL_NUMBERED_REFERENCE_LIST.txt",
  "04_FINAL_REFERENCE_NUMBERING_QC.csv",
  "05_STEP96I_SUMMARY.txt",
  "_STEP96I_COMPLETE.txt"
)

foreach($f in $Required){

  $fp=Join-Path $Out $f

  if(!(Test-Path $fp)){
    throw "Missing Step96I output: $fp"
  }
}

Write-Host ""
Write-Host "STEP96I SUMMARY:" -ForegroundColor Cyan
Get-Content (
  Join-Path $Out "05_STEP96I_SUMMARY.txt"
)

if(Test-Path $Zip){
  Remove-Item $Zip -Force
}

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP96I COMPLETE" -ForegroundColor Green
Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
