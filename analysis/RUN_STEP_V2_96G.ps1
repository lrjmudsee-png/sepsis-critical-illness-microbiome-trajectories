
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script="E:\sepsis_project\code\03_data_processing\V2_96G_discussion_reference_resolution.R"

$Map="E:\sepsis_project\code\03_data_processing\V2_96G_VERIFIED_REFERENCE_POOL.csv"

$Out="E:\sepsis_project\results\V2_96G_DISCUSSION_REFERENCE_RESOLUTION"

$Zip="E:\sepsis_project\results\V2_96G_DISCUSSION_REFERENCE_RESOLUTION_FOR_UPLOAD.zip"

if(!(Test-Path $Map)){
  throw "Missing verified reference pool: $Map"
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

  throw "Step96G R script failed."
}

Copy-Item `
  $Map `
  (Join-Path $Out "V2_96G_VERIFIED_REFERENCE_POOL.csv") `
  -Force

$Required=@(
  "Manuscript_Core_v2_DISCUSSION_REFS_RESOLVED.txt",
  "01_PLACEHOLDER_TO_CITATION_KEY_MAP.csv",
  "02_UNRESOLVED_REFERENCE_PLACEHOLDERS.csv",
  "03_CITATION_KEYS_USED_IN_MANUSCRIPT.csv",
  "04_STEP96G_QC.csv",
  "05_STEP96G_SUMMARY.txt",
  "V2_96G_VERIFIED_REFERENCE_POOL.csv",
  "_STEP96G_COMPLETE.txt"
)

foreach($f in $Required){

  $fp=Join-Path $Out $f

  if(!(Test-Path $fp)){
    throw "Missing Step96G output: $fp"
  }
}

Write-Host ""
Write-Host "STEP96G SUMMARY:" -ForegroundColor Cyan
Get-Content (Join-Path $Out "05_STEP96G_SUMMARY.txt")

if(Test-Path $Zip){
  Remove-Item $Zip -Force
}

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP96G COMPLETE" -ForegroundColor Green
Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
