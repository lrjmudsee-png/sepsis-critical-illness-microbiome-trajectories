
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script="E:\sepsis_project\code\03_data_processing\V2_96F_manuscript_core_assembly_and_gap_audit.R"

$Out="E:\sepsis_project\results\V2_96F_MANUSCRIPT_CORE_ASSEMBLY_AND_GAP_AUDIT"

$Zip="E:\sepsis_project\results\V2_96F_MANUSCRIPT_CORE_ASSEMBLY_AND_GAP_AUDIT_FOR_UPLOAD.zip"

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

  throw "Step96F R script failed."
}

$Required=@(
  "01_MANUSCRIPT_CORE_DRAFT_v1.txt",
  "Main_Table_1_Cohort_Characteristics_MANUSCRIPT.csv",
  "03_FIGURE_ASSET_REGISTRY.csv",
  "05_UNRESOLVED_REFERENCE_PLACEHOLDERS.csv",
  "06_INTERNAL_LANGUAGE_AUDIT.csv",
  "07_MANUSCRIPT_GAP_AUDIT.csv",
  "08_NEXT_STAGE_DECISION.txt",
  "_STEP96F_COMPLETE.txt"
)

foreach($f in $Required){

  $fp=Join-Path $Out $f

  if(!(Test-Path $fp)){
    throw "Missing Step96F output: $fp"
  }
}

Write-Host ""
Write-Host "STEP96F DECISION:" -ForegroundColor Cyan
Get-Content (Join-Path $Out "08_NEXT_STAGE_DECISION.txt")

if(Test-Path $Zip){
  Remove-Item $Zip -Force
}

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP96F COMPLETE" -ForegroundColor Green
Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
