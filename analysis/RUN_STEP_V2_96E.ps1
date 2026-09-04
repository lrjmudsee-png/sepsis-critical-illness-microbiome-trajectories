
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script="E:\sepsis_project\code\03_data_processing\V2_96E_figure3D_cross_cohort_reproducibility_final.R"

$Out="E:\sepsis_project\results\V2_96E_FIGURE3D_CROSS_COHORT_REPRODUCIBILITY_FINAL"

$Zip="E:\sepsis_project\results\V2_96E_FIGURE3D_CROSS_COHORT_REPRODUCIBILITY_FINAL_FOR_UPLOAD.zip"

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

  throw "Step96E R script failed."
}

$Required=@(
  "01_FIGURE3D_MAIN_DATA.csv",
  "02_SUPPLEMENTARY_SENSITIVITY_DATA.csv",
  "Figure3D_Cross_Cohort_Taxonomic_Concordance_MAIN.pdf",
  "Figure3D_Cross_Cohort_Taxonomic_Concordance_MAIN.png",
  "Supplementary_Figure_Cross_Cohort_Taxonomic_Concordance_Sensitivity.pdf",
  "Supplementary_Figure_Cross_Cohort_Taxonomic_Concordance_Sensitivity.png",
  "03_FIGURE3D_LEGEND.txt",
  "04_SUPPLEMENTARY_FIGURE_LEGEND.txt",
  "05_FIGURE_INTEGRATION_MANIFEST.csv",
  "06_FIGURE3D_QC.csv",
  "_STEP96E_COMPLETE.txt"
)

foreach($f in $Required){

  $fp=Join-Path $Out $f

  if(!(Test-Path $fp)){
    throw "Missing Step96E output: $fp"
  }
}

Write-Host ""
Write-Host "STEP96E COMPLETE" -ForegroundColor Green

if(Test-Path $Zip){
  Remove-Item $Zip -Force
}

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
