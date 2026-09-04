$ErrorActionPreference="Stop"

$CodeRoot="E:\sepsis_project\code\03_data_processing"
$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Installer=Join-Path $CodeRoot "GET_VSEARCH_2.31.0.ps1"
$Script=Join-Path $CodeRoot "94B1_CRA002354_fasta_otu97_silva1382.R"

$Vsearch="E:\sepsis_project\tools\vsearch-2.31.0-win-x86_64\bin\vsearch.exe"

if (-not (Test-Path $Vsearch)) {
    powershell -ExecutionPolicy Bypass -File $Installer
}

& $Rscript $Script

exit $LASTEXITCODE
