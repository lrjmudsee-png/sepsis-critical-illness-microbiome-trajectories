$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script="E:\sepsis_project\code\03_data_processing\92B_EII_construction.R"

if(!(Test-Path $Rscript)){
    Write-Host "Rscript missing:"
    Write-Host $Rscript
    exit 1
}

& $Rscript $Script
exit $LASTEXITCODE
