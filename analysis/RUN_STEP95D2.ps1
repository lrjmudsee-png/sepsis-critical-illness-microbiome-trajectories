$ErrorActionPreference = "Stop"

$Rscript = "C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script = "E:\sepsis_project\code\03_data_processing\95D2_publication_figures_sourcepack_only.R"

$LogDir = "E:\sepsis_project\results\V2_35D2_PUBLICATION_FIGURES_SOURCEPACK_ONLY"
New-Item -ItemType Directory -Force -Path $LogDir | Out-Null

$Stdout = Join-Path $LogDir "_RUN_STEP95D2_STDOUT.txt"
$Stderr = Join-Path $LogDir "_RUN_STEP95D2_STDERR.txt"

$p = Start-Process `
    -FilePath $Rscript `
    -ArgumentList "`"$Script`"" `
    -NoNewWindow `
    -Wait `
    -PassThru `
    -RedirectStandardOutput $Stdout `
    -RedirectStandardError $Stderr

Write-Host "Step95D2 exit code:" $p.ExitCode
Write-Host "Stdout:" $Stdout
Write-Host "Stderr:" $Stderr

exit $p.ExitCode
