$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\95E_table1_and_results_draft.R"

$Out="E:\sepsis_project\results\V2_35E_TABLE1_AND_RESULTS_DRAFT"
New-Item -ItemType Directory -Force -Path $Out | Out-Null

$stdout=Join-Path $Out "_RUN_STEP95E_STDOUT.txt"
$stderr=Join-Path $Out "_RUN_STEP95E_STDERR.txt"

$p=Start-Process `
  -FilePath $Rscript `
  -ArgumentList "`"$Script`"" `
  -NoNewWindow `
  -Wait `
  -PassThru `
  -RedirectStandardOutput $stdout `
  -RedirectStandardError $stderr

Write-Host "Step95E exit code:" $p.ExitCode
Write-Host "Stdout:" $stdout
Write-Host "Stderr:" $stderr

exit $p.ExitCode
