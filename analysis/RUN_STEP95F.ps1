$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\95F_results_polish_and_table1_manuscript_freeze.R"

$Out="E:\sepsis_project\results\V2_35F_RESULTS_POLISH_AND_TABLE1_MANUSCRIPT_FREEZE"
New-Item -ItemType Directory -Force -Path $Out | Out-Null

$stdout=Join-Path $Out "_RUN_STEP95F_STDOUT.txt"
$stderr=Join-Path $Out "_RUN_STEP95F_STDERR.txt"

$p=Start-Process `
  -FilePath $Rscript `
  -ArgumentList "`"$Script`"" `
  -NoNewWindow `
  -Wait `
  -PassThru `
  -RedirectStandardOutput $stdout `
  -RedirectStandardError $stderr

Write-Host "Step95F exit code:" $p.ExitCode
exit $p.ExitCode
