$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\95E2_table1_analysis_population_fix_and_results_v2.R"

$Out="E:\sepsis_project\results\V2_35E2_TABLE1_ANALYSIS_POPULATION_FIX_AND_RESULTS_V2"
New-Item -ItemType Directory -Force -Path $Out | Out-Null

$stdout=Join-Path $Out "_RUN_STEP95E2_STDOUT.txt"
$stderr=Join-Path $Out "_RUN_STEP95E2_STDERR.txt"

$p=Start-Process `
  -FilePath $Rscript `
  -ArgumentList "`"$Script`"" `
  -NoNewWindow `
  -Wait `
  -PassThru `
  -RedirectStandardOutput $stdout `
  -RedirectStandardError $stderr

Write-Host "Step95E2 exit code:" $p.ExitCode
Write-Host "Stdout:" $stdout
Write-Host "Stderr:" $stderr

exit $p.ExitCode
