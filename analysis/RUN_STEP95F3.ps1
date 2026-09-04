$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\95F3_results_and_table1_final_freeze_sentence_audit.R"

$Out="E:\sepsis_project\results\V2_35F3_RESULTS_AND_TABLE1_FINAL_FREEZE_SENTENCE_AUDIT"
New-Item -ItemType Directory -Force -Path $Out | Out-Null

$stdout=Join-Path $Out "_RUN_STEP95F3_STDOUT.txt"
$stderr=Join-Path $Out "_RUN_STEP95F3_STDERR.txt"

$p=Start-Process `
  -FilePath $Rscript `
  -ArgumentList "`"$Script`"" `
  -NoNewWindow `
  -Wait `
  -PassThru `
  -RedirectStandardOutput $stdout `
  -RedirectStandardError $stderr

Write-Host "Step95F3 exit code:" $p.ExitCode
Write-Host "Stdout:" $stdout
Write-Host "Stderr:" $stderr

exit $p.ExitCode
