$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\95G_supplementary_tables_and_discussion_v1.R"

$Out="E:\sepsis_project\results\V2_35G_SUPPLEMENTARY_TABLES_AND_DISCUSSION_V1"
New-Item -ItemType Directory -Force -Path $Out | Out-Null

$stdout=Join-Path $Out "_RUN_STEP95G_STDOUT.txt"
$stderr=Join-Path $Out "_RUN_STEP95G_STDERR.txt"

$p=Start-Process `
  -FilePath $Rscript `
  -ArgumentList "`"$Script`"" `
  -NoNewWindow `
  -Wait `
  -PassThru `
  -RedirectStandardOutput $stdout `
  -RedirectStandardError $stderr

Write-Host "Step95G exit code:" $p.ExitCode
Write-Host "Stdout:" $stdout
Write-Host "Stderr:" $stderr

exit $p.ExitCode
