
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Code="E:\sepsis_project\code\03_data_processing\V2_36C_clinical_variable_characterization.R"
$Out="E:\sepsis_project\results\V2_36C_CLINICAL_VARIABLE_CHARACTERIZATION"
$Zip="E:\sepsis_project\results\V2_36C_CLINICAL_VARIABLE_CHARACTERIZATION_FOR_UPLOAD.zip"

$Required=@(
    "clinical_variable_inventory_detailed.csv",
    "enhancement_candidate_table.csv",
    "_STEP36C_COMPLETE.txt"
)

Write-Host "Checking Step36C output..." -ForegroundColor Cyan

$needRun=$false

if (!(Test-Path $Out)) {
    Write-Host "Output folder does not exist." -ForegroundColor Yellow
    $needRun=$true
}
else {
    foreach ($f in $Required) {
        $p=Join-Path $Out $f
        if (!(Test-Path $p)) {
            Write-Host "Missing: $p" -ForegroundColor Yellow
            $needRun=$true
        }
    }
}

if ($needRun) {

    if (!(Test-Path $Code)) {
        throw "Step36C R script not found: $Code"
    }

    Write-Host "Step36C outputs incomplete. Re-running Step36C..." -ForegroundColor Cyan

    New-Item -ItemType Directory -Force -Path $Out | Out-Null

    $stdout=Join-Path $Out "_RUN_STEP36C_STDOUT.txt"
    $stderr=Join-Path $Out "_RUN_STEP36C_STDERR.txt"

    $p=Start-Process `
        -FilePath $Rscript `
        -ArgumentList "`"$Code`"" `
        -NoNewWindow `
        -Wait `
        -PassThru `
        -RedirectStandardOutput $stdout `
        -RedirectStandardError $stderr

    Write-Host "R exit code:" $p.ExitCode

    if ($p.ExitCode -ne 0) {
        Write-Host "STDERR:" -ForegroundColor Red
        if (Test-Path $stderr) {
            Get-Content $stderr
        }
        throw "Step36C failed."
    }
}

Write-Host ""
Write-Host "Final file verification:" -ForegroundColor Cyan

$allOK=$true
foreach ($f in $Required) {
    $p=Join-Path $Out $f
    if (Test-Path $p) {
        $item=Get-Item $p
        Write-Host ("PASS  {0}  {1} bytes" -f $f,$item.Length) -ForegroundColor Green
    }
    else {
        Write-Host ("FAIL  {0}" -f $f) -ForegroundColor Red
        $allOK=$false
    }
}

if (!$allOK) {
    throw "Step36C output verification failed. Do not upload an empty ZIP."
}

Write-Host ""
Write-Host "Previewing enhancement_candidate_table.csv:" -ForegroundColor Cyan
Import-Csv (Join-Path $Out "enhancement_candidate_table.csv") | Format-Table -AutoSize

if (Test-Path $Zip) {
    Remove-Item $Zip -Force
}

Compress-Archive `
    -Path "$Out\*" `
    -DestinationPath $Zip `
    -CompressionLevel Optimal

$zipItem=Get-Item $Zip

Write-Host ""
Write-Host "UPLOAD THIS ZIP:" -ForegroundColor Green
Write-Host $Zip
Write-Host ("ZIP size: {0:N0} bytes" -f $zipItem.Length)
