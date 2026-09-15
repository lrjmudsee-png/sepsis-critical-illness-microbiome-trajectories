param(
    [Parameter(Mandatory = $true)]
    [string]$OutputDirectory
)

$ErrorActionPreference = 'Stop'
$python = 'C:\Users\30657\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
$rscript = 'C:\Program Files\R\R-4.4.0\bin\Rscript.exe'
$code = 'E:\sepsis_project\code\03_data_processing\99_DIRECTION_UPGRADE'
$output = [System.IO.Path]::GetFullPath($OutputDirectory)
$allowedParent = [System.IO.Path]::GetFullPath('E:\sepsis_project\results\V2_UPGRADE_20260913_DIRECTION')

if (-not $output.StartsWith($allowedParent, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Output must stay under $allowedParent"
}
if (Test-Path -LiteralPath $output) {
    throw "Refusing to overwrite existing output directory: $output"
}

New-Item -ItemType Directory -Path $output | Out-Null
try {
    & $python (Join-Path $code '99B0_initialize.py') --output $output
    if ($LASTEXITCODE -ne 0) { throw "99B0_initialize.py failed with exit code $LASTEXITCODE" }

    & $rscript (Join-Path $code '99B1_direction_analysis.R') --output $output
    if ($LASTEXITCODE -ne 0) { throw "99B1_direction_analysis.R failed with exit code $LASTEXITCODE" }

    & $python (Join-Path $code '99B2_finalize.py') --output $output
    if ($LASTEXITCODE -ne 0) { throw "99B2_finalize.py failed with exit code $LASTEXITCODE" }
}
catch {
    $failure = Join-Path $output '_RUN_FAILED.txt'
    $message = @(
        "status=FAILED",
        "output=$output",
        "error=$($_.Exception.Message)",
        "No output was deleted; this directory is retained for diagnosis."
    )
    [System.IO.File]::WriteAllLines($failure, $message, [System.Text.UTF8Encoding]::new($false))
    throw
}

Write-Host "Stage-B result package completed: $output"
