[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$RunRoot,
    [string]$RscriptPath = "C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
)

$ErrorActionPreference = "Stop"

$resolvedRunRoot = [System.IO.Path]::GetFullPath($RunRoot)
$scriptDirectory = Join-Path $resolvedRunRoot "_scripts"
$logDirectory = Join-Path $resolvedRunRoot "logs"
$summaryPath = Join-Path $resolvedRunRoot "reproduction_step_summary.csv"

if (-not (Test-Path -LiteralPath $scriptDirectory -PathType Container)) {
    throw "Script directory not found: $scriptDirectory"
}
if (-not (Test-Path -LiteralPath $RscriptPath -PathType Leaf)) {
    throw "Rscript not found: $RscriptPath"
}

New-Item -ItemType Directory -Force -Path $logDirectory | Out-Null

$rows = [System.Collections.Generic.List[object]]::new()
$scripts = Get-ChildItem -LiteralPath $scriptDirectory -File -Filter "*.R" | Sort-Object Name

$savedLocaleEnvironment = @{}
foreach ($name in @("LANG", "LC_ALL", "LC_CTYPE")) {
    $savedLocaleEnvironment[$name] = [Environment]::GetEnvironmentVariable($name, "Process")
    Remove-Item -LiteralPath ("Env:" + $name) -ErrorAction SilentlyContinue
}

try {
foreach ($script in $scripts) {
    $started = Get-Date
    $stdout = Join-Path $logDirectory ($script.BaseName + ".stdout.log")
    $stderr = Join-Path $logDirectory ($script.BaseName + ".stderr.log")

    Write-Output ("START {0} {1}" -f $script.Name, $started.ToString("s"))
    $process = Start-Process -FilePath $RscriptPath -ArgumentList @("--encoding=UTF-8", $script.FullName) -NoNewWindow -Wait -PassThru -RedirectStandardOutput $stdout -RedirectStandardError $stderr
    $finished = Get-Date

    $rows.Add([pscustomobject]@{
        Script = $script.Name
        Started = $started.ToString("s")
        Finished = $finished.ToString("s")
        ElapsedSeconds = [math]::Round(($finished - $started).TotalSeconds, 3)
        ExitCode = $process.ExitCode
        StdoutLog = [System.IO.Path]::GetRelativePath($resolvedRunRoot, $stdout).Replace("\", "/")
        StderrLog = [System.IO.Path]::GetRelativePath($resolvedRunRoot, $stderr).Replace("\", "/")
    })
    $rows | Export-Csv -LiteralPath $summaryPath -NoTypeInformation -Encoding utf8

    Write-Output ("END {0} EXIT={1} ELAPSED={2}s" -f $script.Name, $process.ExitCode, [math]::Round(($finished - $started).TotalSeconds, 1))
    if ($process.ExitCode -ne 0) {
        throw "Reproduction stopped because $($script.Name) exited with code $($process.ExitCode). See $stderr"
    }
}
}
finally {
    foreach ($name in $savedLocaleEnvironment.Keys) {
        if ($null -eq $savedLocaleEnvironment[$name]) {
            Remove-Item -LiteralPath ("Env:" + $name) -ErrorAction SilentlyContinue
        }
        else {
            [Environment]::SetEnvironmentVariable($name, $savedLocaleEnvironment[$name], "Process")
        }
    }
}

Write-Output "REPRODUCTION_CHAIN_COMPLETE"
Write-Output "Summary: $summaryPath"
