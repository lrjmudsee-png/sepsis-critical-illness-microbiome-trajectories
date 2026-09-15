param(
    [Parameter(Mandatory = $false)]
    [string]$OutputDir = 'E:\sepsis_project\results\V2_UPGRADE_20260913_DIRECTION\A_FEASIBILITY'
)

$ErrorActionPreference = 'Stop'
$scriptRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$pythonCandidates = @(
    'C:\Users\30657\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe',
    'python'
)
$rCandidates = @(
    'C:\Program Files\R\R-4.5.1\bin\Rscript.exe',
    'C:\Program Files\R\R-4.5.0\bin\Rscript.exe',
    'C:\Program Files\R\R-4.4.0\bin\Rscript.exe',
    'Rscript'
)

function Resolve-Executable {
    param([string[]]$Candidates)
    foreach ($candidate in $Candidates) {
        if (Test-Path -LiteralPath $candidate) {
            return $candidate
        }
        $command = Get-Command $candidate -ErrorAction SilentlyContinue
        if ($null -ne $command) {
            return $command.Source
        }
    }
    throw "Required executable not found: $($Candidates -join ', ')"
}

$pythonExe = Resolve-Executable -Candidates $pythonCandidates
$rscriptExe = Resolve-Executable -Candidates $rCandidates
$resolvedOutput = [System.IO.Path]::GetFullPath($OutputDir)

if (Test-Path -LiteralPath $resolvedOutput) {
    $existing = Get-ChildItem -LiteralPath $resolvedOutput -Force -ErrorAction Stop
    if ($existing.Count -gt 0) {
        throw "Output directory is not empty. Use a fresh directory to preserve provenance: $resolvedOutput"
    }
} else {
    New-Item -ItemType Directory -Path $resolvedOutput -Force | Out-Null
}

Write-Output "[1/3] Retrieving and inspecting official supplement"
& $pythonExe (Join-Path $scriptRoot '99A1_download_and_inspect_supplement.py') --output $resolvedOutput
if ($LASTEXITCODE -ne 0) { throw "Source inspection failed with exit code $LASTEXITCODE" }

Write-Output "[2/3] Auditing frozen objects, joins, eligibility, and taxonomy coverage"
& $rscriptExe (Join-Path $scriptRoot '99A2_direction_feasibility_audit.R') --output $resolvedOutput
if ($LASTEXITCODE -ne 0) { throw "R feasibility audit failed with exit code $LASTEXITCODE" }

Write-Output "[3/3] Verifying input guards and finalizing manifest"
& $pythonExe (Join-Path $scriptRoot '99A3_finalize_feasibility.py') --output $resolvedOutput
if ($LASTEXITCODE -ne 0) { throw "Final quality control failed with exit code $LASTEXITCODE" }

Write-Output "A-stage feasibility audit completed: $resolvedOutput"
