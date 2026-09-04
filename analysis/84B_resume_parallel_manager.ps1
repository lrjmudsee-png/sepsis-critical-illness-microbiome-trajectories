param(
    [ValidateSet('Status','Run','Launch','Collect','Unlock')]
    [string]$Mode = 'Status',

    [ValidateSet('PRJEB82425','PRJNA516701','PRJNA578267','PRJNA851469','PRJNA1166732')]
    [string]$Project,

    [int]$MaxParallel = 2,

    [switch]$Force
)

$ErrorActionPreference = 'Stop'

$Rscript = 'C:\Program Files\R\R-4.4.0\bin\Rscript.exe'
$Rscript84B = 'E:\sepsis_project\code\03_data_processing\84B_run_full_cohort_DADA2.R'
$Root = 'E:\sepsis_project'
$Canon = Join-Path $Root 'data\_V2_ANALYSIS_READY'
$ResultRoot = Join-Path $Root 'results\V2_24B_FULL_DADA2'
$ManagerRoot = Join-Path $Root 'results\V2_24B_RESUME_MANAGER'
New-Item -ItemType Directory -Force -Path $ManagerRoot | Out-Null

$Projects = @(
    [pscustomobject]@{ Project='PRJNA578267'; ExpectedRuns=207 },
    [pscustomobject]@{ Project='PRJNA516701'; ExpectedRuns=171 },
    [pscustomobject]@{ Project='PRJNA1166732'; ExpectedRuns=120 },
    [pscustomobject]@{ Project='PRJNA851469'; ExpectedRuns=119 },
    [pscustomobject]@{ Project='PRJEB82425'; ExpectedRuns=92 }
)

function Get-ProjectPaths {
    param([string]$P)

    $Work = Join-Path $Canon "01_PROJECTS\$P\05_work\step84B_full"
    $Filtered = Join-Path $Work 'filtered'

    [pscustomobject]@{
        Work = $Work
        Filtered = $Filtered
        Complete = Join-Path $Work '_STEP84B_FULL_COMPLETE.ok'
        Review = Join-Path $Work '_STEP84B_FULL_QC_REVIEW_REQUIRED.txt'
        Lock = Join-Path $Work '_STEP84B_PARALLEL_RUNNING.lock'
        ErrF = Join-Path $Work "${P}_full_errF.rds"
        ErrR = Join-Path $Work "${P}_full_errR.rds"
        Seqtab = Join-Path $Work "${P}_seqtab_prechim_step84B.rds"
        Nochim = Join-Path $Work "${P}_seqtab_nochim_step84B.rds"
        Summary = Join-Path $ResultRoot "${P}_full_summary.csv"
        Tracking = Join-Path $ResultRoot "${P}_full_tracking.csv"
        LengthQC = Join-Path $ResultRoot "${P}_full_ASV_length_QC.csv"
    }
}

function Get-ProjectStatus {
    param(
        [string]$P,
        [int]$ExpectedRuns
    )

    $x = Get-ProjectPaths $P

    $FilteredCount = 0
    if (Test-Path $x.Filtered) {
        $FilteredCount = @(Get-ChildItem -LiteralPath $x.Filtered -File -Filter '*.fastq.gz' -ErrorAction SilentlyContinue).Count
    }

    $SummaryPass = $false
    $SummaryStatus = ''
    if (Test-Path $x.Summary) {
        try {
            $s = Import-Csv -LiteralPath $x.Summary | Select-Object -First 1
            $SummaryStatus = $s.final_status
            if ($s.full_qc_pass -match '^(TRUE|True|true|1)$') {
                $SummaryPass = $true
            }
        } catch {}
    }

    if (Test-Path $x.Complete) {
        $Stage = 'COMPLETE_PASS'
    }
    elseif (Test-Path $x.Review) {
        $Stage = 'COMPLETE_QC_REVIEW'
    }
    elseif ((Test-Path $x.Summary) -and $SummaryPass -and (Test-Path $x.Nochim)) {
        $Stage = 'RECOVERABLE_NEAR_COMPLETE'
    }
    elseif (Test-Path $x.Nochim) {
        $Stage = 'SEQTAB_EXISTS_NO_FINAL_MARKER'
    }
    elseif ((Test-Path $x.ErrF) -and (Test-Path $x.ErrR)) {
        $Stage = 'ERROR_MODELS_EXIST_PARTIAL_RUN'
    }
    elseif ($FilteredCount -gt 0) {
        $Stage = 'FILTERED_FILES_EXIST_PARTIAL_RUN'
    }
    else {
        $Stage = 'NOT_STARTED'
    }

    [pscustomobject]@{
        Project = $P
        ExpectedRuns = $ExpectedRuns
        Stage = $Stage
        Locked = (Test-Path $x.Lock)
        FilteredFiles = $FilteredCount
        ExpectedFilteredFiles = 2 * $ExpectedRuns
        ErrF = (Test-Path $x.ErrF)
        ErrR = (Test-Path $x.ErrR)
        SeqtabPrechim = (Test-Path $x.Seqtab)
        SeqtabNochim = (Test-Path $x.Nochim)
        Summary = (Test-Path $x.Summary)
        SummaryPass = $SummaryPass
        SummaryStatus = $SummaryStatus
    }
}

function Get-AllStatus {
    $Rows = foreach ($p in $Projects) {
        Get-ProjectStatus -P $p.Project -ExpectedRuns $p.ExpectedRuns
    }
    return @($Rows)
}

function Save-And-ShowStatus {
    $Rows = Get-AllStatus
    $Path = Join-Path $ManagerRoot 'V2_STEP84B_resume_status.csv'
    $Rows | Export-Csv -LiteralPath $Path -NoTypeInformation -Encoding UTF8
    $Rows | Format-Table -AutoSize
    Write-Host "`nStatus CSV: $Path" -ForegroundColor Cyan
}

function Invoke-OneProject {
    param([string]$P)

    $meta = $Projects | Where-Object Project -eq $P | Select-Object -First 1
    if (-not $meta) { throw "Unknown project: $P" }

    $status = Get-ProjectStatus -P $P -ExpectedRuns $meta.ExpectedRuns
    $paths = Get-ProjectPaths $P

    if (-not $Force) {
        if ($status.Stage -eq 'COMPLETE_PASS') {
            Write-Host "$P already COMPLETE_PASS. Skipping." -ForegroundColor Green
            return
        }
        if ($status.Stage -eq 'COMPLETE_QC_REVIEW') {
            Write-Host "$P already completed but requires QC review. Skipping automatic rerun." -ForegroundColor Yellow
            Write-Host 'Use -Force only if you intentionally want to rerun this cohort.'
            return
        }
    }

    New-Item -ItemType Directory -Force -Path $paths.Work | Out-Null

    if (Test-Path $paths.Lock) {
        throw "$P is locked. Another terminal may be running it, or the previous shutdown left a stale lock. Run: .\84B_resume_parallel_manager.ps1 -Mode Unlock -Project $P"
    }

    New-Item -ItemType Directory -Path $paths.Lock | Out-Null
    @(
        "project=$P",
        "pid=$PID",
        "started=$(Get-Date -Format o)",
        "computer=$env:COMPUTERNAME"
    ) | Set-Content -LiteralPath (Join-Path $paths.Lock 'lock_info.txt') -Encoding UTF8

    try {
        Write-Host "Starting/resuming $P ..." -ForegroundColor Cyan
        Write-Host "Existing stage: $($status.Stage)"
        Write-Host 'Important: the original 84B script resumes at PROJECT level.'
        Write-Host 'Completed cohorts are skipped; an interrupted cohort is safely rebuilt from its beginning.'
        Write-Host ''

        & $Rscript $Rscript84B $P
        $Exit = $LASTEXITCODE

        if ($Exit -ne 0) {
            throw "Rscript returned exit code $Exit"
        }

        Write-Host "`n$P worker finished." -ForegroundColor Green
    }
    finally {
        if (Test-Path $paths.Lock) {
            Remove-Item -LiteralPath $paths.Lock -Recurse -Force
        }
    }
}

function Launch-Workers {
    if ($MaxParallel -lt 1) { throw 'MaxParallel must be >= 1.' }
    if ($MaxParallel -gt 3) {
        Write-Warning 'Running more than 3 simultaneous DADA2 cohorts on one Windows workstation is not recommended unless you have abundant RAM and fast NVMe storage.'
    }

    $status = Get-AllStatus
    $Pending = @(
        $status | Where-Object {
            $_.Stage -notin @('COMPLETE_PASS','COMPLETE_QC_REVIEW') -and -not $_.Locked
        }
    )

    if ($Pending.Count -eq 0) {
        Write-Host 'No runnable projects remain.' -ForegroundColor Green
        return
    }

    $Launch = @($Pending | Select-Object -First $MaxParallel)

    foreach ($item in $Launch) {
        $P = $item.Project
        $Manager = $PSCommandPath
        $cmd = "& '$Manager' -Mode Run -Project '$P'"

        Start-Process powershell.exe -ArgumentList @(
            '-NoProfile',
            '-ExecutionPolicy','Bypass',
            '-Command',$cmd
        )

        Write-Host "Launched $P" -ForegroundColor Green
    }

    Write-Host "`nStarted $($Launch.Count) worker window(s)." -ForegroundColor Cyan
    Write-Host 'When one or both finish, run Launch again to start the next pending cohort.'
}

function Collect-Results {
    $Rows = New-Object System.Collections.Generic.List[object]

    foreach ($p in $Projects) {
        $x = Get-ProjectPaths $p.Project
        if (Test-Path $x.Summary) {
            try {
                $s = Import-Csv -LiteralPath $x.Summary | Select-Object -First 1
                if ($s) { $Rows.Add($s) }
            } catch {
                Write-Warning "Could not read $($x.Summary)"
            }
        }
    }

    $Combined = Join-Path $ManagerRoot 'V2_STEP84B_COLLECTED_project_summaries.csv'
    $Rows | Export-Csv -LiteralPath $Combined -NoTypeInformation -Encoding UTF8

    $Status = Get-AllStatus
    $StatusPath = Join-Path $ManagerRoot 'V2_STEP84B_resume_status.csv'
    $Status | Export-Csv -LiteralPath $StatusPath -NoTypeInformation -Encoding UTF8

    $Pass = @($Rows | Where-Object { $_.full_qc_pass -match '^(TRUE|True|true|1)$' }).Count
    $Review = @($Rows | Where-Object { $_.full_qc_pass -notmatch '^(TRUE|True|true|1)$' }).Count

    [pscustomobject]@{
        ApprovedProjects = 5
        SummariesFound = $Rows.Count
        FullQCPass = $Pass
        FullQCReview = $Review
        CompleteMarkers = @($Status | Where-Object Stage -eq 'COMPLETE_PASS').Count
        LockedNow = @($Status | Where-Object Locked).Count
    } | Format-List

    Write-Host "Combined summary: $Combined" -ForegroundColor Cyan
    Write-Host "Status: $StatusPath" -ForegroundColor Cyan
}

switch ($Mode) {
    'Status' {
        Save-And-ShowStatus
    }
    'Run' {
        if (-not $Project) { throw '-Project is required for Mode Run.' }
        Invoke-OneProject -P $Project
    }
    'Launch' {
        Launch-Workers
    }
    'Collect' {
        Collect-Results
    }
    'Unlock' {
        if (-not $Project) { throw '-Project is required for Mode Unlock.' }
        $x = Get-ProjectPaths $Project
        if (Test-Path $x.Lock) {
            Remove-Item -LiteralPath $x.Lock -Recurse -Force
            Write-Host "Removed stale lock for $Project" -ForegroundColor Yellow
        } else {
            Write-Host "No lock exists for $Project"
        }
    }
}
