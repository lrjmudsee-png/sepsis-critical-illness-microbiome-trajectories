$ErrorActionPreference="Stop"

$Root="E:\sepsis_project"
$B1B=Join-Path $Root "results\V2_34B1B_CRA002354_FASTA_WRITE_FIX_AND_OTU97_RESUME"
$Data=Join-Path $Root "data\CRA002354\03_vsearch97_silva1382"
$Out=Join-Path $Root "results\V2_34B1C_B1B_STATUS_AND_RECOVERY_AUDIT"

New-Item -ItemType Directory -Force -Path $Out | Out-Null

$rows = @()

function Add-FileStatus {
    param(
        [string]$Label,
        [string]$Path
    )

    $exists = Test-Path $Path
    $size = $null
    $mtime = $null

    if ($exists) {
        $item = Get-Item $Path
        if (-not $item.PSIsContainer) {
            $size = $item.Length
            $mtime = $item.LastWriteTime
        }
    }

    $script:rows += [pscustomobject]@{
        label = $Label
        path = $Path
        exists = $exists
        size_bytes = $size
        last_write_time = $mtime
    }
}

Add-FileStatus "B1B_complete_flag" (Join-Path $B1B "_STEP94B1B_COMPLETE.ok")
Add-FileStatus "B1B_readiness" (Join-Path $B1B "10_STEP94B1B_ANALYSIS_READINESS.csv")
Add-FileStatus "B1B_interpretation" (Join-Path $B1B "11_STEP94B1B_INTERPRETATION.txt")

Add-FileStatus "pooled_clean_fasta" (Join-Path $Data "CRA002354_publication_aligned_clean_pooled_WIDTH80.fa")
Add-FileStatus "unique_fasta" (Join-Path $Data "CRA002354_clean_uniques.fa")
Add-FileStatus "nonchimera_fasta" (Join-Path $Data "CRA002354_clean_uniques_nonchimeric.fa")
Add-FileStatus "OTU97_centroids" (Join-Path $Data "CRA002354_OTU97_centroids.fa")
Add-FileStatus "OTU97_table" (Join-Path $Data "CRA002354_OTU97_table.tsv")
Add-FileStatus "taxonomy" (Join-Path $Data "CRA002354_OTU97_SILVA1382_taxonomy_minBoot80.csv")
Add-FileStatus "genus_relative_full" (Join-Path $Data "CRA002354_genus_relative_abundance_OTU97_SILVA1382.csv")
Add-FileStatus "genus_relative_primary4000" (Join-Path $Data "CRA002354_genus_relative_abundance_PRIMARY_min4000.csv")
Add-FileStatus "metadata_primary4000" (Join-Path $Data "CRA002354_analysis_metadata_PRIMARY_min4000.csv")
Add-FileStatus "genus_relative_sensitivity2000" (Join-Path $Data "CRA002354_genus_relative_abundance_SENSITIVITY_min2000.csv")
Add-FileStatus "metadata_sensitivity2000" (Join-Path $Data "CRA002354_analysis_metadata_SENSITIVITY_min2000.csv")

$rows | Export-Csv -NoTypeInformation -Encoding UTF8 (Join-Path $Out "01_B1B_FILE_STATUS.csv")

# VSEARCH process status
$vsearch = Get-Process vsearch -ErrorAction SilentlyContinue

if ($null -ne $vsearch) {
    $vsearch |
        Select-Object Id,CPU,WorkingSet64,PrivateMemorySize64,StartTime |
        Export-Csv -NoTypeInformation -Encoding UTF8 (Join-Path $Out "02_VSEARCH_PROCESS_STATUS.csv")
} else {
    [pscustomobject]@{
        Id=""
        CPU=""
        WorkingSet64=""
        PrivateMemorySize64=""
        StartTime=""
    } | Export-Csv -NoTypeInformation -Encoding UTF8 (Join-Path $Out "02_VSEARCH_PROCESS_STATUS.csv")
}

# Detect current completed stage from output files.
$stage = "NOT_STARTED_OR_UNKNOWN"

$pooled = Join-Path $Data "CRA002354_publication_aligned_clean_pooled_WIDTH80.fa"
$uniq = Join-Path $Data "CRA002354_clean_uniques.fa"
$nonchim = Join-Path $Data "CRA002354_clean_uniques_nonchimeric.fa"
$otus = Join-Path $Data "CRA002354_OTU97_centroids.fa"
$otutab = Join-Path $Data "CRA002354_OTU97_table.tsv"
$tax = Join-Path $Data "CRA002354_OTU97_SILVA1382_taxonomy_minBoot80.csv"
$gen = Join-Path $Data "CRA002354_genus_relative_abundance_PRIMARY_min4000.csv"

if ((Test-Path $pooled) -and ((Get-Item $pooled).Length -gt 500MB)) {
    $stage = "POOLED_FASTA_READY"
}
if ((Test-Path $uniq) -and ((Get-Item $uniq).Length -gt 1KB)) {
    $stage = "DEREPLICATION_COMPLETE"
}
if ((Test-Path $nonchim) -and ((Get-Item $nonchim).Length -gt 1KB)) {
    $stage = "CHIMERA_REMOVAL_COMPLETE"
}
if ((Test-Path $otus) -and ((Get-Item $otus).Length -gt 1KB)) {
    $stage = "OTU97_CLUSTERING_COMPLETE"
}
if ((Test-Path $otutab) -and ((Get-Item $otutab).Length -gt 100)) {
    $stage = "OTU_TABLE_COMPLETE"
}
if ((Test-Path $tax) -and ((Get-Item $tax).Length -gt 100)) {
    $stage = "TAXONOMY_COMPLETE"
}
if ((Test-Path $gen) -and ((Get-Item $gen).Length -gt 100)) {
    $stage = "B1B_DATA_PRODUCTS_COMPLETE"
}

$readyFile = Join-Path $B1B "10_STEP94B1B_ANALYSIS_READINESS.csv"
$readyValue = "UNKNOWN"

if (Test-Path $readyFile) {
    try {
        $r = Import-Csv $readyFile
        if ($r.Count -gt 0 -and $null -ne $r[0].ready_for_step94B2) {
            $readyValue = $r[0].ready_for_step94B2
        }
    } catch {
        $readyValue = "READ_ERROR"
    }
}

$completeFlag = Test-Path (Join-Path $B1B "_STEP94B1B_COMPLETE.ok")
$vsearchRunning = $null -ne $vsearch

$decision = if ($completeFlag -and ($readyValue -match "TRUE|True|true")) {
    "B1B_COMPLETE_AND_READY_FOR_B2"
} elseif ($vsearchRunning) {
    "B1B_STILL_RUNNING_DO_NOT_RESTART"
} elseif ($stage -eq "CHIMERA_REMOVAL_COMPLETE" -or
         $stage -eq "OTU97_CLUSTERING_COMPLETE" -or
         $stage -eq "OTU_TABLE_COMPLETE" -or
         $stage -eq "TAXONOMY_COMPLETE" -or
         $stage -eq "B1B_DATA_PRODUCTS_COMPLETE") {
    "B1B_PARTIAL_PRODUCTS_EXIST_PROCESS_NOT_RUNNING_REVIEW_BEFORE_RECOVERY"
} elseif ($stage -eq "DEREPLICATION_COMPLETE") {
    "B1B_STOPPED_AT_OR_DURING_CHIMERA_REMOVAL"
} else {
    "B1B_INCOMPLETE_OR_NOT_STARTED"
}

[pscustomobject]@{
    detected_stage = $stage
    vsearch_running = $vsearchRunning
    B1B_complete_flag = $completeFlag
    ready_for_step94B2 = $readyValue
    decision = $decision
} | Export-Csv -NoTypeInformation -Encoding UTF8 (Join-Path $Out "03_B1B_DECISION.csv")

# Copy readable terminal logs only if not locked.
$logPatterns = @(
    "VSEARCH_01_DEREPLICATION.stderr.txt",
    "VSEARCH_02_CHIMERA_REMOVAL.stderr.txt",
    "VSEARCH_03_OTU97_CLUSTERING.stderr.txt",
    "VSEARCH_04_MAP_READS_TO_OTUS.stderr.txt"
)

foreach ($lf in $logPatterns) {
    $src = Join-Path $B1B $lf
    if (Test-Path $src) {
        try {
            Get-Content $src -Tail 80 |
                Out-File -Encoding UTF8 (Join-Path $Out ("TAIL_" + $lf))
        } catch {
            ("LOCKED_OR_UNREADABLE: " + $_.Exception.Message) |
                Out-File -Encoding UTF8 (Join-Path $Out ("TAIL_" + $lf))
        }
    }
}

$summary = @(
    "STEP94B1C B1B STATUS AUDIT"
    ""
    ("Detected stage: " + $stage)
    ("VSEARCH running: " + $vsearchRunning)
    ("B1B complete flag: " + $completeFlag)
    ("ready_for_step94B2: " + $readyValue)
    ("Decision: " + $decision)
    ""
)

if ($decision -eq "B1B_COMPLETE_AND_READY_FOR_B2") {
    $summary += "NEXT: run Step94B2."
    $summary += 'powershell -ExecutionPolicy Bypass -File "E:\sepsis_project\code\03_data_processing\RUN_STEP94B2.ps1"'
} elseif ($decision -eq "B1B_STILL_RUNNING_DO_NOT_RESTART") {
    $summary += "NEXT: let the current B1B/VSEARCH process continue. Do NOT launch a second B1B."
} else {
    $summary += "NEXT: review this audit before restarting or recovering B1B."
}

$summary | Out-File -Encoding UTF8 (Join-Path $Out "04_STEP94B1C_INTERPRETATION.txt")

"STEP94B1C COMPLETE" | Out-File -Encoding UTF8 (Join-Path $Out "_STEP94B1C_COMPLETE.ok")

Write-Host ""
Write-Host "Detected stage:" $stage
Write-Host "VSEARCH running:" $vsearchRunning
Write-Host "B1B complete flag:" $completeFlag
Write-Host "ready_for_step94B2:" $readyValue
Write-Host "Decision:" $decision
