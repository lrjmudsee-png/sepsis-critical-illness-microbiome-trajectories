<#
Sepsis V2 - Step 01
Build project-level data inventory.
ASCII-only for Windows PowerShell 5.1.

Scans:
- F:\sepsis_V2\02_raw_data
- F:\sepsis
- E:\sepsis_project for project-associated metadata/processed files

Outputs:
- V2_project_inventory_*.csv
- V2_metadata_candidates_*.csv
- V2_processed_candidates_*.csv
- V2_inventory_summary_*.txt

This script does not modify raw data.
#>

[CmdletBinding()]
param(
    [string[]]$RawRoots = @(
        "F:\sepsis_V2\02_raw_data",
        "F:\sepsis"
    ),
    [string]$CodeRoot = "E:\sepsis_project",
    [string]$OutputDir = "E:\sepsis_project\results\V2_00_inventory"
)

$ErrorActionPreference = "Stop"
New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null
$timestamp = Get-Date -Format "yyyyMMdd_HHmmss"

function Get-ProjectId {
    param([string]$Name)

    if ($Name -match '^(PRJ[A-Z0-9]+)$') { return $Matches[1] }
    if ($Name -match '^(CRA[0-9]+)$') { return $Matches[1] }
    if ($Name -match '^(ERP[0-9]+|SRP[0-9]+|DRP[0-9]+)$') { return $Matches[1] }

    return $null
}

function Get-FileCount {
    param([string]$Path,[string[]]$Patterns)

    if (-not (Test-Path $Path)) { return 0 }

    $count = 0
    foreach ($pattern in $Patterns) {
        $count += @(
            Get-ChildItem -Path $Path -Recurse -File -Filter $pattern -ErrorAction SilentlyContinue
        ).Count
    }
    return $count
}

function Get-TotalBytes {
    param([string]$Path,[string[]]$Patterns)

    if (-not (Test-Path $Path)) { return [int64]0 }

    [int64]$sum = 0
    foreach ($pattern in $Patterns) {
        $files = @(
            Get-ChildItem -Path $Path -Recurse -File -Filter $pattern -ErrorAction SilentlyContinue
        )
        foreach ($f in $files) {
            $sum += [int64]$f.Length
        }
    }
    return $sum
}

function Test-MetadataCandidate {
    param([System.IO.FileInfo]$File)

    $name = $File.Name.ToLowerInvariant()
    $ext = $File.Extension.ToLowerInvariant()

    if ($ext -notin @(".csv",".tsv",".txt",".xlsx",".xls",".json",".xml")) {
        return $false
    }

    return (
        $name -match "metadata" -or
        $name -match "manifest" -or
        $name -match "sample" -or
        $name -match "biosample" -or
        $name -match "clinical" -or
        $name -match "phenotype" -or
        $name -match "patient" -or
        $name -match "subject" -or
        $name -match "mapping" -or
        $name -match "runinfo" -or
        $name -match "run_table" -or
        $name -match "ena_run"
    )
}

function Test-ProcessedCandidate {
    param([System.IO.FileInfo]$File)

    $name = $File.Name.ToLowerInvariant()
    $ext = $File.Extension.ToLowerInvariant()

    if ($ext -notin @(".csv",".tsv",".txt",".rds",".rdata",".xlsx",".xls",".biom")) {
        return $false
    }

    return (
        $name -match "asv" -or
        $name -match "otu" -or
        $name -match "feature" -or
        $name -match "taxonomy" -or
        $name -match "taxa" -or
        $name -match "genus" -or
        $name -match "phyloseq" -or
        $name -match "count_table" -or
        $name -match "abundance"
    )
}

# Discover project directories
$projectDirs = @()

foreach ($root in $RawRoots) {
    if (-not (Test-Path $root)) {
        Write-Host "Skip missing root: $root"
        continue
    }

    Write-Host "Scanning root: $root"

    foreach ($dir in @(Get-ChildItem -Path $root -Directory -ErrorAction SilentlyContinue)) {
        $project = Get-ProjectId -Name $dir.Name
        if ($project) {
            $projectDirs += [pscustomobject]@{
                Project = $project
                ProjectPath = $dir.FullName
                Root = $root
            }
        }
    }
}

$projectDirs = @($projectDirs | Sort-Object ProjectPath -Unique)

if ($projectDirs.Count -eq 0) {
    throw "No project directories were found."
}

$inventory = @()
$metadataCandidates = @()
$processedCandidates = @()

foreach ($p in $projectDirs) {
    $project = $p.Project
    $path = $p.ProjectPath

    Write-Host "Inventory: $project"

    $fastqGzCount = Get-FileCount -Path $path -Patterns @("*.fastq.gz","*.fq.gz")
    $fastqPlainCount = Get-FileCount -Path $path -Patterns @("*.fastq","*.fq")
    $fastaGzCount = Get-FileCount -Path $path -Patterns @("*.fa.gz","*.fasta.gz")
    $fastaPlainCount = Get-FileCount -Path $path -Patterns @("*.fa","*.fasta")
    $sraCount = Get-FileCount -Path $path -Patterns @("*.sra")
    $aria2Count = Get-FileCount -Path $path -Patterns @("*.aria2")

    $rawBytes = Get-TotalBytes -Path $path -Patterns @(
        "*.fastq.gz","*.fq.gz","*.fastq","*.fq",
        "*.fa.gz","*.fasta.gz","*.fa","*.fasta","*.sra"
    )

    $allFiles = @(Get-ChildItem -Path $path -Recurse -File -ErrorAction SilentlyContinue)
    $projectMetadata = @($allFiles | Where-Object { Test-MetadataCandidate -File $_ })
    $projectProcessed = @($allFiles | Where-Object { Test-ProcessedCandidate -File $_ })

    foreach ($f in $projectMetadata) {
        $metadataCandidates += [pscustomobject]@{
            Project = $project
            FileName = $f.Name
            Extension = $f.Extension
            SizeKB = [math]::Round($f.Length / 1KB, 2)
            LastWriteTime = $f.LastWriteTime
            FullPath = $f.FullName
        }
    }

    foreach ($f in $projectProcessed) {
        $processedCandidates += [pscustomobject]@{
            Project = $project
            FileName = $f.Name
            Extension = $f.Extension
            SizeKB = [math]::Round($f.Length / 1KB, 2)
            LastWriteTime = $f.LastWriteTime
            FullPath = $f.FullName
        }
    }

    $rawStatus = "NO_RAW"
    $rawCount = $fastqGzCount + $fastqPlainCount + $fastaGzCount + $fastaPlainCount + $sraCount

    if ($rawCount -gt 0) { $rawStatus = "RAW_PRESENT" }
    if ($aria2Count -gt 0) { $rawStatus = "DOWNLOADING_OR_PARTIAL" }

    $inventory += [pscustomobject]@{
        Project = $project
        Root = $p.Root
        ProjectPath = $path
        FastqGz = $fastqGzCount
        FastqPlain = $fastqPlainCount
        FastaGz = $fastaGzCount
        FastaPlain = $fastaPlainCount
        SraFiles = $sraCount
        Aria2Partial = $aria2Count
        RawDataGB = [math]::Round($rawBytes / 1GB, 3)
        MetadataCandidates = $projectMetadata.Count
        ProcessedCandidates = $projectProcessed.Count
        RawStatus = $rawStatus
        HasMetadataCandidate = if ($projectMetadata.Count -gt 0) { "YES" } else { "NO" }
        HasProcessedCandidate = if ($projectProcessed.Count -gt 0) { "YES" } else { "NO" }
    }
}

# Scan CodeRoot for project-associated small files
if (Test-Path $CodeRoot) {
    Write-Host "Scanning CodeRoot for project-associated table files..."

    $smallFiles = @(
        Get-ChildItem -Path $CodeRoot -Recurse -File -ErrorAction SilentlyContinue |
        Where-Object {
            $_.Extension.ToLowerInvariant() -in @(
                ".csv",".tsv",".txt",".xlsx",".xls",".rds",".rdata",".biom"
            )
        }
    )

    $knownProjects = @($inventory.Project | Sort-Object -Unique)

    foreach ($f in $smallFiles) {
        $searchText = $f.FullName

        foreach ($project in $knownProjects) {
            if ($searchText -match [regex]::Escape($project)) {

                if (Test-MetadataCandidate -File $f) {
                    $metadataCandidates += [pscustomobject]@{
                        Project = $project
                        FileName = $f.Name
                        Extension = $f.Extension
                        SizeKB = [math]::Round($f.Length / 1KB, 2)
                        LastWriteTime = $f.LastWriteTime
                        FullPath = $f.FullName
                    }
                }

                if (Test-ProcessedCandidate -File $f) {
                    $processedCandidates += [pscustomobject]@{
                        Project = $project
                        FileName = $f.Name
                        Extension = $f.Extension
                        SizeKB = [math]::Round($f.Length / 1KB, 2)
                        LastWriteTime = $f.LastWriteTime
                        FullPath = $f.FullName
                    }
                }

                break
            }
        }
    }
}

$metadataCandidates = @($metadataCandidates | Sort-Object Project,FullPath -Unique)
$processedCandidates = @($processedCandidates | Sort-Object Project,FullPath -Unique)

$inventoryPath = Join-Path $OutputDir "V2_project_inventory_$timestamp.csv"
$metadataPath = Join-Path $OutputDir "V2_metadata_candidates_$timestamp.csv"
$processedPath = Join-Path $OutputDir "V2_processed_candidates_$timestamp.csv"
$summaryPath = Join-Path $OutputDir "V2_inventory_summary_$timestamp.txt"

$inventory |
    Sort-Object Project |
    Export-Csv -Path $inventoryPath -NoTypeInformation -Encoding UTF8

$metadataCandidates |
    Export-Csv -Path $metadataPath -NoTypeInformation -Encoding UTF8

$processedCandidates |
    Export-Csv -Path $processedPath -NoTypeInformation -Encoding UTF8

$summaryLines = @(
    "Sepsis V2 data inventory",
    "Generated: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')",
    "",
    "Project count: $(@($inventory).Count)",
    "Projects with raw data: $(@($inventory | Where-Object {$_.RawStatus -ne 'NO_RAW'}).Count)",
    "Projects with partial aria2 downloads: $(@($inventory | Where-Object {$_.Aria2Partial -gt 0}).Count)",
    "Projects with metadata candidates: $(@($inventory | Where-Object {$_.HasMetadataCandidate -eq 'YES'}).Count)",
    "Projects with processed-data candidates: $(@($inventory | Where-Object {$_.HasProcessedCandidate -eq 'YES'}).Count)",
    "",
    "Inventory CSV: $inventoryPath",
    "Metadata candidates CSV: $metadataPath",
    "Processed candidates CSV: $processedPath"
)

$summaryLines | Set-Content -Path $summaryPath -Encoding UTF8

Write-Host ""
Write-Host "================ V2 STEP 01 COMPLETE ================"

$inventory |
    Sort-Object Project |
    Format-Table Project,FastqGz,FastqPlain,FastaGz,SraFiles,Aria2Partial,RawDataGB,MetadataCandidates,ProcessedCandidates,RawStatus -AutoSize

Write-Host ""
Write-Host "Inventory:  $inventoryPath"
Write-Host "Metadata:   $metadataPath"
Write-Host "Processed:  $processedPath"
Write-Host "Summary:    $summaryPath"
Write-Host "====================================================="
