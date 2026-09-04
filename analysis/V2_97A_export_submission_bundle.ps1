[CmdletBinding()]
param(
    [string]$ProjectRoot = "",
    [string]$OutputDirectory = ""
)

$ErrorActionPreference = "Stop"

if ([string]::IsNullOrWhiteSpace($ProjectRoot)) {
    $ProjectRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot "..\.."))
}

if ([string]::IsNullOrWhiteSpace($OutputDirectory)) {
    $OutputDirectory = Join-Path $ProjectRoot "V2-SUBMISSION\04_可上传压缩包"
}

$sourceDirectory = Join-Path $ProjectRoot "V2-SUBMISSION\01_投稿用_最终文件"
if (-not (Test-Path -LiteralPath $sourceDirectory -PathType Container)) {
    throw "Submission source directory not found: $sourceDirectory"
}

New-Item -ItemType Directory -Force -Path $OutputDirectory | Out-Null

$timestamp = Get-Date -Format "yyyyMMdd_HHmmss"
$zipPath = Join-Path $OutputDirectory "Sepsis_V2_Submission_$timestamp.zip"
$manifestPath = Join-Path $OutputDirectory "Sepsis_V2_Submission_$timestamp.sha256.csv"

$files = Get-ChildItem -LiteralPath $sourceDirectory -Recurse -File | Sort-Object FullName
$manifest = foreach ($file in $files) {
    $relativePath = [System.IO.Path]::GetRelativePath($sourceDirectory, $file.FullName)
    $hash = Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256
    [pscustomobject]@{
        RelativePath = $relativePath.Replace("\", "/")
        Bytes = $file.Length
        SHA256 = $hash.Hash
    }
}

$manifest | Export-Csv -LiteralPath $manifestPath -NoTypeInformation -Encoding utf8
Compress-Archive -Path (Join-Path $sourceDirectory "*") -DestinationPath $zipPath -CompressionLevel Optimal -Force

$zipHash = Get-FileHash -LiteralPath $zipPath -Algorithm SHA256
Write-Output "Submission bundle created: $zipPath"
Write-Output "Manifest created: $manifestPath"
Write-Output "ZIP SHA256: $($zipHash.Hash)"
