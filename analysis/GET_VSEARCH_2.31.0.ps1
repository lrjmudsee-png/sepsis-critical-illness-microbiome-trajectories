$ErrorActionPreference = "Stop"

$Root = "E:\sepsis_project"
$Tools = Join-Path $Root "tools"
$Version = "2.31.0"
$Folder = Join-Path $Tools "vsearch-$Version-win-x86_64"
$Exe = Join-Path $Folder "bin\vsearch.exe"

if (Test-Path $Exe) {
    Write-Host "VSEARCH already installed:"
    Write-Host $Exe
    & $Exe --version
    exit 0
}

New-Item -ItemType Directory -Force -Path $Tools | Out-Null

$Zip = Join-Path $Tools "vsearch-$Version-win-x86_64.zip"
$Url = "https://github.com/torognes/vsearch/releases/download/v$Version/vsearch-$Version-win-x86_64.zip"

Write-Host "Downloading official VSEARCH $Version Windows binary..."
Invoke-WebRequest -Uri $Url -OutFile $Zip

Write-Host "Extracting..."
Expand-Archive -LiteralPath $Zip -DestinationPath $Tools -Force

if (-not (Test-Path $Exe)) {
    throw "VSEARCH extraction completed but vsearch.exe was not found at $Exe"
}

Write-Host "VSEARCH installed:"
Write-Host $Exe
& $Exe --version
