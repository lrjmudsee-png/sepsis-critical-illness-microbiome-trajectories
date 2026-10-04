
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '04_006_metadata_enhancement_audit.R')

& $Rscript $Script
