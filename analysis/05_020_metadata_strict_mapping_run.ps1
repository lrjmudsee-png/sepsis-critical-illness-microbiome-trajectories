
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_020_metadata_strict_mapping.R')

& $Rscript $Script

exit $LASTEXITCODE
