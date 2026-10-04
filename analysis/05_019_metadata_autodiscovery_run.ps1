
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_019_metadata_autodiscovery.R')

& $Rscript $Script

exit $LASTEXITCODE
