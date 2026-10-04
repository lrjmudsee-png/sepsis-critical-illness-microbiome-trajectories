
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_022_sample_mapping_build.R')

& $Rscript $Script

exit $LASTEXITCODE
