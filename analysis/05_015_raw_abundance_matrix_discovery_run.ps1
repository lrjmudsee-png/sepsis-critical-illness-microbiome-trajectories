$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_015_raw_abundance_matrix_discovery.R')

& $Rscript $Script

exit $LASTEXITCODE
