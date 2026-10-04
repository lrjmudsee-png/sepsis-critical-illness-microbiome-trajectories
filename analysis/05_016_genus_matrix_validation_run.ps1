$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_016_genus_matrix_validation.R')

& $Rscript $Script

exit $LASTEXITCODE
