
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script=(Join-Path $PSScriptRoot '05_031_eii_file_selection_fix.R')

& $Rscript $Script

exit $LASTEXITCODE
