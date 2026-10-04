$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '06_009_cra002354_taxonomy_id_repair.R')

& $Rscript $Script
exit $LASTEXITCODE
