$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '09_003_manuscript_source_pack_assembly.R')

& $Rscript $Script
exit $LASTEXITCODE
