$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '06_008_cra002354_nonchimeric_otu_table_and_rarefaction_fix.R')

& $Rscript $Script

exit $LASTEXITCODE
