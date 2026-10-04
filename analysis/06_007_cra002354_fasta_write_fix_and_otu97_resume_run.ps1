$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script=(Join-Path $PSScriptRoot '06_007_cra002354_fasta_write_fix_and_otu97_resume.R')

& $Rscript $Script

exit $LASTEXITCODE
