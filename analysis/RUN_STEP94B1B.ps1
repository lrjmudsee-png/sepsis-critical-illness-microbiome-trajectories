$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\94B1B_CRA002354_fasta_write_fix_and_otu97_resume.R"

& $Rscript $Script

exit $LASTEXITCODE
