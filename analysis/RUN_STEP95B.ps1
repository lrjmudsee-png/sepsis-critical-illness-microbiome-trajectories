$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\95B_manuscript_source_pack_assembly.R"

& $Rscript $Script
exit $LASTEXITCODE
