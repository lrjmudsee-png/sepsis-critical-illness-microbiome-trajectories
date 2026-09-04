$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\95A_global_evidence_synthesis_and_manuscript_map.R"

& $Rscript $Script
exit $LASTEXITCODE
