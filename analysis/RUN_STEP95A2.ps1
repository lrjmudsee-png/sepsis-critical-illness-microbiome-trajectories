$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\95A2_global_evidence_synthesis_source_anchor_fix.R"

& $Rscript $Script
exit $LASTEXITCODE
