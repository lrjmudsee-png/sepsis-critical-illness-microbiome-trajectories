$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
$Script="E:\sepsis_project\code\03_data_processing\94B2_CRA002354_formal_source_longitudinal_analysis.R"

& $Rscript $Script

exit $LASTEXITCODE
