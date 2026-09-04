$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script="E:\sepsis_project\code\03_data_processing\93C_taxonomic_divergence_vs_EII.R"

& $Rscript $Script

exit $LASTEXITCODE
