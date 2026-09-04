$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script="E:\sepsis_project\code\03_data_processing\93D_taxonomic_heterogeneity_index.R"

& $Rscript $Script

exit $LASTEXITCODE
