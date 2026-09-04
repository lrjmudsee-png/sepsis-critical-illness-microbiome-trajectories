
$ErrorActionPreference="Stop"

$Rscript="C:\Program Files\R\R-4.4.0\bin\Rscript.exe"

$Script="E:\sepsis_project\code\03_data_processing\93O_taxonomic_abundance_metadata_integration.R"

& $Rscript $Script

exit $LASTEXITCODE
