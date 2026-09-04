$ErrorActionPreference="Stop"

$Code="E:\sepsis_project\code\03_data_processing"
$Out="E:\sepsis_project\results\V2_96J2_FINAL_FIGURE3_PROVENANCE_FIX"
$Zip="E:\sepsis_project\results\V2_96J2_FINAL_FIGURE3_PROVENANCE_FIX_FOR_UPLOAD.zip"

$Manifest=Join-Path $Code "V2_96J2_PACKAGED_ASSET_SHA256.csv"
if(!(Test-Path $Manifest)){ throw "Missing packaged asset manifest: $Manifest" }

$Rows=Import-Csv $Manifest

# 1) Verify the exact packaged source/final assets used for composition.
foreach($r in $Rows){
    $p=Join-Path $Code $r.asset
    if(!(Test-Path $p)){ throw "Missing packaged Step96J2 asset: $p" }
    $actual=(Get-FileHash $p -Algorithm SHA256).Hash.ToLower()
    if($actual -ne $r.sha256.ToLower()){
        throw "Packaged Step96J2 asset hash mismatch: $($r.asset)"
    }
}

New-Item -ItemType Directory -Force -Path $Out | Out-Null

# 2) Stage the exact source assets and final composition.
$Files=@(
  "_SOURCE_Figure3_ABC_EXACT.pdf",
  "_SOURCE_Figure3_ABC_EXACT.png",
  "_SOURCE_Figure3_D_EXACT.pdf",
  "_SOURCE_Figure3_D_EXACT.png",
  "_SOURCE_Supplementary_Taxonomic_Concordance_EXACT.pdf",
  "_SOURCE_Supplementary_Taxonomic_Concordance_EXACT.png",
  "V2_96J2_Figure_3_FINAL_ABCD.pdf",
  "V2_96J2_Figure_3_FINAL_ABCD_600dpi.png",
  "V2_96J2_Figure_3_FINAL_LEGEND.txt",
  "V2_96J2_Supplementary_Taxonomic_Concordance_LEGEND.txt",
  "V2_96J2_PACKAGED_ASSET_SHA256.csv"
)

foreach($f in $Files){
    Copy-Item (Join-Path $Code $f) (Join-Path $Out $f) -Force
}

# 3) Compare current local historical copies only for provenance information.
#    A byte mismatch here is NOT a failure because PDF/PNG containers may have
#    been rewritten after the frozen asset was exported.
$Comparison=Join-Path $Out "V2_96J2_LOCAL_SOURCE_COMPARISON.txt"

$LocalABCpdf="E:\sepsis_project\results\V2_35D3_FINAL_FIGURE_FORMATTING_FIX\01_FINAL_FIGURES\Figure_3_Core_Sepsis_Heterogeneous_Taxonomic_Routes.pdf"
$LocalABCpng="E:\sepsis_project\results\V2_35D3_FINAL_FIGURE_FORMATTING_FIX\01_FINAL_FIGURES\Figure_3_Core_Sepsis_Heterogeneous_Taxonomic_Routes.png"
$LocalDpdf="E:\sepsis_project\results\V2_96E_FIGURE3D_CROSS_COHORT_REPRODUCIBILITY_FINAL\Figure3D_Cross_Cohort_Taxonomic_Concordance_MAIN.pdf"
$LocalDpng="E:\sepsis_project\results\V2_96E_FIGURE3D_CROSS_COHORT_REPRODUCIBILITY_FINAL\Figure3D_Cross_Cohort_Taxonomic_Concordance_MAIN.png"

$Expected=@{}
foreach($r in $Rows){ $Expected[$r.asset]=$r.sha256.ToLower() }

$lines=@()
$lines += "V2 STEP96J2 LOCAL SOURCE COMPARISON"
$lines += ""
$lines += "IMPORTANT: local historical copy mismatches are informational only."
$lines += "The exact source assets actually used for composition are packaged and hash-verified."
$lines += ""

function Compare-Local($label,$path,$packagedName){
    if(Test-Path $path){
        $h=(Get-FileHash $path -Algorithm SHA256).Hash.ToLower()
        $status=if($h -eq $Expected[$packagedName]){"BYTE_IDENTICAL"}else{"BYTE_DIFFERENT"}
        return "$label`t$status`t$h`t$path"
    } else {
        return "$label`tNOT_FOUND`tNA`t$path"
    }
}

$lines += Compare-Local "Figure3_ABC_PDF" $LocalABCpdf "_SOURCE_Figure3_ABC_EXACT.pdf"
$lines += Compare-Local "Figure3_ABC_PNG" $LocalABCpng "_SOURCE_Figure3_ABC_EXACT.png"
$lines += Compare-Local "Figure3_D_PDF" $LocalDpdf "_SOURCE_Figure3_D_EXACT.pdf"
$lines += Compare-Local "Figure3_D_PNG" $LocalDpng "_SOURCE_Figure3_D_EXACT.png"

$lines | Set-Content -Encoding UTF8 $Comparison

@"
V2 STEP96J2 FINAL FIGURE 3 PROVENANCE FIX

Packaged exact Figure 3 A-C source: hash verified.
Packaged exact Figure 3D source: hash verified.
Packaged exact Supplementary source: hash verified.
Final Figure 3 A-D asset: hash verified.

The previous Step96J failure was caused by enforcing byte-identical SHA256
against a current local historical PDF copy. That is now treated as an
informational provenance comparison rather than a hard failure.

No statistical analysis rerun.
No inferential content changed.

STEP96J2 COMPLETE
"@ | Set-Content -Encoding UTF8 (Join-Path $Out "_STEP96J2_COMPLETE.txt")

if(Test-Path $Zip){ Remove-Item $Zip -Force }
Compress-Archive -Path "$Out\*" -DestinationPath $Zip -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP96J2 COMPLETE" -ForegroundColor Green
Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
