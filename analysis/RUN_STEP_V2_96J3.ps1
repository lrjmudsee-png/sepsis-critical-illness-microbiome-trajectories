$ErrorActionPreference="Stop"

# Self-contained package:
# always use the folder containing THIS script.
$PackageRoot=$PSScriptRoot

$Out="E:\sepsis_project\results\V2_96J3_FINAL_FIGURE3_SELF_CONTAINED"
$Zip="E:\sepsis_project\results\V2_96J3_FINAL_FIGURE3_SELF_CONTAINED_FOR_UPLOAD.zip"

$Manifest=Join-Path $PackageRoot "V2_96J2_PACKAGED_ASSET_SHA256.csv"

if(!(Test-Path $Manifest)){
    throw "Missing packaged asset manifest beside launcher: $Manifest"
}

$Rows=Import-Csv $Manifest

# ------------------------------------------------------------
# 1. Verify exact packaged assets used for final composition
# ------------------------------------------------------------
foreach($r in $Rows){

    $p=Join-Path $PackageRoot $r.asset

    if(!(Test-Path $p)){
        throw "Missing packaged verified asset: $p"
    }

    $actual=(Get-FileHash $p -Algorithm SHA256).Hash.ToLower()

    if($actual -ne $r.sha256.ToLower()){
        throw "Packaged asset hash mismatch: $($r.asset)"
    }
}

New-Item -ItemType Directory -Force -Path $Out | Out-Null

# ------------------------------------------------------------
# 2. Stage exact verified source assets + final Figure 3
# ------------------------------------------------------------
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

    $src=Join-Path $PackageRoot $f

    if(!(Test-Path $src)){
        throw "Missing required packaged asset: $src"
    }

    Copy-Item $src (Join-Path $Out $f) -Force
}

# ------------------------------------------------------------
# 3. Optional comparison with old local historical copies
#    Informational only; never a hard failure.
# ------------------------------------------------------------
$Comparison=Join-Path $Out "V2_96J3_LOCAL_HISTORICAL_COPY_COMPARISON.txt"

$LocalABCpdf="E:\sepsis_project\results\V2_35D3_FINAL_FIGURE_FORMATTING_FIX\01_FINAL_FIGURES\Figure_3_Core_Sepsis_Heterogeneous_Taxonomic_Routes.pdf"
$LocalABCpng="E:\sepsis_project\results\V2_35D3_FINAL_FIGURE_FORMATTING_FIX\01_FINAL_FIGURES\Figure_3_Core_Sepsis_Heterogeneous_Taxonomic_Routes.png"
$LocalDpdf="E:\sepsis_project\results\V2_96E_FIGURE3D_CROSS_COHORT_REPRODUCIBILITY_FINAL\Figure3D_Cross_Cohort_Taxonomic_Concordance_MAIN.pdf"
$LocalDpng="E:\sepsis_project\results\V2_96E_FIGURE3D_CROSS_COHORT_REPRODUCIBILITY_FINAL\Figure3D_Cross_Cohort_Taxonomic_Concordance_MAIN.png"

$Expected=@{}
foreach($r in $Rows){
    $Expected[$r.asset]=$r.sha256.ToLower()
}

function Compare-Local($label,$path,$packagedName){

    if(Test-Path $path){

        $h=(Get-FileHash $path -Algorithm SHA256).Hash.ToLower()

        if($h -eq $Expected[$packagedName]){
            $status="BYTE_IDENTICAL"
        }else{
            $status="BYTE_DIFFERENT_INFORMATIONAL_ONLY"
        }

        return "$label`t$status`t$h`t$path"

    }else{

        return "$label`tNOT_FOUND_INFORMATIONAL_ONLY`tNA`t$path"
    }
}

$lines=@()
$lines += "V2 STEP96J3 LOCAL HISTORICAL COPY COMPARISON"
$lines += ""
$lines += "These comparisons are informational only."
$lines += "The exact packaged source assets actually used for the final composition are independently hash-verified."
$lines += ""
$lines += Compare-Local "Figure3_ABC_PDF" $LocalABCpdf "_SOURCE_Figure3_ABC_EXACT.pdf"
$lines += Compare-Local "Figure3_ABC_PNG" $LocalABCpng "_SOURCE_Figure3_ABC_EXACT.png"
$lines += Compare-Local "Figure3_D_PDF" $LocalDpdf "_SOURCE_Figure3_D_EXACT.pdf"
$lines += Compare-Local "Figure3_D_PNG" $LocalDpng "_SOURCE_Figure3_D_EXACT.png"

$lines | Set-Content -Encoding UTF8 $Comparison

# ------------------------------------------------------------
# 4. Completion record
# ------------------------------------------------------------
@"
V2 STEP96J3 FINAL FIGURE 3 SELF-CONTAINED STAGING

Launcher resource root:
$PackageRoot

Packaged exact Figure 3 A-C source: hash verified.
Packaged exact Figure 3D source: hash verified.
Packaged exact Supplementary source: hash verified.
Final Figure 3 A-D asset: hash verified.

Old local historical-copy hashes are informational only.

No statistical analysis rerun.
No inferential content changed.

STEP96J3 COMPLETE
"@ | Set-Content -Encoding UTF8 (Join-Path $Out "_STEP96J3_COMPLETE.txt")

if(Test-Path $Zip){
    Remove-Item $Zip -Force
}

Compress-Archive `
  -Path "$Out\*" `
  -DestinationPath $Zip `
  -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP96J3 COMPLETE" -ForegroundColor Green
Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
