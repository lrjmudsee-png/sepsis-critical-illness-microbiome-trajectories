$ErrorActionPreference="Stop"

$Code="E:\sepsis_project\code\03_data_processing"
$Out="E:\sepsis_project\results\V2_96J_FINAL_FIGURE3_COMPOSITION"
$Zip="E:\sepsis_project\results\V2_96J_FINAL_FIGURE3_COMPOSITION_FOR_UPLOAD.zip"

$Base="E:\sepsis_project\results\V2_35D3_FINAL_FIGURE_FORMATTING_FIX\01_FINAL_FIGURES\Figure_3_Core_Sepsis_Heterogeneous_Taxonomic_Routes.pdf"
$PanelD="E:\sepsis_project\results\V2_96E_FIGURE3D_CROSS_COHORT_REPRODUCIBILITY_FINAL\Figure3D_Cross_Cohort_Taxonomic_Concordance_MAIN.pdf"
$Supp="E:\sepsis_project\results\V2_96E_FIGURE3D_CROSS_COHORT_REPRODUCIBILITY_FINAL\Supplementary_Figure_Cross_Cohort_Taxonomic_Concordance_Sensitivity.pdf"

$ExpectedBase="c48532d047db366a60583b3d6f4bd6b92e96658bacd9a827509eb6e97d0b08fa"
$ExpectedD="714bf4aba5ea2eae65677fba3c7d98dba5199f58e68c23aa1da360aeb1bba116"
$ExpectedSupp="034ef6b134d7a048fd33cf4d4ac5c6f687ccde1617242a046a6c1d2b43e7e6a5"

foreach($p in @($Base,$PanelD,$Supp)){
  if(!(Test-Path $p)){ throw "Missing source figure asset: $p" }
}

$ActualBase=(Get-FileHash $Base -Algorithm SHA256).Hash.ToLower()
$ActualD=(Get-FileHash $PanelD -Algorithm SHA256).Hash.ToLower()
$ActualSupp=(Get-FileHash $Supp -Algorithm SHA256).Hash.ToLower()

if($ActualBase -ne $ExpectedBase){ throw "Frozen Figure 3 A-C source hash mismatch." }
if($ActualD -ne $ExpectedD){ throw "Figure 3D source hash mismatch." }
if($ActualSupp -ne $ExpectedSupp){ throw "Supplementary sensitivity source hash mismatch." }

New-Item -ItemType Directory -Force -Path $Out | Out-Null

$Files=@(
  "V2_96J_Figure_3_FINAL_ABCD.pdf",
  "V2_96J_Figure_3_FINAL_ABCD_600dpi.png",
  "V2_96J_Figure_3_FINAL_LEGEND.txt",
  "V2_96J_Supplementary_Taxonomic_Concordance_Sensitivity.pdf",
  "V2_96J_Supplementary_Taxonomic_Concordance_Sensitivity.png",
  "V2_96J_Supplementary_Taxonomic_Concordance_LEGEND.txt",
  "V2_96J_SOURCE_ASSET_SHA256.csv"
)

foreach($f in $Files){
  $src=Join-Path $Code $f
  if(!(Test-Path $src)){ throw "Missing Step96J packaged asset: $src" }
  Copy-Item $src (Join-Path $Out $f) -Force
}

@"
V2 STEP96J FINAL FIGURE 3 COMPOSITION

Source Figure 3 A-C hash verified: TRUE
Source Figure 3D hash verified: TRUE
Supplementary sensitivity hash verified: TRUE

Final Figure 3 A-D staged from the exact verified source assets.
No statistical analysis rerun.
No panel inferential content changed.
Panel A-C retained from Step35D3.
Panel D retained from Step96E.
Only page composition and D title were added.

STEP96J COMPLETE
"@ | Set-Content -Encoding UTF8 (Join-Path $Out "_STEP96J_COMPLETE.txt")

if(Test-Path $Zip){ Remove-Item $Zip -Force }
Compress-Archive -Path "$Out\*" -DestinationPath $Zip -CompressionLevel Optimal

Write-Host ""
Write-Host "STEP96J COMPLETE" -ForegroundColor Green
Write-Host "Upload:"
Write-Host $Zip -ForegroundColor Cyan
