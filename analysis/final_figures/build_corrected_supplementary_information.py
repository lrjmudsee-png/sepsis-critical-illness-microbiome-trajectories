from __future__ import annotations

from pathlib import Path
import argparse

from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Image, PageBreak, Paragraph, SimpleDocTemplate, Spacer


parser = argparse.ArgumentParser(description="Assemble supplementary figures without rerunning inference.")
parser.add_argument("build_directory", type=Path)
parser.add_argument("output_pdf", type=Path)
parser.add_argument("--font-dir", type=Path, required=True, help="Directory containing licensed arial.ttf and arialbd.ttf")
args = parser.parse_args()
BUILD = args.build_directory.resolve(strict=True)
FIGDIR = BUILD / "03_SUPPLEMENTARY_FIGURES"
OUTPUT = args.output_pdf.resolve()
if OUTPUT.exists():
    parser.error("Output PDF already exists; choose a new output path.")
OUTPUT.parent.mkdir(parents=True, exist_ok=True)

font_regular = args.font_dir / "arial.ttf"
font_bold = args.font_dir / "arialbd.ttf"
pdfmetrics.registerFont(TTFont("Arial", str(font_regular)))
pdfmetrics.registerFont(TTFont("Arial-Bold", str(font_bold)))

styles = getSampleStyleSheet()
title = ParagraphStyle(
    "Title",
    parent=styles["Title"],
    fontName="Arial-Bold",
    fontSize=20,
    leading=24,
    alignment=TA_CENTER,
    spaceAfter=12,
)
subtitle = ParagraphStyle(
    "Subtitle",
    parent=styles["Normal"],
    fontName="Arial",
    fontSize=12,
    leading=16,
    alignment=TA_CENTER,
    spaceAfter=18,
)
h2 = ParagraphStyle(
    "H2",
    parent=styles["Heading2"],
    fontName="Arial-Bold",
    fontSize=11.5,
    leading=14,
    textColor=colors.HexColor("#17365D"),
    spaceBefore=7,
    spaceAfter=3,
)
body = ParagraphStyle(
    "Body",
    parent=styles["BodyText"],
    fontName="Arial",
    fontSize=9.2,
    leading=12.3,
    alignment=TA_LEFT,
    spaceAfter=5,
)
caption = ParagraphStyle(
    "Caption",
    parent=styles["BodyText"],
    fontName="Arial",
    fontSize=8.3,
    leading=10.5,
    alignment=TA_LEFT,
    spaceBefore=7,
)


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Arial", 8)
    canvas.setFillColor(colors.HexColor("#666666"))
    canvas.drawCentredString(A4[0] / 2, 0.42 * inch, f"Supplementary Information | {doc.page}")
    canvas.restoreState()


def fit_image(path: Path, max_w: float, max_h: float) -> Image:
    with PILImage.open(path) as im:
        w, h = im.size
    scale = min(max_w / w, max_h / h)
    return Image(str(path), width=w * scale, height=h * scale)


legend_lines = [
    line.strip()
    for line in (BUILD / "FINAL_SUPPLEMENTARY_FIGURE_LEGENDS.txt").read_text(encoding="utf-8").splitlines()
    if line.strip()
]
if len(legend_lines) != 6:
    raise RuntimeError(f"Expected 6 supplementary legends, found {len(legend_lines)}")

figure_paths = [
    FIGDIR / "Supplementary_Figure_S1_Taxonomic_Concordance_Sensitivity.png",
    FIGDIR / "Supplementary_Figure_S2_Pair_Counts.png",
    FIGDIR / "Supplementary_Figure_S3_Bray_Meta_Analysis.png",
    FIGDIR / "Supplementary_Figure_S4_CZM_Aitchison_Meta_Analysis.png",
    FIGDIR / "Supplementary_Figure_S5_Family_Balance.png",
    FIGDIR / "Supplementary_Figure_S6_Healthy_Reference.png",
]

story = [
    Spacer(1, 0.45 * inch),
    Paragraph("Supplementary Information", title),
    Paragraph(
        "Longitudinal gut microbiome displacement across sepsis and critical illness cohorts despite heterogeneous taxonomic trajectories",
        subtitle,
    ),
    Paragraph(
        "This file contains Supplementary Figures S1-S6 and descriptions of Supplementary Data 1-5, which accompany this submission. Detailed methods appear in the main manuscript, Methods M1-M7.",
        body,
    ),
    Paragraph("Supplementary Data 1", h2),
    Paragraph(
        "<b>Supplementary_Data_1.xlsx.</b> Fourteen frozen result sheets: S1 longitudinal evidence; S2 core taxonomic trajectories; S3 external cross-sectional ecological-state comparison; S4 infection-source ecology; S5 source-model sensitivity; S6 source-specific slopes; S7 baseline beta diversity; S8 baseline genus associations; S9 longitudinal genus effects; S10 genus-branch definitions; S11 genus reporting; S12 geometry robustness; S13 independent longitudinal validation; S14 independent longitudinal validation governance sensitivities. Historical sheet labels S1-S14 identify sheets within this single data file, not separate PDF tables.",
        body,
    ),
    Paragraph("Supplementary Data 2", h2),
    Paragraph(
        "<b>Supplementary_Data_2.zip.</b> Eleven-cohort laboratory and processing provenance, and representative source/deposit platform evidence. Missing source facts and the PRJNA1166732 platform-label conflict are retained.",
        body,
    ),
    Paragraph("Supplementary Data 3", h2),
    Paragraph(
        "<b>Supplementary_Data_3.zip.</b> Fixed secondary direction-balance and healthy-reference results, including null findings and analysis definitions. These analyses do not establish mechanism or prognosis.",
        body,
    ),
    Paragraph("Supplementary Data 4", h2),
    Paragraph(
        "<b>Supplementary_Data_4.zip.</b> Post-hoc PRJNA516701 specimen-consistency evidence and sensitivity results. The original primary analysis is unchanged.",
        body,
    ),
    Paragraph("Supplementary Data 5", h2),
    Paragraph(
        "<b>Supplementary_Data_5.zip.</b> Ten selected flat frozen-input tables, data dictionary, input hashes and two parameterized R runners. This reproduces selected statistics only, not full ASV/OTU objects or upstream distance calculations.",
        body,
    ),
    Spacer(1, 6),
    Paragraph(
        "Source data remain subject to their original terms; the MIT software license does not relicense datasets.",
        body,
    ),
]

for i, (figure_path, legend) in enumerate(zip(figure_paths, legend_lines), start=1):
    story.append(PageBreak())
    prefix, rest = legend.split(" | ", 1)
    story.append(Paragraph(prefix, h2))
    max_h = 7.6 * inch
    story.append(fit_image(figure_path, 7.0 * inch, max_h))
    rest = rest.replace("Hedges g_z", "Hedges g<sub>z</sub>")
    story.append(Paragraph(f"<b>{prefix}.</b> {rest}", caption))

doc = SimpleDocTemplate(
    str(OUTPUT),
    pagesize=A4,
    leftMargin=0.62 * inch,
    rightMargin=0.62 * inch,
    topMargin=0.58 * inch,
    bottomMargin=0.62 * inch,
    title="Supplementary Information",
    author="Lu Rongji et al.",
)
doc.build(story, onFirstPage=footer, onLaterPages=footer)
print(OUTPUT)
