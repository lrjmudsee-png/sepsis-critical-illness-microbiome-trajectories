"""Assemble supplemental review pages from native R vector figure PDFs."""
from pathlib import Path
import argparse,html,io
from pypdf import PdfReader,PdfWriter,Transformation
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build_directory',type=Path);p.add_argument('output_pdf',type=Path)
    p.add_argument('--font-dir',type=Path,required=True)
    a=p.parse_args()
    if a.output_pdf.exists():p.error('Choose a new output PDF; existing files are not overwritten.')
    pdfmetrics.registerFont(TTFont('Arial',str(a.font_dir/'arial.ttf')))
    pdfmetrics.registerFont(TTFont('Arial-Bold',str(a.font_dir/'arialbd.ttf')))
    pdfmetrics.registerFontFamily('Arial',normal='Arial',bold='Arial-Bold')
    width,height=letter;left=43.2;inner=width-2*left
    body=ParagraphStyle('body',fontName='Arial',fontSize=11,leading=15.5)
    caption=ParagraphStyle('caption',fontName='Arial',fontSize=12,leading=17.5)
    out=PdfWriter();buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=letter)
    c.setFont('Arial-Bold',16);c.drawString(left,height-53,'Supplemental figures and data sets')
    entries=[
      ('','Longitudinal gut microbiome displacement across sepsis and critical illness cohorts despite heterogeneous taxonomic trajectories'),
      ('Contents','Fig. S1–S6 and a guide to the accompanying supplemental files. Detailed methods, processing parameters, equations and metadata qualifications are in Text S1.'),
      ('Data Set S1','Data_Set_S1.xlsx. Fourteen result sheets (S1–S14) plus a README: discovery evidence, trajectories, external comparisons, infection-source analyses, compositional robustness and independent-validation sensitivities. Sheet labels identify parts of this workbook.'),
      ('Data Set S2','Data_Set_S2.zip. Laboratory and processing provenance for eleven cohorts, including source/deposit platform evidence and unrecovered source information.'),
      ('Data Set S3','Data_Set_S3.zip. Fixed secondary family-balance and healthy-reference analyses, their definitions, exclusions and sensitivity results.'),
      ('Data Set S4','Data_Set_S4.zip. Post-hoc PRJNA516701 specimen-consistency audit and sensitivity results.'),
      ('Data Set S5','Data_Set_S5.zip. Selected flat frozen statistical inputs, a data dictionary, checksums and parameterized runners. These reproduce selected contrasts rather than full upstream ASV/OTU distance calculations.'),
      ('Data Set S6','Data_Set_S6.zip. Final figure-rendering scripts, display settings, all 24 frozen figure-input CSVs, checksums and reproduction instructions. Data and software are separated within the archive; the software license does not relicense input data.'),
      ('Data terms','Source data remain subject to their original terms. The MIT software license does not relicense third-party data. Historical filenames inside data archives retain their original provenance labels.')]
    y=height-77
    for heading,text in entries:
        paragraph=Paragraph(('<b>'+heading+'.</b> ' if heading else '')+html.escape(text),body)
        _,h=paragraph.wrap(inner,height);paragraph.drawOn(c,left,y-h);y-=h+12
    if y<38:raise RuntimeError('Cover text does not fit.')
    c.showPage();c.save();buf.seek(0);out.add_page(PdfReader(buf).pages[0])
    lines=[s.strip() for s in (a.build_directory/'FINAL_SUPPLEMENTARY_FIGURE_LEGENDS.txt').read_text(encoding='utf-8').splitlines() if s.strip()]
    figures=sorted((a.build_directory/'03_SUPPLEMENTARY_FIGURES').glob('*.pdf'))
    if len(figures)!=6 or len(lines)!=6:raise RuntimeError('Exactly six figures and legends are required.')
    for i,(file,legend) in enumerate(zip(figures,lines),1):
        legend=legend.replace('Supplementary Figure','Fig.').replace(' | ','. ',1)
        markup=html.escape(legend).replace('Hedges g_z','Hedges g<sub>z</sub>')
        paragraph=Paragraph(markup,caption);_,h=paragraph.wrap(inner,height)
        fig=PdfReader(file).pages[0]
        scale=min(inner/float(fig.mediabox.width),(height-150-h)/float(fig.mediabox.height))
        if scale<.86:raise RuntimeError('Legend leaves insufficient figure space; use an additional legend page.')
        fw,fh=float(fig.mediabox.width)*scale,float(fig.mediabox.height)*scale
        top=height-70;buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=letter)
        c.setFont('Arial-Bold',12);c.drawString(left,height-42,f'Fig. S{i}')
        paragraph.drawOn(c,left,top-fh-14-h)
        c.setFont('Arial',9);c.drawCentredString(width/2,24,f'Supplemental figures | {i+1}')
        c.showPage();c.save();buf.seek(0);page=PdfReader(buf).pages[0]
        page.merge_transformed_page(fig,Transformation().scale(scale).translate((width-fw)/2,top-fh));out.add_page(page)
    out.add_metadata({'/Title':'Supplemental figures S1–S6','/Author':'Lu et al.'})
    a.output_pdf.parent.mkdir(parents=True,exist_ok=True)
    with a.output_pdf.open('wb') as f:out.write(f)
    print(a.output_pdf)

if __name__=='__main__':main()
