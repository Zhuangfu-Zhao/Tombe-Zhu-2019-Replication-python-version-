from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.lib import colors
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
import os

cwd = r"E:\mine\project\structural_models\Replication-to-TombeZhu\Code\20150811_backup"
output_path = os.path.join(cwd, "raw_data_description.pdf")

doc = SimpleDocTemplate(output_path, pagesize=letter,
                        leftMargin=0.75*inch, rightMargin=0.75*inch,
                        topMargin=0.75*inch, bottomMargin=0.75*inch)

styles = getSampleStyleSheet()
title_style = ParagraphStyle('Title', parent=styles['Title'], fontSize=18, spaceAfter=24)
section_style = ParagraphStyle('Section', parent=styles['Heading1'], fontSize=14, spaceBefore=20, spaceAfter=12)
subsection_style = ParagraphStyle('Subsection', parent=styles['Heading2'], fontSize=11, spaceBefore=10, spaceAfter=6)
text_style = ParagraphStyle('Text', parent=styles['Normal'], fontSize=10, spaceAfter=4)
tabletext_style = ParagraphStyle('TableText', parent=styles['Normal'], fontSize=9, spaceAfter=2)

story = []

story.append(Paragraph("Raw Data Description", title_style))
story.append(Spacer(1, 8))
story.append(Paragraph("Trade, Migration and Productivity: A Quantitative Analysis of China", text_style))
story.append(Paragraph("Replication Files (20150811_backup)", text_style))
story.append(Spacer(1, 20))

sections = [
    ("1. migration_data.dta", 
        "Full province-sector matrix of migration shares, bilateral distances, and instrumental variables",
        "3,540 rows x 14 columns",
        ["Variable", "Description"],
        [("importer", "Province name (destination)"),
        ("exporter", "Province name (origin)"),
        ("j", "Sector code (1=agriculture, 2=non-agriculture)"),
        ("i", "Sector code (1=agriculture, 2=non-agriculture)"),
        ("mij_mii2000", "Migration share relative to staying in 2000"),
        ("mij_mii2005", "Migration share relative to staying in 2005"),
        ("Vj", "Real income in destination province-sector in 2000"),
        ("Vi", "Real income in origin province-sector in 2000"),
        ("mij2000", "Bilateral migration share in 2000"),
        ("mij2005", "Bilateral migration share in 2005"),
        ("distance", "Bilateral distance between provinces (km)"),
        ("rYLj_bartik_sector", "Bartik IV: weighted average of neighbor income by sector"),
        ("rYLj_bartik_sector_uhs2000", "Bartik IV: unweighted harmonic mean specification"),
        ("rYLj_IVneighbour", "Neighbor income as IV for migration elasticity")]
    ),
    ("2. trade_data.dta",
        "Bilateral trade flows between 8 regions for 2002 and 2007",
        "144 rows x 25 columns",
        ["Variable", "Description"],
        [("i", "Origin region code"),
        ("j", "Destination region code"),
        ("importer", "Region name (destination)"),
        ("exporter", "Region name (origin)"),
        ("year", "Year (2002 or 2007)"),
        ("Eij_ag", "Exports from i to j in agriculture"),
        ("Eij_na", "Exports from i to j in non-agriculture"),
        ("Ei_ag", "Total exports from i in agriculture"),
        ("Ei_na", "Total exports from i in non-agriculture"),
        ("Eii_ag", "Home consumption in agriculture"),
        ("Eii_na", "Home consumption in non-agriculture"),
        ("PIij_ag", "Trade share agriculture"),
        ("PIij_na", "Trade share non-agriculture"),
        ("PIii_ag", "Home share agriculture"),
        ("PIii_na", "Home share non-agriculture"),
        ("ti_ag", "Exporter-specific cost agriculture"),
        ("ti_na", "Exporter-specific cost non-agriculture"),
        ("tj_ag", "Importer-specific cost agriculture"),
        ("tj_na", "Importer-specific cost non-agriculture"),
        ("tij_ag", "Bilateral trade cost agriculture"),
        ("tij_na", "Bilateral trade cost non-agriculture")]
    ),
    ("3. employment_realGDP_data.csv",
        "Provincial real GDP per worker and employment in 2000 and 2005",
        "30 rows x 9 columns (30 provinces)",
        ["Variable", "Description"],
        [("province", "Province name"),
        ("rYLa2000", "Real GDP per worker in agriculture, 2000"),
        ("rYLn2000", "Real GDP per worker in non-agriculture, 2000"),
        ("La2000", "Employment in agriculture, 2000 (millions)"),
        ("Ln2000", "Employment in non-agriculture, 2000 (millions)"),
        ("rYLa2005", "Real GDP per worker in agriculture, 2005"),
        ("rYLn2005", "Real GDP per worker in non-agriculture, 2005"),
        ("drYLa_data", "Observed change in agriculture real income"),
        ("drYLn_data", "Observed change in non-agriculture real income")]
    ),
    ("4. mij2000.csv & mij2005.csv",
        "Bilateral migration shares in vector form (to be reshaped in MATLAB)",
        "3,600 rows x 2 columns each",
        ["Variable", "Description"],
        [("mij_mii", "Migration share relative to staying (mij/mii)"),
        ("mij", "Raw migration share (bilateral)"),
        ("", "Matrix structure: 60x60 (30 provinces x 2 sectors)"),
        ("", "Flattened: [prov1_ag, prov1_na, prov2_ag, ..., prov30_na]")]
    ),
    ("5. trade_ag.csv & trade_na.csv",
        "Bilateral trade matrix between 30 provinces in 2002",
        "30 rows x 30 columns each",
        ["", ""],
        [("Matrix Structure", "30x30 bilateral trade flow matrix"),
        ("Row", "Exporter province"),
        ("Column", "Importer province"),
        ("Value", "Trade share (proportion)"),
        ("Note", "Provinces ordered by 'region' variable")]
    ),
    ("6. tauhat.csv",
        "Measured change in migration costs",
        "961 rows x 6 columns",
        ["Variable", "Description"],
        [("dni_ag", "Change in migration cost (agriculture, symmetric)"),
        ("dni_na", "Change in migration cost (non-agriculture, symmetric)"),
        ("dni_asym_ag", "Change in asymmetric cost (agriculture)"),
        ("dni_asym_na", "Change in asymmetric cost (non-agriculture)"),
        ("dni_dist_ag", "Change in distance cost component (agriculture)"),
        ("dni_dist_na", "Change in distance cost component (non-agriculture)")]
    ),
    ("7. dT.mat",
        "Estimated productivity changes for model to match observed real income changes",
        "MATLAB binary file",
        ["", ""],
        [("Content", "Productivity changes for each province-sector"),
        ("Dimensions", "60x1 vector (30 provinces x 2 sectors)")]
    )
]

for title, desc, shape, header, *rows in sections:
    story.append(Paragraph(title, section_style))
    story.append(Paragraph(desc, text_style))
    story.append(Paragraph(f"Shape: {shape}", text_style))
    story.append(Spacer(1, 6))
    
    table_data = [header]
    if title in ["4. mij2000.csv & mij2005.csv", "5. trade_ag.csv & trade_na.csv", "7. dT.mat"]:
        for r in rows[0]:
            table_data.append(list(r))
    else:
        for r in rows[0]:
            table_data.append(list(r))
    
    col_widths = [2*inch, 4*inch]
    table = Table(table_data, colWidths=col_widths, repeatRows=1)
    table.setStyle(TableStyle([
        ('FONTSIZE', (0, 0), (-1, -1), 9),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
        ('TEXTTOP', (0, 0), (-1, -1), 2),
        ('LEFTPADDING', (0, 0), (-1, -1), 4),
        ('RIGHTPADDING', (0, 0), (-1, -1), 4),
        ('GRID', (0, 0), (-1, -1), 0.5, colors.black),
        ('BACKGROUND', (0, 0), (-1, 0), colors.lightgrey),
        ('BACKGROUND', (0, 1), (-1, -1), colors.white),
    ]))
    story.append(table)
    story.append(Spacer(1, 12))

story.append(PageBreak())

story.append(Paragraph("Additional Information", section_style))
notes_data = [
    ["Sector Codes", "1 = Agriculture, 2 = Non-Agriculture"],
    ["Years", "Migration: 2000, 2005 | Trade: 2002, 2007"],
    ["Regions (8)", "Central coast, Eastern coast, Far east, North, North coast, Northwest, Southcentral, Southwest"],
    ["Provinces (30)", "Anhui, Beijing, Chongqing, Fujian, Gansu, Guangdong, Guangxi, Guizhou, Hainan, Hebei, Heilongjiang, Henan, Hubei, Hunan, Inner Mongolia, Jiangsu, Jiangxi, Jilin, Liaoning, Ningxia, Qinghai, Shandong, Shanghai, Shannxi, Shanxi, Sichuan, Tianjin, Xinjiang, Yunnan, Zhejiang"],
]

table = Table(notes_data, colWidths=[1.5*inch, 5*inch])
table.setStyle(TableStyle([
    ('FONTSIZE', (0, 0), (-1, -1), 9),
    ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
    ('TEXTTOP', (0, 0), (-1, -1), 2),
    ('LEFTPADDING', (0, 0), (-1, -1), 4),
    ('RIGHTPADDING', (0, 0), (-1, -1), 4),
    ('GRID', (0, 0), (-1, -1), 0.5, colors.black),
]))
story.append(table)

story.append(Spacer(1, 20))
story.append(Paragraph("Note: mij2000.csv and mij2005.csv are vector forms that must be reshaped to 60x60 matrices in MATLAB.", text_style))
story.append(Paragraph("Data Source: China National Bureau of Statistics", text_style))
story.append(Paragraph("Citation: Tombe, T. and Zhu, X. (2015) 'Trade, Migration and Productivity: A Quantitative Analysis of China'", text_style))

doc.build(story)
print(f"PDF created: {output_path}")