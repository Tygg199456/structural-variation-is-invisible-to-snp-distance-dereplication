#!/usr/bin/env python3
"""Compose Supplementary.pdf from tables/S*.tsv.

One file, as the journal's checklist asks for.  Tables are typeset with reportlab so
no LaTeX is needed.  Run from the repository root:

    python3 scripts/build_supplementary.py
"""
import csv, glob, os

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (LongTable, PageBreak, Paragraph, SimpleDocTemplate,
                                Spacer, TableStyle)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "Supplementary.pdf")

TITLES = {
 "S1_software": ("Table S1.", "Complete software version manifest, including the reason for each pinned dependency."),
 "S2_recurrent_blocks": ("Table S2a.", "The 161 distinct >=10 kb uncalled coordinates and how many genomes share each (5 kb-binned start)."),
 "S2_accessions": ("Table S2b.", "All genome assemblies analysed (accessions), in the order used by the pipeline."),
 "S3_sweep_tierA": ("Table S3a.", "Threshold sweep on Tier A: retained genomes, largest cluster diameter and size, and same-strain pairs merged."),
 "S3_sweep_tierB": ("Table S3b.", "Threshold sweep on Tier B, including the two off-grid integers (t = 18, 19) that bracket the arithmetic floor."),
 "S4_plasmid_three_state": ("Table S4.", "Plasmid three-state reclassification of the 258 closed genomes from each assembly's own submitted sequence inventory."),
 "S5_baseline_layer": ("Table S5.", "Per-genome classification of the baseline layer, with the S1/S2 decompositions."),
 "S6_pseudogene_entanglement": ("Table S6.", "Pseudogene-deletion entanglement in the curated CO92 annotation against the automated annotation."),
 "S7_interpretation": ("Table S7.", "Interpretation-layer guidance for the four settings of Section 4.7."),
 "S8_second_reference": ("Table S8.", "The dereplication analysis repeated against a second reference for the same 256 genomes. Coordinates are not comparable between the two references; read each column within itself."),
}

styles = getSampleStyleSheet()
H = ParagraphStyle("H", parent=styles["Heading2"], spaceAfter=4)
CAP = ParagraphStyle("CAP", parent=styles["BodyText"], fontSize=8.5, leading=11, spaceAfter=6, alignment=TA_LEFT)
NOTE = ParagraphStyle("NOTE", parent=styles["BodyText"], fontSize=8.5, leading=12, spaceAfter=8)
CELL = ParagraphStyle("CELL", parent=styles["BodyText"], fontSize=6.2, leading=7.4)
CELLH = ParagraphStyle("CELLH", parent=CELL, fontName="Helvetica-Bold")


def read_tsv(path):
    src, header, rows = "", None, []
    for line in open(path, newline=""):
        line = line.rstrip("\n")
        if line.startswith("# source:"):
            src = line.split(":", 1)[1].strip()
            continue
        if not line:
            continue
        parts = line.split("\t")
        if header is None:
            header = parts
        else:
            rows.append(parts)
    return src, header, rows


def build():
    doc = SimpleDocTemplate(OUT, pagesize=A4,
                            leftMargin=15 * mm, rightMargin=15 * mm,
                            topMargin=15 * mm, bottomMargin=15 * mm,
                            title="Supplementary material", author="")
    story = [Paragraph("Supplementary material", styles["Title"]),
             Spacer(1, 6),
             Paragraph("Structural variation is invisible to SNP-distance dereplication: the "
                       "reference strain of <i>Yersinia pestis</i> merges with its own attenuated "
                       "derivative.", styles["Heading3"]),
             Spacer(1, 8),
             Paragraph("All tables below are generated from the archived derived tables by "
                       "<font face='Courier'>scripts/make_tables.py</font> and typeset by "
                       "<font face='Courier'>scripts/build_supplementary.py</font>. Each table "
                       "carries the file it comes from. Supplementary Note 1 and Supplementary "
                       "Note 2 are part of the manuscript text.", NOTE),
             PageBreak()]
    for path in sorted(glob.glob(os.path.join(ROOT, "tables", "S*.tsv"))):
        stem = os.path.basename(path)[:-4]
        label, cap = TITLES.get(stem, (stem, ""))
        src, header, rows = read_tsv(path)
        story += [Paragraph("%s %s" % (label, cap), CAP)]
        ncol = 1
        if len(header) == 1 and len(rows) > 200:
            ncol = 4                                    # long single-column lists: pack them
            per = (len(rows) + ncol - 1) // ncol
            cols = [rows[i * per:(i + 1) * per] for i in range(ncol)]
            rows = [[str(c[i]) if i < len(c) else "" for c in cols] for i in range(per)]
            header = [header[0]] * ncol
        data = [[Paragraph(h.replace("_", " "), CELLH) for h in header]]
        for r in rows:
            data.append([Paragraph(str(c), CELL) for c in r])
        avail = doc.width
        if len(header) == 1:
            widths = [avail]
        else:
            first = min(avail * 0.45, 0.45 * avail)
            rest = (avail - first) / (len(header) - 1)
            widths = [first] + [rest] * (len(header) - 1)
        t = LongTable(data, colWidths=widths, repeatRows=1)
        t.setStyle(TableStyle([
            ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#bbbbbb")),
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#eeeeee")),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 2),
            ("RIGHTPADDING", (0, 0), (-1, -1), 2),
            ("TOPPADDING", (0, 0), (-1, -1), 1.2),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 1.2),
        ]))
        story += [t, Paragraph("source: %s" % src, ParagraphStyle("src", parent=CELL, fontSize=6,
                                                                  textColor=colors.HexColor("#666666"))),
                  Spacer(1, 10)]
    doc.build(story)
    print("wrote %s (%.1f MB)" % (os.path.relpath(OUT, ROOT), os.path.getsize(OUT) / 1e6))


if __name__ == "__main__":
    build()
