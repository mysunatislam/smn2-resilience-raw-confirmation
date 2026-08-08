from __future__ import annotations

import argparse
import csv
import hashlib
import zipfile
from datetime import date
from pathlib import Path

from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "manuscript"
SUBMISSION = MANUSCRIPT / "submission"
PDF_OUTPUT = ROOT / "output" / "pdf" / "supplementary_material_review.pdf"
TABLE_OUTPUT = ROOT / "output" / "submission" / "supplementary_tables_S1-S15.zip"
MANIFEST_OUTPUT = ROOT / "output" / "submission" / "submission_artifact_manifest.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def metadata_is_complete() -> bool:
    rows = read_tsv(SUBMISSION / "metadata_required.tsv")
    required = [row for row in rows if not row["status"].startswith("OPTIONAL")]
    return all(row["status"] == "PASS" for row in required)


def figure_number(figure_id: str) -> int:
    return int(figure_id.removeprefix("S"))


def draw_footer(pdf: canvas.Canvas, page_width: float, page_number: int) -> None:
    pdf.setStrokeColor(colors.HexColor("#C7CDD4"))
    pdf.line(18 * mm, 14 * mm, page_width - 18 * mm, 14 * mm)
    pdf.setFillColor(colors.HexColor("#4E5964"))
    pdf.setFont("Helvetica", 8)
    pdf.drawString(18 * mm, 9 * mm, "Supplementary material - review build")
    pdf.drawRightString(page_width - 18 * mm, 9 * mm, str(page_number))


def draw_paragraph(
    pdf: canvas.Canvas,
    text: str,
    style: ParagraphStyle,
    x: float,
    y_top: float,
    width: float,
) -> float:
    paragraph = Paragraph(text, style)
    _, height = paragraph.wrap(width, 1000 * mm)
    paragraph.drawOn(pdf, x, y_top - height)
    return height


def build_pdf(final_mode: bool) -> None:
    if final_mode and not metadata_is_complete():
        raise RuntimeError(
            "Final mode requires completed author and declaration metadata."
        )

    PDF_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figures = read_tsv(MANUSCRIPT / "supplementary_figure_index.tsv")
    captions = {
        row["figure_id"]: row["caption"]
        for row in read_tsv(SUBMISSION / "supplementary_figure_captions.tsv")
    }
    tables = read_tsv(MANUSCRIPT / "supplementary_table_index.tsv")
    if len(figures) != 25 or len(captions) != 25 or len(tables) != 15:
        raise RuntimeError("Expected 25 figures, 25 captions, and 15 tables.")

    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        "FigureTitle",
        parent=styles["Heading2"],
        fontName="Helvetica-Bold",
        fontSize=11,
        leading=14,
        textColor=colors.HexColor("#18212A"),
        alignment=TA_LEFT,
        spaceAfter=4,
    )
    caption_style = ParagraphStyle(
        "Caption",
        parent=styles["BodyText"],
        fontName="Helvetica",
        fontSize=9,
        leading=12,
        textColor=colors.HexColor("#26313B"),
        alignment=TA_LEFT,
    )
    cover_style = ParagraphStyle(
        "Cover",
        parent=styles["Title"],
        fontName="Helvetica-Bold",
        fontSize=19,
        leading=23,
        alignment=TA_CENTER,
        textColor=colors.HexColor("#17212B"),
    )

    pdf = canvas.Canvas(str(PDF_OUTPUT), pagesize=A4, pageCompression=1)
    pdf.setTitle("Supplementary material - motor-neuron resilience in SMA")
    page_number = 1
    width, height = A4
    pdf.setFillColor(colors.HexColor("#A52A2A"))
    pdf.rect(0, height - 18 * mm, width, 18 * mm, fill=1, stroke=0)
    pdf.setFillColor(colors.white)
    pdf.setFont("Helvetica-Bold", 10)
    pdf.drawCentredString(width / 2, height - 11.5 * mm, "REVIEW BUILD - NOT FOR JOURNAL UPLOAD")
    y = height - 50 * mm
    y -= draw_paragraph(
        pdf,
        "Cross-model analysis prioritizes reproducible human motor-neuron resilience candidates in spinal muscular atrophy",
        cover_style,
        25 * mm,
        y,
        width - 50 * mm,
    )
    pdf.setFillColor(colors.HexColor("#26313B"))
    pdf.setFont("Helvetica", 11)
    pdf.drawCentredString(width / 2, y - 15 * mm, "Supplementary material")
    pdf.drawCentredString(width / 2, y - 23 * mm, "Target journal: Journal of Neurology")
    pdf.setFont("Helvetica", 9)
    pdf.drawCentredString(width / 2, y - 34 * mm, f"Generated {date.today().isoformat()}")
    note = (
        "This review artifact contains the frozen computational results. It is not a final submission file: "
        "author metadata, declarations, and corresponding-author details remain unresolved, and no functional "
        "or RT-PCR results are represented as completed evidence."
    )
    draw_paragraph(pdf, note, caption_style, 28 * mm, y - 52 * mm, width - 56 * mm)
    draw_footer(pdf, width, page_number)
    pdf.showPage()

    for row in sorted(figures, key=lambda item: figure_number(item["figure_id"])):
        image_path = MANUSCRIPT / row["filename"]
        if not image_path.exists():
            raise FileNotFoundError(image_path)
        with PILImage.open(image_path) as image:
            pixel_width, pixel_height = image.size
        page_size = landscape(A4) if pixel_width / pixel_height > 1.15 else A4
        pdf.setPageSize(page_size)
        width, height = page_size
        page_number += 1
        margin_x = 18 * mm
        top = height - 18 * mm
        title = f"Supplementary Figure {row['figure_id']}. {row['title']}."
        title_height = draw_paragraph(
            pdf, title, title_style, margin_x, top, width - 2 * margin_x
        )
        caption_top = top - title_height - 2 * mm
        caption_height = draw_paragraph(
            pdf,
            captions[row["figure_id"]],
            caption_style,
            margin_x,
            caption_top,
            width - 2 * margin_x,
        )
        image_top = caption_top - caption_height - 5 * mm
        image_bottom = 20 * mm
        available_width = width - 2 * margin_x
        available_height = image_top - image_bottom
        scale = min(
            available_width / pixel_width,
            available_height / pixel_height,
        )
        draw_width = pixel_width * scale
        draw_height = pixel_height * scale
        x = (width - draw_width) / 2
        y = image_bottom + (available_height - draw_height) / 2
        pdf.drawImage(
            str(image_path),
            x,
            y,
            width=draw_width,
            height=draw_height,
            preserveAspectRatio=True,
            mask="auto",
        )
        draw_footer(pdf, width, page_number)
        pdf.showPage()

    pdf.setPageSize(A4)
    width, height = A4
    page_number += 1
    y = height - 20 * mm
    y -= draw_paragraph(
        pdf,
        "Supplementary Tables S1-S15",
        ParagraphStyle(
            "TableSection",
            parent=title_style,
            fontSize=15,
            leading=18,
        ),
        20 * mm,
        y,
        width - 40 * mm,
    )
    y -= 5 * mm
    y -= draw_paragraph(
        pdf,
        "The complete tables are supplied separately as machine-readable TSV files in supplementary_tables_S1-S15.zip.",
        caption_style,
        20 * mm,
        y,
        width - 40 * mm,
    )
    y -= 6 * mm
    for row in tables:
        entry = f"<b>Supplementary Table {row['table_id']}.</b> {row['title']}."
        paragraph = Paragraph(entry, caption_style)
        _, item_height = paragraph.wrap(width - 40 * mm, 1000 * mm)
        if y - item_height < 22 * mm:
            draw_footer(pdf, width, page_number)
            pdf.showPage()
            page_number += 1
            y = height - 20 * mm
        paragraph.drawOn(pdf, 20 * mm, y - item_height)
        y -= item_height + 4 * mm
    draw_footer(pdf, width, page_number)
    pdf.save()


def build_table_archive() -> None:
    TABLE_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    tables = read_tsv(MANUSCRIPT / "supplementary_table_index.tsv")
    with zipfile.ZipFile(TABLE_OUTPUT, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for row in tables:
            source = MANUSCRIPT / row["filename"]
            if not source.exists():
                raise FileNotFoundError(source)
            archive.write(source, arcname=source.name)
        archive.write(
            MANUSCRIPT / "supplementary_table_index.tsv",
            arcname="supplementary_table_index.tsv",
        )
        archive.writestr(
            "README.txt",
            "Supplementary Tables S1-S15 for the frozen SMA motor-neuron resilience analysis.\n"
            "Files are tab-delimited UTF-8 text and should be cited by supplementary table number.\n",
        )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_manifest() -> None:
    MANIFEST_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    files = [PDF_OUTPUT, TABLE_OUTPUT]
    with MANIFEST_OUTPUT.open("w", encoding="ascii", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["artifact", "bytes", "sha256", "submission_status"])
        for path in files:
            writer.writerow(
                [
                    path.relative_to(ROOT).as_posix(),
                    path.stat().st_size,
                    sha256(path),
                    "REVIEW_ONLY",
                ]
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("review", "final"), default="review")
    args = parser.parse_args()
    build_pdf(final_mode=args.mode == "final")
    build_table_archive()
    write_manifest()
    print(PDF_OUTPUT)
    print(TABLE_OUTPUT)
    print(MANIFEST_OUTPUT)


if __name__ == "__main__":
    main()
