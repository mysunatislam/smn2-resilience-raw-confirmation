from __future__ import annotations

import csv
from pathlib import Path

from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
FIGURE_ROOT = ROOT / "manuscript" / "figures"
SUBMISSION_ROOT = ROOT / "manuscript" / "submission"
AUDIT_PATH = SUBMISSION_ROOT / "figure_asset_audit.tsv"
MARKER_PATH = SUBMISSION_ROOT / "FIGURE_ASSET_AUDIT_COMPLETE.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def png_metrics(path: Path) -> tuple[int, int, float, float]:
    with Image.open(path) as image:
        width, height = image.size
        dpi = image.info.get("dpi", (0.0, 0.0))
    return width, height, float(dpi[0]), float(dpi[1])


def main() -> None:
    SUBMISSION_ROOT.mkdir(parents=True, exist_ok=True)
    assets = [
        ("Figure 1", FIGURE_ROOT / "figure_1_raw_vs_processed_concordance.png"),
        ("Figure 2", FIGURE_ROOT / "figure_2_candidate_evidence_matrix.png"),
        ("Figure 3", FIGURE_ROOT / "figure_3_raw_splice_confirmation.png"),
    ]
    for row in read_tsv(ROOT / "manuscript" / "supplementary_figure_index.tsv"):
        assets.append((row["figure_id"], ROOT / "manuscript" / row["filename"]))

    records: list[dict[str, object]] = []
    for asset_id, png_path in assets:
        if not png_path.exists():
            records.append(
                {
                    "asset_id": asset_id,
                    "png_path": png_path.relative_to(ROOT).as_posix(),
                    "pixel_width": 0,
                    "pixel_height": 0,
                    "dpi_x": 0,
                    "dpi_y": 0,
                    "vector_pdf": False,
                    "publication_basis": "missing",
                    "status": "FAIL",
                }
            )
            continue
        width, height, dpi_x, dpi_y = png_metrics(png_path)
        pdf_path = png_path.with_suffix(".pdf")
        vector_pdf = pdf_path.exists() and pdf_path.stat().st_size > 1000
        high_resolution_png = min(dpi_x, dpi_y) >= 590 and min(width, height) >= 2400
        if vector_pdf:
            basis = "vector_pdf"
        elif high_resolution_png:
            basis = "png_600_dpi"
        else:
            basis = "insufficient_raster_without_pdf"
        records.append(
            {
                "asset_id": asset_id,
                "png_path": png_path.relative_to(ROOT).as_posix(),
                "pixel_width": width,
                "pixel_height": height,
                "dpi_x": f"{dpi_x:.1f}",
                "dpi_y": f"{dpi_y:.1f}",
                "vector_pdf": vector_pdf,
                "publication_basis": basis,
                "status": "PASS" if vector_pdf or high_resolution_png else "FAIL",
            }
        )

    fields = list(records[0])
    with AUDIT_PATH.open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)

    failures = [row for row in records if row["status"] != "PASS"]
    with MARKER_PATH.open("w", encoding="ascii", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["status", "assets_checked", "assets_passed", "assets_failed"])
        writer.writerow(
            ["COMPLETE" if not failures else "FAILED", len(records), len(records) - len(failures), len(failures)]
        )
    if failures:
        names = ", ".join(str(row["asset_id"]) for row in failures)
        raise SystemExit(f"Figure asset audit failed: {names}")
    print(f"Figure asset audit passed for {len(records)} assets")


if __name__ == "__main__":
    main()
